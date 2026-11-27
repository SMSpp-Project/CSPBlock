/*--------------------------------------------------------------------------*/
/*---------------------- File FlatZincCSPSolver.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the FlatZincCSPSolver class, the Solver of a CSPBlock
 * through an external FlatZinc solver run as a separate process.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- IMPLEMENTATION -----------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <thread>

#include <cctype>
#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <cstring>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "FlatZincCSPSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register FlatZincCSPSolver in the Solver factory

SMSpp_insert_in_factory_cpp_0( FlatZincCSPSolver );

/*--------------------------------------------------------------------------*/
/*---------------------------- LOCAL FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

/// true if \p file is found and executable, on the PATH if it has no slash

bool executable( const std::string & file )
{
 if( file.find( '/' ) != std::string::npos )
  return( access( file.c_str() , X_OK ) == 0 );

 const char * path = std::getenv( "PATH" );
 if( ! path )
  return( false );
 std::istringstream in( path );
 std::string dir;
 while( std::getline( in , dir , ':' ) ) {
  if( dir.empty() )
   dir = ".";
  if( access( ( dir + "/" + file ).c_str() , X_OK ) == 0 )
   return( true );
  }
 return( false );
 }

/*--------------------------------------------------------------------------*/
/// the template of a temporary file \p name in the temporary directory of
/// the system, empty if there is none

std::string temp_template( const char * name )
{
 std::error_code ec;
 const auto dir = std::filesystem::temp_directory_path( ec );
 return( ec ? std::string() : ( dir / name ).string() );
 }

}  // namespace

/*--------------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/

void FlatZincCSPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 CSPBlock * csp = nullptr;
 if( block && ! ( csp = dynamic_cast< CSPBlock * >( block ) ) )
  throw( std::invalid_argument( "FlatZincCSPSolver::set_Block: the Block is "
				"not a CSPBlock" ) );

 Solver::set_Block( block );
 f_csp = csp;
 f_reload = true;
 f_status = kUnEval;
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR HANDLING PARAMETERS ---------------------*/
/*--------------------------------------------------------------------------*/

void FlatZincCSPSolver::set_par( idx_type par , double value )
{
 if( par == dblMaxTime ) {
  MaxTime = value;
  f_reload = true;  // the cached solution was found under another time limit
  }
 else
  Solver::set_par( par , value );
 }

/*--------------------------------------------------------------------------*/

void FlatZincCSPSolver::set_par( idx_type par , std::string && value )
{
 if( par == strFlatZincSolver ) {
  FlatZincSolver = std::move( value );
  f_reload = true;  // another backend may find another solution
  }
 else
  Solver::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & FlatZincCSPSolver::get_dflt_str_par( idx_type par ) const
{
 static const std::string empty;
 if( par == strFlatZincSolver )
  return( empty );
 return( Solver::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & FlatZincCSPSolver::get_str_par( idx_type par ) const
{
 if( par == strFlatZincSolver )
  return( FlatZincSolver );
 return( Solver::get_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type FlatZincCSPSolver::str_par_str2idx( const std::string & name )
 const
{
 if( name == "strFlatZincSolver" )
  return( strFlatZincSolver );
 return( Solver::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & FlatZincCSPSolver::str_par_idx2str( idx_type idx ) const
{
 static const std::string name = "strFlatZincSolver";
 if( idx == strFlatZincSolver )
  return( name );
 return( Solver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

void FlatZincCSPSolver::set_par( idx_type par ,
				 std::vector< std::string > && value )
{
 if( par == vstrFlatZincArgs ) {
  FlatZincArgs = std::move( value );
  f_reload = true;  // other arguments may find another solution
  }
 else
  Solver::set_par( par , std::move( value ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< std::string > &
 FlatZincCSPSolver::get_dflt_vstr_par( idx_type par ) const
{
 static const std::vector< std::string > empty;
 if( par == vstrFlatZincArgs )
  return( empty );
 return( Solver::get_dflt_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

const std::vector< std::string > &
 FlatZincCSPSolver::get_vstr_par( idx_type par ) const
{
 if( par == vstrFlatZincArgs )
  return( FlatZincArgs );
 return( Solver::get_vstr_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type
 FlatZincCSPSolver::vstr_par_str2idx( const std::string & name ) const
{
 if( name == "vstrFlatZincArgs" )
  return( vstrFlatZincArgs );
 return( Solver::vstr_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & FlatZincCSPSolver::vstr_par_idx2str( idx_type idx ) const
{
 static const std::string name = "vstrFlatZincArgs";
 if( idx == vstrFlatZincArgs )
  return( name );
 return( Solver::vstr_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int FlatZincCSPSolver::compute( bool changedvars )
{
 if( ! f_csp )
  throw( std::logic_error( "FlatZincCSPSolver::compute: no CSPBlock "
			   "attached" ) );

 // the backend is given by the configuration, with no default in the code:
 // an empty strFlatZincSolver is a misconfiguration, not a reason to invent
 // one
 if( FlatZincSolver.empty() )
  throw( std::invalid_argument( "FlatZincCSPSolver::compute: no backend; "
				"strFlatZincSolver is empty, set it in the "
				"ComputeConfig" ) );

 lock();

 // read the Modification: with nothing changed since the last compute() and
 // a result already there, the cached one stands [see the class]
 process_outstanding_Modification();
 if( ( ! f_reload ) && ( f_status != kUnEval ) ) {
  unlock();
  return( f_status );
  }

 f_status = kUnEval;
 f_has_sol = false;
 f_complete = false;
 f_value = Inf< double >();
 v_sol.clear();

 if( ! backend_found() ) {  // the backend is set but not found on the PATH
  f_status = kError;
  unlock();
  return( f_status );
  }

 // write the FlatZinc of the CSPBlock to a temporary file
 std::string tmpl = temp_template( "fzncsp_XXXXXX.fzn" );
 const int fd = tmpl.empty() ? -1 : mkstemps( tmpl.data() , 4 );
 if( fd < 0 ) {
  f_status = kError;
  unlock();
  return( f_status );
  }
 const std::string model = flatzinc();
 const bool written =
  write( fd , model.data() , model.size() ) == ssize_t( model.size() );
 close( fd );
 const std::string fzn = tmpl;

 if( ! written ) {
  unlink( fzn.c_str() );
  f_status = kError;
  unlock();
  return( f_status );
  }

 std::string out;
 const int rc = run_backend( fzn , out );
 unlink( fzn.c_str() );

 if( rc < 0 ) {  // the backend could not be run
  f_status = kError;
  unlock();
  return( f_status );
  }

 f_status = parse_output( out );

 // the backend enforces a soft constraint as a hard one, its weight being an
 // annotation, and ignores the costs of a "solve satisfy": unsatisfiability
 // then proves nothing, and a solution is not proven optimal
 if( ( f_status == kInfeasible ) && ( ! f_csp->all_hard() ) )
  f_status = kError;
 if( f_status == kOK ) {
  f_complete = f_complete && f_csp->all_hard() &&
               ( ( ! f_csp->has_costs() ) ||
                 ( model.find( "\nsolve satisfy;" ) == std::string::npos ) );
  if( ! f_complete )
   f_status = kLowPrecision;
  }

 if( rc == 1 ) {  // terminated because the time limit was reached
  f_complete = false;
  if( f_status != kInfeasible )
   f_status = kStopTime;
  }

 f_reload = false;  // this result holds until the next Modification
 unlock();
 return( f_status );
 }

/*--------------------------------------------------------------------------*/

void FlatZincCSPSolver::process_outstanding_Modification( void )
{
 // the backend is run on the textual FlatZinc of the whole CSPBlock,
 // rewritten at each compute(): there is no incremental interface, so any
 // Modification (a change of the physical representation, a CSPBlockChange of
 // a branching) means the model is written and solved again
 while( pop() )
  f_reload = true;
 }

/*--------------------------------------------------------------------------*/

std::string FlatZincCSPSolver::flatzinc( void ) const
{
 std::ostringstream full;
 f_csp->print( full , 'C' );

 // mark every variable as an output_var, so that the backend prints its
 // value: print() writes one "var ...: x<i> ... ;" per line, and adding
 // " :: output_var" before the ";" leaves the other annotations in place;
 // an inactive constraint, which the backend would not know to be so from
 // its annotation, is left out
 static const std::string inactive = " :: inactive;";
 std::ostringstream out;
 std::istringstream in( full.str() );
 std::string line;
 while( std::getline( in , line ) ) {
  if( ( line.rfind( "constraint " , 0 ) == 0 ) &&
      ( line.size() >= inactive.size() ) &&
      ( line.compare( line.size() - inactive.size() , inactive.size() ,
		      inactive ) == 0 ) )
   continue;
  if( ( line.rfind( "var " , 0 ) == 0 ) && ( line.back() == ';' ) )
   line.insert( line.size() - 1 , " :: output_var" );
  out << line << "\n";
  }
 return( out.str() );
 }

/*--------------------------------------------------------------------------*/

int FlatZincCSPSolver::run_backend( const std::string & fzn ,
				    std::string & out ) const
{
 if( FlatZincSolver.empty() )
  return( -1 );

 // the command: the executable, its arguments, and the FlatZinc file last
 std::vector< std::string > args;
 args.reserve( FlatZincArgs.size() + 2 );
 args.push_back( FlatZincSolver );
 for( const auto & a : FlatZincArgs )
  args.push_back( a );
 args.push_back( fzn );

 std::vector< char * > argv;
 argv.reserve( args.size() + 1 );
 for( auto & a : args )
  argv.push_back( const_cast< char * >( a.c_str() ) );
 argv.push_back( nullptr );

 // the standard output of the backend goes to a temporary file
 std::string otmpl = temp_template( "fzncsp_out_XXXXXX" );
 const int ofd = otmpl.empty() ? -1 : mkstemp( otmpl.data() );
 if( ofd < 0 )
  return( -1 );

 const pid_t pid = fork();
 if( pid < 0 ) {
  close( ofd );
  unlink( otmpl.c_str() );
  return( -1 );
  }

 if( pid == 0 ) {  // the child runs the backend
  dup2( ofd , STDOUT_FILENO );
  close( ofd );
  setpgid( 0 , 0 );  // own process group, to terminate the whole tree
  execvp( argv[ 0 ] , argv.data() );
  _exit( 127 );  // execvp failed
  }

 setpgid( pid , pid );  // also here, whichever call wins the race

 const bool finite = MaxTime < Inf< double >();
 const auto start = std::chrono::steady_clock::now();
 int wstatus = 0;
 int terminated = 0;
 for( ; ; ) {
  const pid_t r = waitpid( pid , & wstatus , WNOHANG );
  if( r == pid )  // the backend ended by itself
   break;
  if( ( r < 0 ) && ( errno != EINTR ) )  // waitpid failed
   break;
  if( finite ) {
   const std::chrono::duration< double > el =
    std::chrono::steady_clock::now() - start;
   if( el.count() >= MaxTime ) {  // the time limit is reached
    kill( - pid , SIGTERM );
    std::this_thread::sleep_for( std::chrono::milliseconds( 200 ) );
    kill( - pid , SIGKILL );
    waitpid( pid , & wstatus , 0 );
    terminated = 1;
    break;
    }
   }
  std::this_thread::sleep_for( std::chrono::milliseconds( 10 ) );
  }

 // read the output of the backend back
 lseek( ofd , 0 , SEEK_SET );
 std::ifstream read( otmpl );
 std::ostringstream buf;
 buf << read.rdbuf();
 out = buf.str();
 close( ofd );
 unlink( otmpl.c_str() );

 if( ( ! terminated ) && WIFEXITED( wstatus ) &&
     ( WEXITSTATUS( wstatus ) == 127 ) )
  return( -1 );  // execvp failed in the child

 return( terminated );
 }

/*--------------------------------------------------------------------------*/

int FlatZincCSPSolver::parse_output( const std::string & out )
{
 const auto n = f_csp->get_number_variables();
 std::vector< int > latest;      // the last complete assignment
 std::vector< int > current( n , 0 );
 std::vector< unsigned char > set( n , 0 );
 bool infeasible = false;
 bool optimal = false;

 std::istringstream in( out );
 std::string line;
 while( std::getline( in , line ) ) {
  if( line.find( "UNSATISFIABLE" ) != std::string::npos ) {
   infeasible = true;
   continue;
   }
  if( ( line.find( "UNKNOWN" ) != std::string::npos ) ||
      ( line.find( "UNBOUNDED" ) != std::string::npos ) ||
      ( line.rfind( "=====ERROR" , 0 ) != std::string::npos ) )
   continue;
  if( line.rfind( "----------" , 0 ) == 0 ) {  // end of a solution
   bool whole = true;
   for( std::size_t i = 0 ; whole && ( i < n ) ; ++i )
    whole = set[ i ];
   if( whole )
    latest = current;
   std::fill( set.begin() , set.end() , 0 );
   continue;
   }
  if( line.rfind( "==========" , 0 ) == 0 ) {  // search complete
   optimal = true;
   continue;
   }

  // an assignment "x<i> = <value> ;"
  const auto eq = line.find( '=' );
  if( eq == std::string::npos )
   continue;
  std::string name = line.substr( 0 , eq );
  while( ( ! name.empty() ) && std::isspace( ( unsigned char ) name.back() ) )
   name.pop_back();
  if( ( name.size() < 2 ) || ( name[ 0 ] != 'x' ) )
   continue;
  std::size_t idx = 0;
  try { idx = std::stoul( name.substr( 1 ) ); }
  catch( ... ) { continue; }
  if( idx >= n )
   continue;
  std::string val = line.substr( eq + 1 );
  const auto semi = val.find( ';' );
  if( semi != std::string::npos )
   val = val.substr( 0 , semi );
  try { current[ idx ] = std::stoi( val ); }
  catch( ... ) { continue; }
  set[ idx ] = 1;
  }

 if( infeasible && latest.empty() )
  return( kInfeasible );
 if( latest.empty() )
  return( kError );

 v_sol = std::move( latest );
 f_has_sol = true;
 try { f_value = f_csp->get_objective_value( v_sol ); }
 catch( ... ) {  // a value outside its universe: the solution is not usable
  f_has_sol = false;
  v_sol.clear();
  f_value = Inf< double >();
  return( kError );
  }
 // a solution of a problem with no costs and only hard constraints is
 // optimal, whether the backend proves the search complete or not
 f_complete = optimal || ( ( ! f_csp->has_costs() ) && f_csp->all_hard() );
 return( kOK );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

void FlatZincCSPSolver::get_var_solution( Configuration * )
{
 if( ! f_has_sol )
  throw( std::logic_error( "FlatZincCSPSolver::get_var_solution: no "
			   "solution" ) );

 // the ColVariable are there only if the abstract representation has been
 // generated, and then there is one per variable
 auto & var = f_csp->get_abstract_variables();
 if( var.size() != v_sol.size() )
  return;
 for( std::size_t i = 0 ; i < v_sol.size() ; ++i )
  var[ i ].set_value( v_sol[ i ] );
 }

/*--------------------------------------------------------------------------*/

bool FlatZincCSPSolver::backend_found( void ) const
{
 return( ( ! FlatZincSolver.empty() ) && executable( FlatZincSolver ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------- End File FlatZincCSPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/

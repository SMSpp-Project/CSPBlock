/*--------------------------------------------------------------------------*/
/*------------------- File FlatZincCSPSolver_test.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test for FlatZincCSPSolver on a 4 queens CSPBlock: the BlockSolverConfig
 * FlatZincCSPSolverBSCfg.txt attaches the Solver, compute() runs the backend
 * on the FlatZinc of the CSPBlock, and the solution it reads back is checked
 * for feasibility [see CSPBlock::is_feasible()]. If the backend of the
 * configuration is not installed the test exits with 77, which CTest reads as
 * "skipped", so that the module builds and tests without a FlatZinc solver.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <sstream>
#include <string>

#include "BlockSolverConfig.h"
#include "CSPBlock.h"
#include "FlatZincCSPSolver.h"

// the checks hold in every build type, Release included
#undef NDEBUG
#include <cassert>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// the 4 queens in FlatZinc, as MiniZinc would flatten it [see battery.cpp]

static std::string queens4( void )
{
 std::ostringstream o;
 for( int i = 0 ; i < 4 ; ++i )
  o << "var 1..4: q" << i << " :: output_var;\n";
 o << "array [1..4] of var int: q = [q0, q1, q2, q3];\n";
 o << "constraint fzn_all_different_int(q);\n";
 for( int i = 0 ; i < 4 ; ++i )
  for( int j = i + 1 ; j < 4 ; ++j ) {
   o << "constraint int_lin_ne([1, -1], [q" << i << ", q" << j << "], "
     << j - i << ");\n";
   o << "constraint int_lin_ne([1, -1], [q" << i << ", q" << j << "], "
     << i - j << ");\n";
   }
 o << "solve satisfy;\n";
 return( o.str() );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 CSPBlock b;
 std::istringstream in( queens4() );
 b.load( in );

 auto cfg = Configuration::deserialize( "FlatZincCSPSolverBSCfg.txt" );
 auto bsc = dynamic_cast< BlockSolverConfig * >( cfg );
 assert( bsc );
 bsc->apply( & b );

 assert( ! b.get_registered_solvers().empty() );
 auto solver = dynamic_cast< FlatZincCSPSolver * >(
			       b.get_registered_solvers().front() );
 assert( solver );

 // the configuration bites: a malformed ComputeConfig is discarded whole and
 // in silence, so I check that the parameters are those of the file and not
 // the defaults (an empty backend, no arguments, no time limit)
 assert( solver->get_str_par( FlatZincCSPSolver::strFlatZincSolver ) ==
	 "minizinc" );
 const auto & args =
  solver->get_vstr_par( FlatZincCSPSolver::vstrFlatZincArgs );
 assert( ( args.size() == 2 ) && ( args[ 0 ] == "--solver" ) &&
	 ( args[ 1 ] == "cp-sat" ) );
 assert( solver->get_dbl_par( Solver::dblMaxTime ) == 60 );

 int rv = 0;
 if( ! solver->backend_found() ) {
  std::cout << "FlatZincCSPSolver: the backend \""
	    << solver->get_str_par( FlatZincCSPSolver::strFlatZincSolver )
	    << "\" is not installed, skipping" << std::endl;
  rv = 77;  // CTest reads 77 as "skipped"
  }
 else {
  const int status = solver->compute();
  const bool ok = ( status == Solver::kOK ) && solver->has_var_solution() &&
		  b.is_feasible( solver->get_solution() );
  std::cout << "FlatZincCSPSolver: 4 queens, status " << status
	    << ( ok ? ": OK" : ": FAILED" ) << std::endl;
  assert( ok );

  // reoptimization: restrict the first variable to its value in the
  // solution, which issues a CSPBlockMod the Solver reads, and solve again
  const CSPBlock::v_int sol( solver->get_solution().begin() ,
			     solver->get_solution().end() );
  b.chg_domains( Block::Subset{ 0 } , CSPBlock::v_v_int{ { sol[ 0 ] } } );
  const int status2 = solver->compute();
  const bool ok2 = ( status2 == Solver::kOK ) &&
		   solver->has_var_solution() &&
		   b.is_feasible( solver->get_solution() ) &&
		   ( solver->get_solution()[ 0 ] == sol[ 0 ] );
  std::cout << "FlatZincCSPSolver: reoptimization, status " << status2
	    << ( ok2 ? ": OK" : ": FAILED" ) << std::endl;
  assert( ok2 );
  }

 bsc->clear();
 bsc->apply( & b );
 delete cfg;

 // the backend of the configuration bites on the behaviour: one that is not
 // installed makes compute() return kError (a valid one gives kOK), so the
 // parameter is not silently ignored [see FlatZincCSPSolverBSCfg-missing.txt]
 auto cfgm =
  Configuration::deserialize( "FlatZincCSPSolverBSCfg-missing.txt" );
 auto bscm = dynamic_cast< BlockSolverConfig * >( cfgm );
 assert( bscm );
 bscm->apply( & b );
 auto solverm = dynamic_cast< FlatZincCSPSolver * >(
				b.get_registered_solvers().front() );
 assert( solverm );
 assert( solverm->get_str_par( FlatZincCSPSolver::strFlatZincSolver ) ==
	 "fzn-no-such-backend" );
 assert( solverm->compute() == Solver::kError );
 bscm->clear();
 bscm->apply( & b );
 delete cfgm;

 return( rv );
 }

/*--------------------------------------------------------------------------*/
/*----------------- End File FlatZincCSPSolver_test.cpp --------------------*/
/*--------------------------------------------------------------------------*/

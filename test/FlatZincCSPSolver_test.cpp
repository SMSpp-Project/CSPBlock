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
  }

 bsc->clear();
 bsc->apply( & b );
 delete cfg;
 return( rv );
 }

/*--------------------------------------------------------------------------*/
/*----------------- End File FlatZincCSPSolver_test.cpp --------------------*/
/*--------------------------------------------------------------------------*/

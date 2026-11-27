/*--------------------------------------------------------------------------*/
/*----------------------- File FlatZincCSPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class FlatZincCSPSolver, which implements the
 * Solver concept [see Solver.h] for a CSPBlock [see CSPBlock.h] by writing its
 * FlatZinc and running the MiniZinc driver on it as a separate process,
 * reading the solution back from the output of the process. The executable
 * and its arguments are parameters, so that the same class drives OR-Tools
 * CP-SAT, Gecode, Chuffed or any other solver that MiniZinc knows.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __FlatZincCSPSolver
 #define __FlatZincCSPSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <string>
#include <vector>

#include "CSPBlock.h"
#include "Solver.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup FlatZincCSPSolver_CLASSES Classes in FlatZincCSPSolver.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS FlatZincCSPSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/// a Solver of a CSPBlock through an external FlatZinc solver
/** The FlatZincCSPSolver class implements the Solver concept [see Solver.h]
 * for a CSPBlock [see CSPBlock.h] by delegating to an external FlatZinc
 * solver, run as a separate process: compute() writes the CSPBlock in
 * FlatZinc [see CSPBlock::print()], preceded by the few definitions that make
 * it a MiniZinc model (see below), to a temporary file, runs the MiniZinc
 * driver on it, and reads the values of the variables back from its output.
 *
 * The backend is not fixed: the string parameter strFlatZincSolver is its
 * executable, i.e., "minizinc" (or a path to it), looked for on the PATH,
 * and the vector-of-string parameter vstrFlatZincArgs is its arguments, each
 * an entry (e.g., [ "--solver" , "com.google.ortools.sat" ] for the CP-SAT
 * of OR-Tools, "org.gecode.gecode" for Gecode, "org.chuffed.chuffed" for
 * Chuffed, the identifiers that do not change from a version of MiniZinc to
 * another); the command run is the executable, the arguments, and the path
 * of the model last. The executable has no
 * default in the code, and it is given by a configuration file: compute()
 * throws if strFlatZincSolver is empty, and returns kError if it is not
 * found. The arguments are an entry each, rather than one string, because a
 * ComputeConfig reads a string parameter as a single token, while a
 * vector-of-string one is a count and then its entries.
 *
 * compute() writes every variable of the CSPBlock as an output_var, so that
 * the backend prints its value, and leaves the inactive constraints out; it
 * reads from the output the last complete assignment (for "solve
 * minimize" the backend prints the improving ones, the last being the best),
 * sets the status to kOK if a solution is found and proven optimal,
 * kLowPrecision if a solution is found but the backend could not prove it
 * optimal (a soft constraint, or costs that are not a linear objective in
 * one variable [see CSPBlock::print()]), kStopTime if the time limit
 * dblMaxTime is reached first (the process is then terminated), kInfeasible
 * if the backend reports =====UNSATISFIABLE===== and all the constraints are
 * hard, and kError if it fails for any other reason or cannot be run. After
 * kOK, kLowPrecision and kStopTime get_solution() is the assignment found,
 * if any, one value per variable, which is_feasible() of the CSPBlock checks
 * and get_objective_value() evaluates; get_ub() is that value and get_lb()
 * is it too after kOK, -INF otherwise.
 *
 * The solution is read through get_solution(), one value per variable, as
 * CSPBlock::is_feasible( std::span ) takes it; get_var_solution() also
 * writes it into the ColVariable of the abstract representation of the
 * CSPBlock, if this has been generated [see
 * CSPBlock::generate_abstract_variables()].
 *
 * A FlatZincCSPSolver reads the Modification of the CSPBlock [see
 * add_Modification()], the CSPBlockMod of a change of the physical
 * representation and those a CSPBlockChange of a branching issues: the backend
 * is run on the FlatZinc of the whole CSPBlock, rewritten at each compute(),
 * so there is no incremental interface to it, and any Modification only marks
 * the cached solution stale, the next compute() rewriting the FlatZinc and
 * running the backend again; a compute() with no Modification since the last
 * one returns the cached result. A warm start of the backend is not possible
 * through the textual FlatZinc, so the reoptimization is "rewrite the
 * changed model and solve it again".
 *
 * The model is the FlatZinc of print() read as MiniZinc, which the driver
 * flattens for the backend with the library of the backend, so that each
 * global constraint is native or decomposed as the backend wants: it
 * includes globals.mzn, for all_different and cumulative, and defines the
 * positive table of FlatZinc (one array of tuples, while MiniZinc wants a 2D
 * one), the negative table and the element of an offset other than 1, which
 * print() writes as predicates that are not in FlatZinc. The annotations of
 * print() that MiniZinc does not know (universe, group, cost, weight) are
 * removed. The weight of a soft constraint and the costs of the values, but
 * for a linear objective in one variable, are thus lost: the backend
 * enforces a soft constraint as a hard one and looks for any feasible
 * solution. The hard constraint satisfaction problems, with at most a linear
 * objective in one variable, are what it solves exactly. */

class FlatZincCSPSolver : public Solver
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor: it does nothing, the backend is run by compute()

 FlatZincCSPSolver( void ) : Solver() {}

/*------------------------------- DESTRUCTOR -------------------------------*/
 /// destructor, does nothing

 ~FlatZincCSPSolver() override = default;

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// sets the CSPBlock to solve; throws if the Block is not a CSPBlock

 void set_Block( Block * block ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets the double parameters; dblMaxTime is the time limit of compute()

 void set_par( idx_type par , double value ) override;

 using Solver::set_par;  // keep the other set_par() visible

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the double parameters, dblMaxTime included

 [[nodiscard]] double get_dbl_par( idx_type par ) const override {
  return( par == dblMaxTime ? MaxTime : Solver::get_dbl_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the string parameters of FlatZincCSPSolver, on top of those of Solver

 enum str_par_type_FZCSP {
  strFlatZincSolver = strLastAlgPar ,  ///< the executable of the backend
                               /**< The executable of the FlatZinc backend
				* (i.e., "minizinc" or a path), looked for
				* on the PATH; empty (the default) makes
				* compute() throw. */
  strLastAlgParFZCSP           ///< first new string parameter of derived
                               ///< classes
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets the string parameters, strFlatZincSolver included

 void set_par( idx_type par , std::string && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of string parameters

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return( idx_type( strLastAlgParFZCSP ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the default of the string parameters

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the string parameters

 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of the string parameter with the given name

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the name of the string parameter with the given index

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the vector-of-string parameters, on top of those of Solver

 enum vstr_par_type_FZCSP {
  vstrFlatZincArgs = vstrLastAlgPar ,  ///< the arguments of the backend
                               /**< The arguments passed to the executable of
				* strFlatZincSolver, one per entry (e.g.,
				* [ "--solver" , "org.gecode.gecode" ]); empty by
				* default. */
  vstrLastAlgParFZCSP          ///< first new vector-of-string parameter of
                               ///< derived classes
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets the vector-of-string parameters, vstrFlatZincArgs included

 void set_par( idx_type par , std::vector< std::string > && value ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of vector-of-string parameters

 [[nodiscard]] idx_type get_num_vstr_par( void ) const override {
  return( idx_type( vstrLastAlgParFZCSP ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the default of the vector-of-string parameters

 [[nodiscard]] const std::vector< std::string > &
  get_dflt_vstr_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the vector-of-string parameters

 [[nodiscard]] const std::vector< std::string > &
  get_vstr_par( idx_type par ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of the vector-of-string parameter with the given name

 [[nodiscard]] idx_type vstr_par_str2idx( const std::string & name )
  const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the name of the vector-of-string parameter with the given index

 [[nodiscard]] const std::string & vstr_par_idx2str( idx_type idx )
  const override;

/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
 /// solves the CSPBlock by running the backend on its FlatZinc
 /** Writes the CSPBlock in FlatZinc, runs the backend on it and reads the
  * solution back; the returned status is described in the comments to the
  * class. */

 int compute( bool changedvars = true ) override;

/*---------------------- METHODS FOR READING RESULTS -----------------------*/
 /// the lower bound: the value of the solution if complete, -INF otherwise

 [[nodiscard]] OFValue get_lb( void ) override {
  return( f_complete && f_has_sol ? OFValue( f_value ) : -Inf< OFValue >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the value of the solution found, +INF if none

 [[nodiscard]] OFValue get_ub( void ) override {
  return( f_has_sol ? OFValue( f_value ) : Inf< OFValue >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// true if the last compute() has found a solution

 [[nodiscard]] bool has_var_solution( void ) override { return( f_has_sol ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// writes the solution into the abstract ColVariable of the CSPBlock
 /** Writes the solution found into the integer ColVariable of the abstract
  * representation of the CSPBlock, if it has been generated, and does
  * nothing otherwise, the solution being read by get_solution(); throws if
  * there is no solution [see has_var_solution()]. */

 void get_var_solution( Configuration * solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the solution found, one value per variable, empty if none

 [[nodiscard]] const std::vector< int > & get_solution( void ) const {
  return( v_sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// true if the backend of strFlatZincSolver is found on the PATH
 /** Returns true if strFlatZincSolver is not empty and is found and
  * executable, on the PATH if it has no slash and as a path otherwise; a
  * test checks it to skip when no backend is installed. */

 [[nodiscard]] bool backend_found( void ) const;

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// reads the Modification of the queue, marking the cached solution stale
 /** Drains the Modification queue [see add_Modification()]: since the backend
  * has no incremental interface [see the class], any Modification sets
  * f_reload, so that the next compute() rewrites the FlatZinc and runs the
  * backend again. */

 void process_outstanding_Modification( void );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the FlatZinc of the CSPBlock, every variable an output_var, the
 /// inactive constraints left out

 [[nodiscard]] std::string flatzinc( void ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// runs the backend on the file \p fzn, its output into \p out
 /** Runs the command of strFlatZincSolver on the FlatZinc file \p fzn, the
  * standard output of the process collected into \p out; returns 0 if the
  * process ended by itself, 1 if it was terminated because the time limit
  * dblMaxTime was reached, and -1 if it could not be run. */

 int run_backend( const std::string & fzn , std::string & out ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// reads the FlatZinc output \p out into the solution, returning the status
 /** Reads the values of the variables from the FlatZinc output \p out, the
  * last complete assignment, into v_sol, setting f_has_sol, f_complete and
  * f_value; returns kOK if a complete solution is there, kInfeasible if \p
  * out reports unsatisfiability, and kError otherwise [see the class]. */

 int parse_output( const std::string & out );

/*---------------------------- PROTECTED FIELDS ----------------------------*/

 CSPBlock * f_csp = nullptr;   ///< the CSPBlock being solved

 std::string FlatZincSolver;   ///< the parameter strFlatZincSolver

 std::vector< std::string > FlatZincArgs;  ///< the parameter vstrFlatZincArgs

 double MaxTime = Inf< double >();  ///< the time limit of compute() (s)

 bool f_reload = true;         ///< if the model has to be rewritten and solved

 int f_status = kUnEval;       ///< status of the last compute()

 std::vector< int > v_sol;     ///< the solution, one value per variable

 bool f_has_sol = false;       ///< if the last compute() found a solution

 bool f_complete = false;      ///< if the backend ended the search

 double f_value = Inf< double >();  ///< the objective value of the solution

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*---------------------------- PRIVATE METHODS -----------------------------*/

 SMSpp_insert_in_factory_h;  // insert FlatZincCSPSolver in the Solver factory

/*--------------------------------------------------------------------------*/

 };  // end( class( FlatZincCSPSolver ) )

/** @} end( group( FlatZincCSPSolver_CLASSES ) ) */

/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* FlatZincCSPSolver.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File FlatZincCSPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/

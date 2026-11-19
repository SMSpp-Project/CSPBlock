/*--------------------------------------------------------------------------*/
/*----------------------- File FlatZincCSPSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class FlatZincCSPSolver, which implements the
 * Solver concept [see Solver.h] for a CSPBlock [see CSPBlock.h] by writing its
 * FlatZinc and running an external FlatZinc solver on it as a separate
 * process, reading the solution back from the output of the process. The
 * executable of the backend, and its fixed arguments, are a string parameter,
 * so that the same class drives OR-Tools CP-SAT (minizinc --solver cp-sat, or
 * fzn-cp-sat), Gecode, Chuffed or any other FlatZinc solver.
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
 * FlatZinc [see CSPBlock::print()] to a temporary file, runs the backend on
 * it, and reads the values of the variables back from its output, in the same
 * way as the :MILPSolver hand a MILP to CPLEX, Gurobi, SCIP or HiGHS.
 *
 * The backend is not fixed: the string parameter strFlatZincSolver is its
 * executable (e.g., "minizinc", "fzn-cp-sat", "fzn-gecode"), looked for on
 * the PATH, and the vector-of-string parameter vstrFlatZincArgs is its
 * arguments, each an entry (e.g., [ "--solver" , "cp-sat" ] for the CP-SAT of
 * OR-Tools through MiniZinc); the command run is the executable, the
 * arguments, and the path of the FlatZinc file last. The executable has no
 * default in the code, which is why a default lives in a configuration file
 * [see the test]: an empty strFlatZincSolver makes compute() return kError.
 * The arguments are an entry each, rather than one string, because a
 * ComputeConfig reads a string parameter as a single token, while a
 * vector-of-string one is a count and then its entries.
 *
 * compute() writes every variable of the CSPBlock as an output_var, so that
 * the backend prints its value; it reads from the FlatZinc output the last
 * complete assignment (for "solve minimize" the backend prints the improving
 * ones, the last being the best), sets the status to kOK if a solution is
 * found and the search is complete, kStopTime if the time limit dblMaxTime is
 * reached first (the process is then terminated), kInfeasible if the backend
 * reports =====UNSATISFIABLE=====, and kError if it fails for any other
 * reason or cannot be run. After kOK and kStopTime get_solution() is the
 * assignment found, one value per variable, which is_feasible() of the
 * CSPBlock checks and get_objective_value() evaluates; get_ub() is that value
 * and get_lb() is it too if the search is complete, -INF otherwise.
 *
 * Since the CSPBlock has no abstract representation [see CSPBlock.h], there
 * are no Variable to write the solution into: get_var_solution() has nothing
 * to do, and the solution is read through get_solution(), a std::span of the
 * values, as CSPBlock::is_feasible( std::span ) takes it.
 *
 * A FlatZincCSPSolver reads the Modification of the CSPBlock [see
 * add_Modification()], the CSPBlockMod of a change of the physical
 * representation and those a CSPBlockChange of a branching issues: the backend
 * is run on the FlatZinc of the whole CSPBlock, rewritten at each compute(),
 * so there is no incremental interface to it, and any Modification only marks
 * the cached solution stale, the next compute() rewriting the FlatZinc and
 * running the backend again; a compute() with no Modification since the last
 * one returns the cached result. A warm start of the backend is not possible
 * through the textual FlatZinc (minizinc in batch), so the reoptimization here
 * is "rewrite the changed model and solve it again", not a warm start of the
 * backend; a warm start is left to the decomposition Solver to come.
 *
 * The writer emits standard FlatZinc predicates for the positive table, the
 * linear constraints, all_different, the element of offset 1 and cumulative,
 * which a stock backend knows; the negative table, the element of another
 * offset and the soft (weighted) constraints are written as the predicates
 * and the annotations that CSPBlock::load() reads, which a stock backend does
 * not define. A CSPBlock using those therefore needs a FlatZinc library that
 * defines them to be solved this way, which is out of the scope of this
 * class; the hard constraint satisfaction problems of the standard predicates,
 * with at most a linear objective in one variable, are what it solves. */

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
				* (e.g., "minizinc", "fzn-cp-sat"), looked for
				* on the PATH; empty (the default) makes
				* compute() return kError. */
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
				* [ "--solver" , "cp-sat" ]); empty by
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
 /// the CSPBlock has no abstract representation, so this does nothing
 /** The CSPBlock has no Variable to write the solution into [see the class]:
  * this throws if there is no solution [see has_var_solution()] and does
  * nothing otherwise, the solution being read by get_solution(). */

 void get_var_solution( Configuration * solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the solution found, one value per variable, empty if none

 [[nodiscard]] const std::vector< int > & get_solution( void ) const {
  return( v_sol );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// true if the backend of strFlatZincSolver is found on the PATH
 /** Returns true if strFlatZincSolver is not empty and its first token, the
  * executable of the backend, is found and executable on the PATH or as a
  * path; what a test checks to skip when no backend is installed. */

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
 /// the FlatZinc of the CSPBlock, every variable an output_var

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

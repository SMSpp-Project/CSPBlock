/*--------------------------------------------------------------------------*/
/*---------------------------- File CSPBlock.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class CSPBlock, which implements the
 * Block concept [see Block.h] for the constraint satisfaction problems of
 * constraint programming (CP), i.e., finding values of a set of variables,
 * each with a finite set of integer values, satisfying a set of constraints
 * of a catalogue, and for their weighted version, where some constraints
 * are soft and the values have costs. Also the Modification and the Change
 * classes of the CSPBlock are here.
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

#ifndef __CSPBlock
 #define __CSPBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <span>
#include <vector>

#include "Block.h"
#include "Change.h"
#include "CSPGlobalConstraints.h"
#include "FRowConstraint.h"
#include "Objective.h"
#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup CSPBlock_CLASSES Classes in CSPBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS CSPBlock -----------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for constraint programming (CP)
/** The CSPBlock class implements the Block concept [see Block.h] for the
 * constraint satisfaction problems of constraint programming (CP): given n
 * variables v_0, ..., v_{n - 1}, each with a finite set of integer values,
 * and m constraints of the catalogue below, find a value of each variable
 * satisfying all the constraints.
 *
 * Each variable v has a *universe* U_v, the sorted set of the values it may
 * ever take, and a *current domain* D_v, the subset of U_v it may take now:
 * keeping both makes every later restriction of a domain (e.g., by a
 * branching) a change of D_v, with the indices of the values unchanged. The
 * values of all the universes are stored one variable after the other, so
 * that each value of each variable has a *slot*, its position there [see
 * value_index()], to which the current domain and the costs are aligned.
 *
 * Each constraint has a kind [see constraint_kind], a *scope*, i.e., the
 * indices of the variables it is on (a variable may appear more than once),
 * and integer *parameters*, whose meaning depends on the kind:
 *
 * - kExtension, the (positive) table: the parameters are tuples of as many
 *   values as the scope, one after the other, and the values of the scope
 *   must be one of them; the scope cannot be empty;
 *
 * - kConflicts, the negative table: as kExtension, but the values of the
 *   scope must be none of the tuples;
 *
 * - kLinear: with k variables in the scope, the parameters are k
 *   coefficients a_i, then an operator op [see linear_op], and then a
 *   right-hand side b, and \f$\sum_i a_i v_i \; op \; b\f$ must hold (with
 *   an empty scope, 0 op b);
 *
 * - kAllDifferent: no parameters, the variables of the scope must take
 *   pairwise different values;
 *
 * - kElement: the scope is ( i , z , x_0 , ... , x_{k - 1} ), with k at
 *   least 1, and the parameter is an offset o: \f$o \leq i < o + k\f$ and
 *   \f$x_{i - o} = z\f$ must hold;
 *
 * - kCumulative: the scope holds the start times s_0, ..., s_{k - 1} of k
 *   tasks, and the parameters are their durations d_0, ..., d_{k - 1}, their
 *   resource usages r_0, ..., r_{k - 1}, all non-negative, and the capacity
 *   b: at every time t, the sum of the r_j of the tasks with
 *   \f$s_j \leq t < s_j + d_j\f$ must be at most b.
 *
 * The scopes and the parameters of all the constraints are stored one
 * constraint after the other [see RaggedArray].
 *
 * Each constraint has a weight: one of infinite weight is *hard*, i.e., it
 * must be satisfied, while one of finite (non-negative) weight is *soft*,
 * i.e., it may be violated at the cost of its weight. A constraint may also
 * be *inactive*, in which case it is as if it were not there; it is not
 * removed, so that the indices of the others do not change. Each value d
 * of each variable v may have a *cost* c_{v,d}, of any sign, paid if v takes
 * d: this covers a linear objective \f$\sum_v c_v v\f$ as well as the unary
 * cost functions of the weighted CSP, and is where the Lagrangian terms of
 * a decomposition land. The problem is then to find values of the variables
 * in their current domains satisfying all the active hard constraints and
 * minimizing the costs of the values taken plus the weights of the active
 * soft constraints violated. A non-linear objective (e.g., a makespan) is
 * an auxiliary variable, with its constraints, and costs on its values.
 *
 * All this is the "physical representation" of the CSPBlock, together with
 * the groups the variables may be dealt out to [see set_variable_groups()],
 * which are what a structure of the CSPBlock is to be made of. The
 * "abstract representation" is not defined yet: the generate_abstract_*()
 * methods of Block do nothing, and the CSPBlock has no Variable, Constraint
 * or Objective. Hence, its solutions are vectors of values, one per variable,
 * which is_feasible( std::span< const int > ) and get_objective_value()
 * check against the physical representation.
 *
 * A CSPBlock can be load()-ed from a subset of the FlatZinc format, and
 * deserialize()-d and serialize()-d out of and into a netCDF group [see
 * serialize()]. FlatZinc is chosen as the text format because it is what
 * MiniZinc produces for a solver once it has flattened a model, keeping
 * native the global constraints the solver declares, and because it is read
 * one statement at a time with no further library, whereas XCSP3 is XML. */

class CSPBlock : public Block
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/

 using v_int = std::vector< int >;        ///< a vector of int

 using v_v_int = std::vector< v_int >;    ///< a vector of vectors of int

 using v_dbl = std::vector< double >;     ///< a vector of double

 using c_v_dbl = const v_dbl;             ///< a const vector of double

 using v_flag = std::vector< unsigned char >;  ///< a vector of 0/1 flags

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// rows of ints of different length, one after the other
 /** Rows of int of different length, stored as a single vector of the
  * entries of all of them, one row after the other, and the index where
  * each one starts in it, which takes a fraction of the memory of a
  * std::vector per row; the i-th row is a std::span of its entries. This
  * holds the universes of the variables, the scopes and the parameters of
  * the constraints. */

 class RaggedArray {
  public:
  /// a row: the span of its entries
  using value_type = std::span< const int >;

  /// returns the number of rows
  [[nodiscard]] std::size_t size( void ) const {
   return( v_start.size() - 1 );
   }

  /// returns true if there is no row
  [[nodiscard]] bool empty( void ) const { return( size() == 0 ); }

  /// returns the i-th row, the span of its entries
  [[nodiscard]] value_type operator[]( std::size_t i ) const {
   return( value_type( v_val.data() + v_start[ i ] ,
		       v_start[ i + 1 ] - v_start[ i ] ) );
   }

  /// returns true if the rows are the same, in the same order
  [[nodiscard]] bool operator==( const RaggedArray & o ) const {
   return( ( v_val == o.v_val ) && ( v_start == o.v_start ) );
   }

  /// returns the entries of all the rows, one row after the other
  [[nodiscard]] const v_int & values( void ) const { return( v_val ); }

  /// returns where each row starts in values(), and then its size
  [[nodiscard]] const std::vector< std::size_t > & starts( void ) const {
   return( v_start );
   }

  /// adds a row after the existing ones
  void push_back( std::span< const int > r ) {
   v_val.insert( v_val.end() , r.begin() , r.end() );
   v_start.push_back( v_val.size() );
   }

  /// replaces the rows with those given by their entries and starts
  /** Replaces the rows with those whose entries are \p val and where the
   * i-th starts at \p start[ i ], \p start having one more element than
   * the rows, the last being the size of \p val; an exception is thrown if
   * \p start is not so. */
  void assign( v_int && val , std::vector< std::size_t > && start );

  /// removes all the rows, and frees their memory
  void clear( void ) {
   v_int().swap( v_val );
   std::vector< std::size_t >( 1 , 0 ).swap( v_start );
   }

  private:
  v_int v_val;                                ///< the entries
  std::vector< std::size_t > v_start = { 0 }; ///< where each row starts
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the kinds of the constraints [see the class]

 enum constraint_kind {
  kExtension = 0 ,  ///< positive table
  kConflicts ,      ///< negative table
  kLinear ,         ///< linear (in)equality
  kAllDifferent ,   ///< all different
  kElement ,        ///< x_{i - o} = z
  kCumulative ,     ///< tasks sharing a resource of given capacity
  kNumberKinds      ///< the number of kinds
  };

 /// the operators of a kLinear constraint

 enum linear_op {
  kLE = 0 ,  ///< \f$\sum_i a_i v_i \leq b\f$
  kEQ ,      ///< \f$\sum_i a_i v_i = b\f$
  kNE        ///< \f$\sum_i a_i v_i \neq b\f$
  };

/*------------------------------ CONSTRUCTOR -------------------------------*/
 /// constructor of CSPBlock, taking a pointer to the father Block
 /** Constructor of CSPBlock. It accepts a pointer to the father Block
  * (defaulting to nullptr, both because the root Block has no father and so
  * that this can also be used as the void constructor required by the Block
  * factory). The CSPBlock is empty: no variable and no constraint. */

 explicit CSPBlock( Block * father = nullptr ) : Block( father ) {}

/*------------------------------- DESTRUCTOR -------------------------------*/
 /// destructor of CSPBlock

 ~CSPBlock() override = default;

/*------------------------- OTHER INITIALIZATIONS --------------------------*/
 /// loads the CSPBlock out of the given universes and constraints
 /** Loads the CSPBlock out of the universes of the variables in \p
  * universes, one per variable, and of the constraints whose kinds,
  * scopes and parameters are in \p kinds, \p scopes and \p params, with
  * the weights in \p weights, replacing whatever the CSPBlock held. The
  * values of a universe are sorted, repeated ones being kept once, and a
  * universe cannot be empty; the current domains are the universes, the
  * costs are 0 and all the constraints are active. If \p weights is empty
  * all the constraints are hard, otherwise it must have one weight per
  * constraint, +INF for a hard one and a finite non-negative value for a
  * soft one. An exception is thrown, and nothing changes, if anything is
  * wrong [see the class]. */

 void load( v_v_int && universes , v_int && kinds , v_v_int && scopes ,
	    v_v_int && params , v_dbl && weights = {} );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// loads the CSPBlock out of an istream in the FlatZinc format
 /** Loads the CSPBlock out of an istream in a subset of the FlatZinc
  * format, \p frmt being 'F' or 0. The statements read are the parameter
  * declarations of int and arrays of int, the declarations of the
  * variables of int (with a range or a set of values) and bool (the
  * values 0 and 1) type and of the arrays of them, the constraints below
  * and the solve item; the predicate declarations are skipped, and anything
  * else (e.g., a variable of type int with no range) throws exception. A
  * variable declared equal to another is the same variable, one declared
  * equal to a value has that value only in its current domain; a value in
  * an array of variables, where a variable is expected, is a variable of
  * that value only, made once per value.
  *
  * The constraints are int_lin_le, int_lin_eq and int_lin_ne (kLinear);
  * int_le, int_lt, int_eq, int_ne, int_ge, int_gt and bool2int (kLinear
  * over two variables); fzn_all_different_int and all_different_int
  * (kAllDifferent); array_var_int_element (kElement of offset 1);
  * array_int_element and set_in (kExtension over the index and the value,
  * and over the variable); fzn_table_int and table_int (kExtension); and
  * fzn_cumulative and cumulative, with durations, usages and capacity that
  * are parameters (kCumulative). Moreover, three constraints are read that
  * are not in FlatZinc but are what print() writes: csp_table_int_neg( x ,
  * t ) (kConflicts), csp_element_offset( i , o , x , z ) (kElement of
  * offset o) and csp_lin_le/csp_lin_eq/csp_lin_ne, as the int_lin_* ones.
  *
  * The solve item "minimize v" ("maximize v") gives each value d of the
  * variable v the cost d (- d). The following annotations, which print()
  * writes, are also read: universe( S ) on a variable, its universe being
  * the set S and its current domain the declared one; cost( [ c ... ] ) on
  * a variable, the costs of the values of its universe; group( g ) on a
  * variable, its group [see set_variable_groups()]; weight( w ) on a
  * constraint, which is soft with weight w; and inactive on a constraint.
  * The other annotations are skipped. */

 void load( std::istream & input , char frmt = 0 ) override;

 using Block::load;  ///< the load() from a file name

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the format of
  * CSPBlock [see serialize()]. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// sets the groups of the variables, part of the physical representation
 /** Sets the group of each variable, \p groups having one non-negative
  * entry per variable, or none to say that there are no groups; the groups
  * are numbered from 0 to the largest entry, and a group may be empty. The
  * groups are what a structure of the CSPBlock is to be built out of, as
  * in SATBlock; the structure itself waits for the abstract
  * representation. */

 void set_variable_groups( v_int && groups );

/*--------------------- Methods for handling Variable ----------------------*/
 /// builds the abstract Variable: one integer ColVariable per variable
 /** Builds the abstract representation of the variables [see
  * generate_abstract_constraints() for the constraints], one integer
  * ColVariable per variable of the CSPBlock, and registers them with
  * add_static_variable(). Each ColVariable is made integer; its numeric
  * bounds, i.e., the min and the max of its current domain, are a box
  * OneVarConstraint built by generate_abstract_constraints(), since a
  * ColVariable carries only sign and unit bounds in its type, not an
  * arbitrary min..max [see ColVariable.h]. The method guards against being
  * called twice, as the Block protocol asks. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*-------------------- Methods for handling Constraint ---------------------*/
 /// builds the abstract Constraint: the bounds and one per active hard one
 /** Builds the abstract representation of the constraints, calling
  * generate_abstract_variables() first if it has not been called. It builds
  * (i) a box OneVarConstraint per variable, holding the min..max of its
  * current domain, and (ii) one abstract Constraint per active and hard
  * constraint of the physical catalogue: a kLinear with the kLE or the kEQ
  * operator is an FRowConstraint over a LinearFunction, the globals and a
  * kLinear with the kNE operator are the module-local :CSPConstraint [see
  * CSPGlobalConstraints.h], each reusing is_satisfied() as its verifier. A
  * soft constraint is not built here, since its penalty belongs to the
  * objective, which phase 1 does not generate [see generate_objective()]. */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------- Methods for handling Objective ---------------------*/
 /// the abstract objective is not defined yet: does nothing
 /** The abstract objective (the costs of the values and the weights of the
  * soft constraints) is left to a later phase [see doc/cp-design.md]; phase 1
  * only generates the Variable and the Constraint. */

 void generate_objective( Configuration * = nullptr ) override {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the abstract integer ColVariable, one per variable
 /** Returns the integer ColVariable of the abstract representation, one per
  * variable of the CSPBlock, which generate_abstract_variables() builds; it
  * is empty until then. It is not const, so that a Solver, or a test, can set
  * their values before checking the abstract representation [see
  * is_feasible( bool )]. */

 [[nodiscard]] std::vector< ColVariable > & get_abstract_variables( void ) {
  return( v_var );
  }

/*------------------ Methods for reading the data of the CSPBlock ----------*/
 /// the sense of the objective, which is minimization

 [[nodiscard]] int get_objective_sense( void ) const override {
  return( Objective::eMin );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of variables

 [[nodiscard]] unsigned int get_number_variables( void ) const {
  return( v_universe.size() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of values of all the universes, i.e., of slots

 [[nodiscard]] unsigned int get_number_values( void ) const {
  return( v_universe.values().size() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the number of constraints, the inactive ones included

 [[nodiscard]] unsigned int get_number_constraints( void ) const {
  return( v_kind.size() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the universes of the variables, sorted

 [[nodiscard]] const RaggedArray & get_universes( void ) const {
  return( v_universe );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the universe of the i-th variable, sorted

 [[nodiscard]] std::span< const int > get_universe( Index i ) const {
  return( v_universe[ i ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the slot of the value d of the i-th variable
 /** Returns the slot of the value d of the i-th variable, i.e., its position
  * in get_universes().values(), or Inf< Index >() if d is not in its
  * universe. */

 [[nodiscard]] Index value_index( Index i , int d ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the current domains: 1 for each slot whose value is there

 [[nodiscard]] const v_flag & get_domain_mask( void ) const {
  return( v_domain );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the current domain of the i-th variable, sorted

 [[nodiscard]] v_int get_domain( Index i ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns true if d is in the current domain of the i-th variable

 [[nodiscard]] bool in_domain( Index i , int d ) const {
  const auto s = value_index( i , d );
  return( ( s != Inf< Index >() ) && v_domain[ s ] );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the costs of the values, one per slot

 [[nodiscard]] c_v_dbl & get_costs( void ) const { return( v_costs ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns true if some value has a nonzero cost

 [[nodiscard]] bool has_costs( void ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the kinds of the constraints, constraint_kind values

 [[nodiscard]] const v_int & get_kinds( void ) const { return( v_kind ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the scopes of the constraints

 [[nodiscard]] const RaggedArray & get_scopes( void ) const {
  return( v_scope );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the parameters of the constraints

 [[nodiscard]] const RaggedArray & get_parameters( void ) const {
  return( v_param );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the weights of the constraints, +INF for the hard ones

 [[nodiscard]] c_v_dbl & get_weights( void ) const { return( v_weights ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns true if the k-th constraint is hard

 [[nodiscard]] bool is_hard( Index k ) const {
  return( v_weights[ k ] == Inf< double >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns true if all the constraints are hard

 [[nodiscard]] bool all_hard( void ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns 1 for each active constraint and 0 for each inactive one

 [[nodiscard]] const v_flag & get_active( void ) const { return( v_active ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns true if the k-th constraint is active

 [[nodiscard]] bool is_active( Index k ) const { return( v_active[ k ] ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the groups of the variables, empty if there are none

 [[nodiscard]] const v_int & get_variable_groups( void ) const {
  return( v_group );
  }

/*--------------------------------------------------------------------------*/
/*------------- Methods for checking the state of the CSPBlock -------------*/
 /// returns true if the values in \p a satisfy the k-th constraint
 /** Returns true if the values in \p a, one per variable of the CSPBlock,
  * satisfy the k-th constraint, whatever its weight and whether it is
  * active or not, the current domains being left out [see the class]; an
  * exception is thrown if \p a does not have one value per variable. */

 [[nodiscard]] bool is_satisfied( Index k , std::span< const int > a ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns true if the values in \p a are a feasible solution
 /** Returns true if each value in \p a, one per variable of the CSPBlock,
  * is in the current domain of its variable and the values satisfy every
  * active hard constraint; an exception is thrown if \p a does not have
  * one value per variable. */

 [[nodiscard]] bool is_feasible( std::span< const int > a ) const;

 using Block::is_feasible;  ///< the one of the abstract representation

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns true if the abstract representation is feasible at its values
 /** With \p useabstract false this is the inherited Block check, i.e., the
  * sub-Blocks, of which the CSPBlock has none. With \p useabstract true it
  * checks the abstract representation [see generate_abstract_constraints()]
  * at the current values of its ColVariable: it compute()s and feasible()s
  * the box of each variable and each generated abstract Constraint, and
  * returns their conjunction; an exception is thrown if the abstract
  * constraints are not generated yet. The two differ from
  * is_feasible( std::span ) only where the physical representation is richer
  * than the abstract one, i.e., a current domain with holes, which the box
  * widens to its min..max, and a soft constraint, which is not generated. */

 bool is_feasible( bool useabstract = false ,
		   Configuration * fsbc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the sum of the weights of the active soft constraints violated
 /** Returns the sum of the weights of the active soft constraints that the
  * values in \p a, one per variable, violate, the hard ones being left
  * out whether they are satisfied or not [see is_feasible()]. */

 [[nodiscard]] double get_violated_weight( std::span< const int > a ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the value of the objective at the values in \p a
 /** Returns the value of the objective at the values in \p a, one per
  * variable, i.e., the sum of the costs of the values taken plus
  * get_violated_weight( a ); an exception is thrown if a value is not in
  * the universe of its variable, which has no cost for it. */

 [[nodiscard]] double get_objective_value( std::span< const int > a ) const;

/*-------------------------- Methods for printing --------------------------*/
 /// prints the CSPBlock on an ostream with the given verbosity
 /** Prints the CSPBlock: the number of variables and of constraints, and,
  * with \p vlvl == 'C', the whole instance in the FlatZinc format, which is
  * what load() reads back. The i-th variable is named x<i>, its current
  * domain is declared and its universe, if larger, is the annotation
  * universe(); if the only costs are those of a linear objective in a
  * single variable, with coefficient 1 or - 1, the solve item minimizes or
  * maximizes that variable, otherwise the costs are the annotations cost()
  * of the variables and the solve item is satisfy. A kExtension and a
  * kAllDifferent constraint are written as fzn_table_int and
  * fzn_all_different_int, a kElement of offset 1 as array_var_int_element,
  * a kCumulative as fzn_cumulative, a kLinear as int_lin_le, int_lin_eq or
  * int_lin_ne; a kConflicts and a kElement of another offset as the
  * constraints of CSPBlock that load() reads [see load()]. */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/*----------- METHODS FOR MODIFYING THE PHYSICAL REPRESENTATION ------------*/
/** @name Methods for modifying the physical representation
 *
 * These methods change the physical representation of the CSPBlock,
 * issuing the CSPBlockMod [see below] as \p issueMod says [see
 * Observer::make_par()]. Since the CSPBlock has no abstract representation
 * yet, there is nothing else to change. Nothing is ever removed, so that
 * the indices of the variables, of the slots and of the constraints never
 * change for those that are there, and a Modification can be read after
 * other ones have been issued.
 *  @{ */

 /// changes the current domains of the variables in the Subset
 /** Changes the current domains of the variables in the Subset \p nms, the
  * new one of the variable nms[ k ] being the values in \p domains[ k ],
  * which must be in its universe (repeated ones being kept once); \p ordered
  * tells if \p nms is ordered increasingly, which it is made to be
  * otherwise, together with \p domains. The CSPBlockSbstMod issued has
  * type eChgDomain and the variables whose domains do change, and it says
  * that the feasible region shrinks if a value has been removed and that
  * it grows if one has been added. An exception is thrown, and nothing
  * changes, if a variable is not there, is repeated, or a value is not in
  * its universe. */

 void chg_domains( Subset && nms , v_v_int && domains , bool ordered = false ,
		   ModParam issueMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the current domains of the variables in the Range
 /** Same as chg_domains( Subset ), for the variables in the Range \p rng,
  * the new domain of the variable i being \p domains[ i - rng.first ]. */

 void chg_domains( v_v_int && domains , Range rng = INFRange ,
		   ModParam issueMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// adds values to the universes of the variables in the Subset
 /** Adds the values in \p values[ k ] to the universe of the variable
  * nms[ k ] of the Subset \p nms, and to its current domain, with cost 0;
  * the values already in the universe are left as they are. The slots of
  * the values after them move, the slot of a value being its position in
  * the universes [see value_index()]. The CSPBlockSbstMod issued has type
  * eAddValues and the variables whose universes do grow. An exception is
  * thrown, and nothing changes, if a variable is not there or is repeated,
  * or if \p values is shorter than \p nms. */

 void add_values( Subset && nms , v_v_int && values ,
		  ModParam issueMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the costs of the values of the slots in the Range
 /** Changes the costs of the values of the slots in the Range \p rng [see
  * value_index()], the new ones being in \p NCost. The CSPBlockRngdMod
  * issued has type eChgCost and the Range of the slots. An exception is
  * thrown if \p NCost is shorter than \p rng or a cost is not finite. */

 void chg_costs( MF_dbl_sp NCost , Range rng = INFRange ,
		 ModParam issueMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the costs of the values of the slots in the Subset
 /** Same as chg_costs( Range ), for the slots in the Subset \p nms; \p
  * ordered tells if \p nms is ordered increasingly, which it is made to be
  * otherwise, together with \p NCost. */

 void chg_costs( MF_dbl_sp NCost , Subset && nms , bool ordered = false ,
		 ModParam issueMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the weights of the constraints in the Range
 /** Changes the weights of the constraints in the Range \p rng, the new
  * ones being in \p NWeight, +INF making a constraint hard. The
  * CSPBlockRngdMod issued has type eChgWeight and the Range of the
  * constraints, and it says that the feasible region shrinks if a
  * constraint turns hard and that it grows if one turns soft. An exception
  * is thrown if \p NWeight is shorter than \p rng or a weight is wrong
  * [see load()]. */

 void chg_weights( MF_dbl_sp NWeight , Range rng = INFRange ,
		   ModParam issueMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the weights of the constraints in the Subset
 /** Same as chg_weights( Range ), for the constraints in the Subset \p
  * nms; \p ordered tells if \p nms is ordered increasingly, which it is
  * made to be otherwise, together with \p NWeight. */

 void chg_weights( MF_dbl_sp NWeight , Subset && nms , bool ordered = false ,
		   ModParam issueMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// makes the constraints in the Range active or inactive
 /** Makes the constraints in the Range \p rng active if \p active is true
  * and inactive otherwise. The CSPBlockRngdMod issued has type eChgActive
  * and the Range of the constraints, and it says what does change: the
  * feasible region shrinks (grows) if a hard constraint turns active
  * (inactive), the objective may go up (down) if a soft one does. */

 void chg_active( bool active , Range rng = INFRange ,
		  ModParam issueMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// makes the constraints in the Subset active or inactive
 /** Same as chg_active( Range ), for the constraints in the Subset \p nms;
  * \p ordered tells if \p nms is ordered increasingly, which it is made to
  * be otherwise. */

 void chg_active( bool active , Subset && nms , bool ordered = false ,
		  ModParam issueMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// adds the given constraints after the existing ones
 /** Adds the constraints whose kinds, scopes and parameters are in \p
  * kinds, \p scopes and \p params after the existing ones, active and with
  * the weights in \p weights, all hard if it is empty [see load()]. The
  * CSPBlockRngdMod issued has type eAddConstraints and the Range of the
  * indices of the new constraints. An exception is thrown, and nothing is
  * added, if any of them is wrong. */

 void add_constraints( v_int && kinds , v_v_int && scopes ,
		       v_v_int && params , v_dbl && weights = {} ,
		       ModParam issueMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// adds variables after the existing ones
 /** Adds one variable after the existing ones per element of \p universes,
  * its universe, which is also its current domain, with costs 0 and in no
  * constraint, which the constraints added after can then use [see
  * add_constraints()]. If the variables have groups, \p groups must give
  * the group of each new one, otherwise it must be empty. The
  * CSPBlockRngdMod issued has type eAddVariables and the Range of the
  * indices of the new variables. An exception is thrown, and nothing is
  * added, if a universe is empty or the groups are wrong. */

 void add_variables( v_v_int && universes , v_int && groups = {} ,
		     ModParam issueMod = eNoBlck );

/** @} ---------------------------------------------------------------------*/

/*------------------------- Methods for serializing ------------------------*/
 /// serializes the CSPBlock into a netCDF group
 /** Serializes the CSPBlock into a netCDF group, besides the "type"
  * attribute written by Block::serialize(), as:
  *
  * - the dimension "NumberVariables", the number n of the variables, and
  *   the dimension "NumberValues", the number of the slots;
  *
  * - the int variable "Universe", over "NumberValues", the values of the
  *   universes one variable after the other, and the uint variable
  *   "UniverseStart", over "NumberVariables", the slot of the first value
  *   of each variable;
  *
  * - the ubyte variable "Domain", over "NumberValues", 1 for the values in
  *   the current domains and 0 for the others; it is absent if each current
  *   domain is the universe;
  *
  * - the double variable "Costs", over "NumberValues", the costs of the
  *   values; it is absent if they are all 0;
  *
  * - the int variable "VariableGroups", over "NumberVariables", the groups
  *   of the variables [see set_variable_groups()]; it is absent if there are
  *   none;
  *
  * - the dimension "NumberConstraints", the number m of the constraints,
  *   and the int variable "ConstraintKind" over it, the constraint_kind of
  *   each one;
  *
  * - the int variable "Scope", over the dimension "NumberScopeEntries", the
  *   scopes one constraint after the other, and the uint variable
  *   "ScopeStart", over "NumberConstraints", where each one starts; both
  *   are absent if all the scopes are empty;
  *
  * - the int variable "Parameters", over the dimension
  *   "NumberParameterEntries", and the uint variable "ParametersStart", the
  *   same for the parameters;
  *
  * - the double variable "Weights", over "NumberConstraints", the weights
  *   of the constraints, +INF for the hard ones; it is absent if all the
  *   constraints are hard;
  *
  * - the ubyte variable "Active", over "NumberConstraints", 1 for the
  *   active constraints and 0 for the inactive ones; it is absent if all
  *   the constraints are active.
  *
  * Everything about the constraints is absent if there are none. */

 void serialize( netCDF::NcGroup & group ) const override;

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// checks a constraint against the variables, throwing if it is wrong

 void check_constraint( int kind , std::span< const int > scope ,
			std::span< const int > param ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// checks the weights from the first-th on, throwing if one is wrong

 void check_weights( const v_dbl & weights , Index first = 0 ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// checks that \p a has one value per variable, throwing otherwise

 void check_assignment( std::span< const int > a ,
			const char * method ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// changes the activity of the k-th constraint, returning what changes

 Modification::ModConcern set_active( Index k , bool active );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// empties the CSPBlock

 void guts_of_destructor( void );

/*---------------------------- PROTECTED FIELDS ----------------------------*/

 RaggedArray v_universe;  ///< the universes of the variables, sorted

 v_flag v_domain;         ///< 1 for each slot in the current domain

 v_dbl v_costs;           ///< the cost of each slot

 v_int v_group;           ///< the group of each variable, if any

 v_int v_kind;            ///< the kind of each constraint

 RaggedArray v_scope;     ///< the scopes of the constraints

 RaggedArray v_param;     ///< the parameters of the constraints

 v_dbl v_weights;         ///< the weights of the constraints, +INF if hard

 v_flag v_active;         ///< 1 for each active constraint

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 // the abstract representation [see generate_abstract_variables() and
 // generate_abstract_constraints()]; v_var is declared before the
 // constraints so that, on destruction, the constraints (destroyed first,
 // in reverse order of declaration) unregister from ColVariable still alive

 unsigned char f_AR = 0;  ///< what of the abstract representation is there

 static constexpr unsigned char HasVar = 1;  ///< the Variable are there
 static constexpr unsigned char HasCon = 2;  ///< the Constraint are there

 std::vector< ColVariable > v_var;  ///< one integer ColVariable per variable

 std::vector< BoxConstraint > v_box;  ///< the min..max bound of each variable

 std::vector< FRowConstraint > v_lin;  ///< the kLinear kLE/kEQ constraints

 std::vector< AllDifferentConstraint > v_alldiff;  ///< the kAllDifferent ones

 std::vector< TableConstraint > v_table;  ///< the kExtension/kConflicts ones

 std::vector< ElementConstraint > v_element;  ///< the kElement ones

 std::vector< CumulativeConstraint > v_cumul;  ///< the kCumulative ones

 std::vector< NotEqualConstraint > v_ne;  ///< the kLinear kNE constraints

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*---------------------------- PRIVATE METHODS -----------------------------*/

 SMSpp_insert_in_factory_h;  // insert CSPBlock in the Block factory

/*--------------------------------------------------------------------------*/

 };  // end( class( CSPBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS CSPBlockMod ----------------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Modification for the changes of a CSPBlock
/** Derived class from Modification for the changes of the physical
 * representation of a CSPBlock: its type says what has changed, and the
 * derived classes CSPBlockRngdMod and CSPBlockSbstMod which variables (for
 * eChgDomain, eAddValues and eAddVariables), slots (for eChgCost) or
 * constraints (for the others). */

class CSPBlockMod : public Modification
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/
 /// the types of CSPBlockMod

 enum CSPBlock_mod_type {
  eChgDomain = 0 ,   ///< the current domains of some variables have changed
  eAddValues ,       ///< values have been added to some universes
  eChgCost ,         ///< the costs of some values have changed
  eChgWeight ,       ///< the weights of some constraints have changed
  eAddConstraints ,  ///< some constraints have been added
  eChgActive ,       ///< some constraints have turned (in)active
  eAddVariables      ///< some variables have been added
  };

 /// what a Modification of the given type changes
 /** What a Modification of the given type changes of the model [see
  * Modification::ModConcern], which is what changes() says and what the
  * methods that issue it ask before issuing it [see Observer::issue_pmod()];
  * the types whose effect depends on the change say all that may happen,
  * and the methods that issue them say which does: a change of the domains
  * makes the feasible region shrink and/or grow, adding values makes it
  * grow, the costs are the objective, the weights are the objective but a
  * constraint turning hard (soft) makes the region shrink (grow), adding
  * constraints makes it shrink, a constraint turning active or inactive may
  * do anything, and so may adding variables. */

 static ModConcern changes_of( int type ) {
  switch( type ) {
   case( eChgDomain ):
    return( eModPhys | eRegnShrink | eRegnGrow );
   case( eAddValues ):
    return( eModPhys | eRegnGrow );
   case( eChgCost ):
    return( ePhysObj );
   case( eChgWeight ):
    return( ePhysObj | eModCnsSide | eRegnShrink | eRegnGrow );
   case( eAddConstraints ):
    return( eModPhys | eModCnsSet | eRegnShrink | eObjUp );
   case( eChgActive ):
    return( eModPhys | eModCnsSide | eEffAny );
   default:
    return( eModPhys | eEffAny );
   }
  }

/*------------------------ CONSTRUCTOR & DESTRUCTOR ------------------------*/
 /// constructor: takes the CSPBlock and the type

 CSPBlockMod( CSPBlock * const fblock , int type )
  : Modification( changes_of( type ) ) , f_Block( fblock ) , f_type( type )
 {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// destructor, does nothing

 ~CSPBlockMod() override = default;

/*----------------------- PUBLIC METHODS OF THE CLASS ----------------------*/
 /// returns the CSPBlock the CSPBlockMod refers to

 Block * get_Block( void ) const override { return( f_Block ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the type of the CSPBlockMod

 [[nodiscard]] int type( void ) const { return( f_type ); }

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// prints the CSPBlockMod

 void print( std::ostream & output ) const override {
  output << "CSPBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eChgDomain ):      output << "change domains "; break;
   case( eAddValues ):      output << "add values "; break;
   case( eChgCost ):        output << "change costs "; break;
   case( eChgWeight ):      output << "change weights "; break;
   case( eAddConstraints ): output << "add constraints "; break;
   case( eChgActive ):      output << "change activity "; break;
   case( eAddVariables ):   output << "add variables "; break;
   default:                 output << "type " << f_type << " ";
   }
  }

/*---------------------------- PROTECTED FIELDS ----------------------------*/

 CSPBlock * f_Block;  ///< the CSPBlock the CSPBlockMod refers to

 int f_type;          ///< the type of the CSPBlockMod

/*--------------------------------------------------------------------------*/

 };  // end( class( CSPBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS CSPBlockRngdMod --------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from CSPBlockMod for the changes of a Range

class CSPBlockRngdMod : public CSPBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*------------------------ CONSTRUCTOR & DESTRUCTOR ------------------------*/
 /// constructor: takes the CSPBlock, the type and the Range

 CSPBlockRngdMod( CSPBlock * const fblock , int type , Block::Range rng )
  : CSPBlockMod( fblock , type ) , f_rng( rng ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// destructor, does nothing

 ~CSPBlockRngdMod() override = default;

/*----------------------- PUBLIC METHODS OF THE CLASS ----------------------*/
 /// returns the Range

 [[nodiscard]] Block::c_Range & rng( void ) const { return( f_rng ); }

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// prints the CSPBlockRngdMod

 void print( std::ostream & output ) const override {
  CSPBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

/*---------------------------- PROTECTED FIELDS ----------------------------*/

 Block::Range f_rng;  ///< the Range

/*--------------------------------------------------------------------------*/

 };  // end( class( CSPBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS CSPBlockSbstMod --------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from CSPBlockMod for the changes of a Subset

class CSPBlockSbstMod : public CSPBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*------------------------ CONSTRUCTOR & DESTRUCTOR ------------------------*/
 /// constructor: takes the CSPBlock, the type and the Subset, "consumed"

 CSPBlockSbstMod( CSPBlock * const fblock , int type , Block::Subset && nms )
  : CSPBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// destructor, does nothing

 ~CSPBlockSbstMod() override = default;

/*----------------------- PUBLIC METHODS OF THE CLASS ----------------------*/
 /// returns the Subset, ordered increasingly

 [[nodiscard]] Block::c_Subset & nms( void ) const { return( f_nms ); }

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// prints the CSPBlockSbstMod

 void print( std::ostream & output ) const override {
  CSPBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

/*---------------------------- PROTECTED FIELDS ----------------------------*/

 Block::Subset f_nms;  ///< the Subset

/*--------------------------------------------------------------------------*/

 };  // end( class( CSPBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS CSPBlockChange --------------------------*/
/*--------------------------------------------------------------------------*/
/// a Change restricting or restoring the current domains of a CSPBlock
/** A CSPBlockChange changes the current domains of some variables of a
 * CSPBlock [see CSPBlock::chg_domains()], as a branching does: with the
 * value d_k given for the k-th variable v_k, the type eAssign makes the
 * domain of v_k the value d_k only (v_k = d_k), eRemove takes d_k out of it
 * (v_k != d_k), eSplitLE keeps only its values not larger than d_k
 * (v_k <= d_k) and eSplitGT only those larger than d_k (v_k > d_k); a
 * value d_k that is not in the domain leaves it empty under eAssign and
 * unchanged under eRemove. The type eSetDomain makes the domain of each v_k
 * the values given for it, which is what the undo of the others is.
 *
 * With doUndo, apply() returns the CSPBlockChange of type eSetDomain that
 * brings the domains back. A variable cannot be repeated.
 *
 * In netCDF a CSPBlockChange is a group with, besides the "type" of every
 * Change, the attribute "CSPBlockChange_type", the uint variable "Index"
 * along the dimension "NumVar", the variables, and the int variable "Value"
 * along the dimension "NumValues", the values: one per variable, or for
 * eSetDomain the domains one after the other, with the uint variable
 * "ValueStart" along "NumVar" saying where each starts; "Value" is absent
 * if there are no values. */

class CSPBlockChange : public Change
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/
 /// the types of CSPBlockChange

 enum CSPBlock_chg_type {
  eAssign = 0 ,  ///< v_k = d_k
  eRemove ,      ///< v_k != d_k
  eSplitLE ,     ///< v_k <= d_k
  eSplitGT ,     ///< v_k > d_k
  eSetDomain ,   ///< the domain of v_k is the values given
  eEmpty         ///< nothing, for the factory
  };

/*------------------------ CONSTRUCTOR & DESTRUCTOR ------------------------*/
 /// constructor: takes the type, the variables and the values, "consumed"
 /** Constructor: takes the type, the variables in \p nms and the values in
  * \p values, one per variable except for eSetDomain, whose values are the
  * domains one after the other, the k-th having \p sizes[ k ] values. */

 explicit CSPBlockChange( int type = eEmpty , Block::Subset && nms = {} ,
			  std::vector< int > && values = {} ,
			  std::vector< unsigned int > && sizes = {} )
  : f_type( type ) , f_nms( std::move( nms ) ) ,
    f_values( std::move( values ) ) , f_sizes( std::move( sizes ) ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// destructor, does nothing

 ~CSPBlockChange() override = default;

/*----------------------- PUBLIC METHODS OF THE CLASS ----------------------*/
 /// applies the CSPBlockChange to a CSPBlock, returning the undo if asked

 Change * apply( Block * block , bool doUndo = false ,
		 ModParam issueMod = eNoBlck ,
		 ModParam issueAMod = eNoBlck ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// reads the CSPBlockChange out of a netCDF group

 void deserialize( const netCDF::NcGroup & group ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// writes the CSPBlockChange into a netCDF group

 void serialize( netCDF::NcGroup & group ) const override;

 using Change::serialize;  // keep the other serialize() visible

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the type of the CSPBlockChange

 [[nodiscard]] int type( void ) const { return( f_type ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the variables

 [[nodiscard]] Block::c_Subset & nms( void ) const { return( f_nms ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the values

 [[nodiscard]] const std::vector< int > & values( void ) const {
  return( f_values );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the sizes of the domains (for eSetDomain)

 [[nodiscard]] const std::vector< unsigned int > & sizes( void ) const {
  return( f_sizes );
  }

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// prints the CSPBlockChange

 void print( std::ostream & output ) const override {
  static const char * name[] = { "assign " , "remove " , "split <= " ,
				 "split > " , "set domain " , "empty " };
  output << "CSPBlockChange[" << this << "]: "
	 << name[ ( f_type >= eAssign ) && ( f_type <= eEmpty ) ? f_type
				                                : eEmpty ]
	 << f_nms.size() << " variables" << std::endl;
  }

/*---------------------------- PROTECTED FIELDS ----------------------------*/

 int f_type;                         ///< the type

 Block::Subset f_nms;                ///< the variables

 std::vector< int > f_values;        ///< the values

 std::vector< unsigned int > f_sizes;  ///< the sizes of the domains

/*------------------------ PRIVATE PART OF THE CLASS -----------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( CSPBlockChange ) )

/** @} end( group( CSPBlock_CLASSES ) ) */

/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* CSPBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File CSPBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/

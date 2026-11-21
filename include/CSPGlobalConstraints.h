/*--------------------------------------------------------------------------*/
/*----------------------- File CSPGlobalConstraints.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *abstract* representation of the global constraints of
 * the CSPBlock [see CSPBlock.h]. The physical representation of the CSPBlock
 * is a compact catalogue of constraints, each a kind [see
 * CSPBlock::constraint_kind], a scope and some parameters; the abstract
 * representation of a global, following section 3 of doc/cp-design.md, is a
 * concrete :Constraint over the integer ColVariable of the CSPBlock, giving
 * two things: a verifier feasible(), by which any Solver checks an assignment
 * cheaply, and a linear relaxation, by which a Solver that cannot digest the
 * global natively still has a weak but valid linear view of it.
 *
 * The verifier is the same one the physical representation already has: each
 * class here reuses CSPBlock::is_satisfied() for its constraint index, so the
 * abstract feasible() and the physical is_feasible() agree by construction.
 * The linear case of the catalogue (kLinear) is already an FRowConstraint and
 * is not here; the classes here are the global ones, i.e., AllDifferent,
 * Table (positive and negative), Element and Cumulative, plus NotEqual for a
 * kLinear with the kNE operator, which is not a single row.
 *
 * These are module-local classes, in the CSPBlock module and not in the
 * SMS++ core: whether the abstract globals belong in the core, as the CP
 * analogue of RowConstraint, or in the module, is left open in section 3 of
 * the design note, and placing them here keeps phase 1 additive and away
 * from the core.
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

#ifndef __CSPGlobalConstraints
 #define __CSPGlobalConstraints
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <list>
#include <stdexcept>
#include <vector>

#include "ColVariable.h"
#include "Constraint.h"
#include "FRowConstraint.h"
#include "LinearFunction.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------------------- CLASSES ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup CSPGlobalConstraints_CLASSES Classes in CSPGlobalConstraints.h
 *  @{ */

class CSPBlock;  // forward declaration: the owner of the constraints

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS CSPConstraint --------------------------*/
/*--------------------------------------------------------------------------*/
/// base class of the abstract global constraints of a CSPBlock
/** The CSPConstraint class is the common base of the abstract global
 * constraints of a CSPBlock [see CSPBlock.h]. It is a concrete :Constraint
 * over a fixed set of integer ColVariable, the scope of the global, and it
 * implements the :Constraint interface by delegating to the physical
 * representation: compute() reads the values of the ColVariable of the
 * CSPBlock, hands them to CSPBlock::is_satisfied() for the constraint index
 * this stands for, and stores the answer, which feasible() returns. The
 * abstract verifier is thus the physical one, so the two never disagree.
 *
 * Each derived class differs only in its linear relaxation, lin_relaxation()
 * [see the derived classes], which is best effort: a cheap, standard
 * linearization where one is at hand, an empty list with a TODO otherwise.
 * The verifier is in the base, since is_satisfied() already dispatches on the
 * kind; a derived class exists so that a Solver can tell the globals apart
 * and so that each carries its own relaxation.
 *
 * The scope is given once, at set_data(), and is not meant to change: the
 * single pure virtual of ThinVarDepInterface that would change it,
 * remove_variable(), throws. */

class CSPConstraint : public Constraint
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*------------------------------ PUBLIC TYPES ------------------------------*/
 /// the type of the vector of (pointers to) the ColVariable of the scope

 using v_var_ptr = std::vector< ColVariable * >;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// virtualized concrete iterator over the "active" ColVariable
 /** A concrete class deriving from ThinVarDepInterface::v_iterator for
  * sifting through the "active" ColVariable of the constraint, i.e., the
  * distinct variables of its scope, held as a vector of pointers. */

 class v_iterator : public ThinVarDepInterface::v_iterator {
  public:

  explicit v_iterator( v_var_ptr::iterator itr ) : itr_( itr ) {}
  v_iterator * clone( void ) override { return( new v_iterator( itr_ ) ); }
  void operator++( void ) final { ++itr_; }
  reference operator*( void ) const final { return( **itr_ ); }
  pointer operator->( void ) const final { return( *itr_ ); }
  bool operator==( const ThinVarDepInterface::v_iterator & rhs ) const final {
   const auto * tmp =
   #ifdef NDEBUG
    static_cast< const CSPConstraint::v_iterator * >( & rhs );
    return( itr_ == tmp->itr_ );
   #else
    dynamic_cast< const CSPConstraint::v_iterator * >( & rhs );
    return( tmp ? itr_ == tmp->itr_ : false );
   #endif
   }
  bool operator!=( const ThinVarDepInterface::v_iterator & rhs ) const final {
   const auto * tmp =
   #ifdef NDEBUG
    static_cast< const CSPConstraint::v_iterator * >( & rhs );
    return( itr_ != tmp->itr_ );
   #else
    dynamic_cast< const CSPConstraint::v_iterator * >( & rhs );
    return( tmp ? itr_ != tmp->itr_ : true );
   #endif
   }

  private:

  v_var_ptr::iterator itr_;
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// virtualized concrete const_iterator over the "active" ColVariable

 class v_const_iterator : public ThinVarDepInterface::v_const_iterator {
  public:

  explicit v_const_iterator( v_var_ptr::const_iterator itr ) : itr_( itr ) {}
  v_const_iterator * clone( void ) override {
   return( new v_const_iterator( itr_ ) );
   }
  void operator++( void ) final { ++itr_; }
  reference operator*( void ) const final { return( **itr_ ); }
  pointer operator->( void ) const final { return( *itr_ ); }
  bool operator==( const ThinVarDepInterface::v_const_iterator & rhs )
   const final {
   const auto * tmp =
   #ifdef NDEBUG
    static_cast< const CSPConstraint::v_const_iterator * >( & rhs );
    return( itr_ == tmp->itr_ );
   #else
    dynamic_cast< const CSPConstraint::v_const_iterator * >( & rhs );
    return( tmp ? itr_ == tmp->itr_ : false );
   #endif
   }
  bool operator!=( const ThinVarDepInterface::v_const_iterator & rhs )
   const final {
   const auto * tmp =
   #ifdef NDEBUG
    static_cast< const CSPConstraint::v_const_iterator * >( & rhs );
    return( itr_ != tmp->itr_ );
   #else
    dynamic_cast< const CSPConstraint::v_const_iterator * >( & rhs );
    return( tmp ? itr_ != tmp->itr_ : true );
   #endif
   }

  private:

  v_var_ptr::const_iterator itr_;
  };

/*------------------------ CONSTRUCTOR & DESTRUCTOR ------------------------*/
 /// constructor: takes the Block the constraint belongs to

 explicit CSPConstraint( Block * fblock = nullptr ) : Constraint( fblock ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// destructor: unregisters from the ColVariable of the scope

 ~CSPConstraint() override { unregister_vars(); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// "rough destructor": does not unregister from the ColVariable
 /** As in the other :Constraint, clear() empties the scope without taking
  * the constraint out of the active list of its ColVariable, which is what
  * is wanted when the ColVariable are being destroyed anyway. */

 void clear( void ) override { f_var.clear(); }

/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
 /// sets the CSPBlock, the constraint index and the scope, and registers
 /** Sets the CSPBlock \p csp this constraint stands for, the index \p k of
  * the constraint in its physical catalogue and the distinct ColVariable of
  * its scope, in \p vars, registering the constraint as active on each of
  * them. It is meant to be called once, by the CSPBlock that builds the
  * abstract representation. */

 void set_data( CSPBlock * csp , Index k , v_var_ptr && vars ) {
  unregister_vars();
  f_csp = csp;
  f_k = k;
  f_var = std::move( vars );
  for( auto v : f_var )
   v->add_active( this );
  }

/*------------ METHODS DESCRIBING THE BEHAVIOR OF THE Constraint -----------*/
 /// evaluates the constraint through the physical verifier of the CSPBlock
 /** Reads the values of all the ColVariable of the CSPBlock, rounds them to
  * integers, and stores whether CSPBlock::is_satisfied() holds for the
  * constraint index this stands for; feasible() then returns it. Defined in
  * CSPBlock.cpp, where the CSPBlock is a complete type. */

 int compute( bool changedvars = true ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns true if the constraint was satisfied at the last compute()

 [[nodiscard]] bool feasible( void ) const override { return( f_feasible ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the linear relaxation of the constraint, as FRowConstraint
 /** Returns the linear relaxation of the constraint as a list of newly
  * allocated FRowConstraint over the ColVariable of the scope, which become
  * property of the caller, i.e., the caller deletes them. The base returns
  * an empty list, i.e., no linear view; the derived classes override it
  * where a cheap, valid linearization is at hand [see them]. Defined in
  * CSPBlock.cpp. */

 [[nodiscard]] virtual std::list< FRowConstraint * > lin_relaxation( void )
  const { return( std::list< FRowConstraint * >() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// returns the index of the constraint in the physical catalogue

 [[nodiscard]] Index get_constraint_index( void ) const { return( f_k ); }

/*------- METHODS FOR HANDLING THE "ACTIVE" Variable IN THE Constraint -----*/

 [[nodiscard]] Index get_num_active_var( void ) const override {
  return( f_var.size() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 Index is_active( const Variable * var ) const override {
  const auto it = std::find( f_var.begin() , f_var.end() , var );
  return( it != f_var.end() ? Index( it - f_var.begin() ) : Inf< Index >() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 [[nodiscard]] Variable * get_active_var( Index i ) const override {
  return( i < f_var.size() ? f_var[ i ] : nullptr );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 v_iterator * v_begin( void ) override {
  return( new CSPConstraint::v_iterator( f_var.begin() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 [[nodiscard]] v_const_iterator * v_begin( void ) const override {
  return( new CSPConstraint::v_const_iterator( f_var.begin() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 v_iterator * v_end( void ) override {
  return( new CSPConstraint::v_iterator( f_var.end() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 [[nodiscard]] v_const_iterator * v_end( void ) const override {
  return( new CSPConstraint::v_const_iterator( f_var.end() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the scope is fixed: removing a variable is not supported, it throws

 void remove_variable( Index i , ModParam issueMod = eModBlck ) override {
  throw( std::logic_error( "CSPConstraint::remove_variable: the scope of an "
			   "abstract global constraint is fixed" ) );
  }

/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

/*--------------------------- PROTECTED METHODS ----------------------------*/
 /// takes the constraint out of the active list of each of its ColVariable

 void unregister_vars( void ) {
  for( auto v : f_var )
   v->remove_active( this );
  f_var.clear();
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// prints the constraint

 void print( std::ostream & output ) const override {
  output << "CSPConstraint [" << this << "] of Block [" << get_Block()
	 << "], index " << f_k << ", " << f_var.size() << " active variables"
	 << std::endl;
  }

/*---------------------------- PROTECTED FIELDS ----------------------------*/

 CSPBlock * f_csp = nullptr;  ///< the CSPBlock this constraint stands for

 Index f_k = 0;               ///< the index of the constraint in the block

 v_var_ptr f_var;             ///< the distinct ColVariable of the scope

 bool f_feasible = false;     ///< the answer of the last compute()

/*--------------------------------------------------------------------------*/

 };  // end( class( CSPConstraint ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS AllDifferentConstraint ---------------------*/
/*--------------------------------------------------------------------------*/
/// the abstract form of a kAllDifferent constraint
/** The variables of the scope must take pairwise different values. The
 * verifier is the base one; the linear relaxation is the cheap, valid "sum"
 * bound (Williams and Yan): if the k distinct variables must take distinct
 * values, their sum is at least the sum of the k smallest and at most the sum
 * of the k largest distinct values that any of them may take, which is read
 * off the current domains. This is a single FRowConstraint lb <= sum x_i <=
 * ub; it is weak but valid, and far cheaper than the pairwise big-M or the
 * clique inequalities. */

class AllDifferentConstraint : public CSPConstraint
{
 public:

 explicit AllDifferentConstraint( Block * fblock = nullptr )
  : CSPConstraint( fblock ) {}

 ~AllDifferentConstraint() override = default;

 [[nodiscard]] std::list< FRowConstraint * > lin_relaxation( void )
  const override;

 };  // end( class( AllDifferentConstraint ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS TableConstraint -------------------------*/
/*--------------------------------------------------------------------------*/
/// the abstract form of a kExtension or kConflicts constraint
/** The values of the scope must be one of the tuples (kExtension) or none of
 * them (kConflicts). The verifier is the base one. The linear relaxation of
 * a table is its convex hull, which is not cheap; a prototype weak view would
 * be the per-column value bounds (for a positive table) or nothing (for a
 * negative one). TODO: emit the per-column bounds for kExtension and a
 * convex-hull or big-M view in general; for now the relaxation is empty. */

class TableConstraint : public CSPConstraint
{
 public:

 explicit TableConstraint( Block * fblock = nullptr )
  : CSPConstraint( fblock ) {}

 ~TableConstraint() override = default;

 // TODO: lin_relaxation() left as the base (empty); a positive table admits
 //       the cheap per-column min..max bounds, a general one its convex hull.

 };  // end( class( TableConstraint ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS ElementConstraint ------------------------*/
/*--------------------------------------------------------------------------*/
/// the abstract form of a kElement constraint
/** With the scope ( i , z , x_0 , ... , x_{k - 1} ) and offset o, the index
 * variable i selects x_{i - o} and forces it equal to z. The verifier is the
 * base one. A valid linear relaxation needs an indicator or big-M per entry
 * (z - x_j is 0 when i = o + j), which is not cheap. TODO: emit the big-M
 * form; for now the relaxation is empty. */

class ElementConstraint : public CSPConstraint
{
 public:

 explicit ElementConstraint( Block * fblock = nullptr )
  : CSPConstraint( fblock ) {}

 ~ElementConstraint() override = default;

 // TODO: lin_relaxation() left as the base (empty); the standard view is the
 //       big-M/indicator linking z to x_{i - o}.

 };  // end( class( ElementConstraint ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS CumulativeConstraint ----------------------*/
/*--------------------------------------------------------------------------*/
/// the abstract form of a kCumulative constraint
/** At every time t the resource used by the tasks running at t is at most the
 * capacity. The verifier is the base one. A valid linear relaxation is the
 * time-indexed or the energetic-reasoning formulation, neither of which is
 * cheap and both of which need auxiliary variables the prototype does not
 * add. TODO: emit the time-indexed relaxation; for now the relaxation is
 * empty. */

class CumulativeConstraint : public CSPConstraint
{
 public:

 explicit CumulativeConstraint( Block * fblock = nullptr )
  : CSPConstraint( fblock ) {}

 ~CumulativeConstraint() override = default;

 // TODO: lin_relaxation() left as the base (empty); the standard views are
 //       the time-indexed and the energetic one, both needing extra
 //       variables.

 };  // end( class( CumulativeConstraint ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS NotEqualConstraint -----------------------*/
/*--------------------------------------------------------------------------*/
/// the abstract form of a kLinear constraint with the kNE operator
/** A kLinear constraint with the kLE or the kEQ operator is a single
 * FRowConstraint and is built directly by the CSPBlock, without a class here;
 * one with the kNE operator, \f$\sum_i a_i v_i \neq b\f$, is not a single row
 * (its feasible region is not convex), so it is kept as an abstract verifier
 * here, with no linear view. The verifier is the base one, over is_satisfied.
 * TODO: emit the two-sided big-M disjunction (sum <= b - 1 OR sum >= b + 1)
 * with an auxiliary binary, when the abstract objective and the auxiliaries
 * are available; for now the relaxation is empty. */

class NotEqualConstraint : public CSPConstraint
{
 public:

 explicit NotEqualConstraint( Block * fblock = nullptr )
  : CSPConstraint( fblock ) {}

 ~NotEqualConstraint() override = default;

 // TODO: lin_relaxation() left as the base (empty); the kNE is a non-convex
 //       disjunction, linearized only with an auxiliary binary and big-M.

 };  // end( class( NotEqualConstraint ) )

/** @} end( group( CSPGlobalConstraints_CLASSES ) ) */

/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* CSPGlobalConstraints.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File CSPGlobalConstraints.h --------------------*/
/*--------------------------------------------------------------------------*/

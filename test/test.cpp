/*--------------------------------------------------------------------------*/
/*------------------------------ File test.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Tests for CSPBlock on hand-made cases: the factory, the loader from
 * vectors and from FlatZinc (and what each refuses), the check of each kind
 * of constraint at its limits, the feasibility of an assignment and its
 * objective, the round trips through netCDF and through the writer, the
 * physical changes and what each Modification says it changes [see
 * Modification::changes()], and the CSPBlockChange with its undo.
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

#include <cstdio>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CSPBlock.h"
#include "FakeSolver.h"

// the checks hold in every build type, Release included
#undef NDEBUG
#include <cassert>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using v_int = CSPBlock::v_int;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// loads a CSPBlock out of a FlatZinc string

static void load_string( CSPBlock & b , const std::string & text )
{
 std::istringstream in( text );
 b.load( in , 'F' );
 }

/*--------------------------------------------------------------------------*/
/// true if the call throws std::invalid_argument

static bool throws( const std::function< void( void ) > & call )
{
 try {
  call();
  }
 catch( std::invalid_argument & ) {
  return( true );
  }
 return( false );
 }

/*--------------------------------------------------------------------------*/
/// true if the two CSPBlock have the same physical representation

static bool same( const CSPBlock & a , const CSPBlock & b )
{
 return( ( a.get_universes() == b.get_universes() ) &&
	 ( a.get_domain_mask() == b.get_domain_mask() ) &&
	 ( a.get_costs() == b.get_costs() ) &&
	 ( a.get_variable_groups() == b.get_variable_groups() ) &&
	 ( a.get_kinds() == b.get_kinds() ) &&
	 ( a.get_scopes() == b.get_scopes() ) &&
	 ( a.get_parameters() == b.get_parameters() ) &&
	 ( a.get_weights() == b.get_weights() ) &&
	 ( a.get_active() == b.get_active() ) );
 }

/*--------------------------------------------------------------------------*/
/// the CSPBlock written and read back through netCDF and through FlatZinc

static void round_trips( const CSPBlock & b )
{
 const char * file = "CSPBlock_unit_test.nc4";
 {
  netCDF::NcFile f( file , netCDF::NcFile::replace );
  auto g = f.addGroup( "Block" );
  b.serialize( g );
  }
 {
  netCDF::NcFile f( file , netCDF::NcFile::read );
  auto d = Block::new_Block( f.getGroup( "Block" ) );
  auto c = dynamic_cast< CSPBlock * >( d );
  assert( c && same( *c , b ) );
  delete d;
  }
 std::remove( file );

 std::ostringstream full;
 b.print( full , 'C' );
 CSPBlock p;
 load_string( p , full.str() );
 assert( same( p , b ) );
 }

/*--------------------------------------------------------------------------*/

static void test_factory( void )
{
 auto block = Block::new_Block( "CSPBlock" );
 assert( block && dynamic_cast< CSPBlock * >( block ) );
 delete block;
 auto change = Change::new_Change( "CSPBlockChange" );
 assert( change && dynamic_cast< CSPBlockChange * >( change ) );
 delete change;
 }

/*--------------------------------------------------------------------------*/

static void test_load( void )
{
 CSPBlock b;

 // the universes are sorted and their repeated values kept once
 b.load( { { 3 , 1 , 3 } , { -2 } } , { CSPBlock::kAllDifferent } ,
	 { { 0 , 1 } } , { {} } , { 5 } );
 assert( ( b.get_number_variables() == 2 ) &&
	 ( b.get_number_values() == 3 ) );
 assert( ( v_int( b.get_universe( 0 ).begin() , b.get_universe( 0 ).end() )
	   == v_int{ 1 , 3 } ) );
 assert( ( b.get_domain( 0 ) == v_int{ 1 , 3 } ) && ( ! b.is_hard( 0 ) ) );
 assert( ( b.value_index( 0 , 3 ) == 1 ) && ( b.value_index( 1 , -2 ) == 2 ) &&
	 ( b.value_index( 0 , 2 ) == Inf< Block::Index >() ) &&
	 ( b.value_index( 2 , 1 ) == Inf< Block::Index >() ) );

 // whatever is wrong throws and leaves the CSPBlock as it was
 const CSPBlock::v_v_int u = { { 0 , 1 } , { 0 , 1 } };
 const auto bad = [ & ]( int kind , v_int scope , v_int param ,
			 double w = Inf< double >() ) {
  return( throws( [ & ]() {
   b.load( CSPBlock::v_v_int( u ) , { kind } , { scope } , { param } ,
	   { w } ); } ) );
  };
 assert( throws( [ & ]() { b.load( { { 1 } , {} } , {} , {} , {} ); } ) );
 assert( bad( CSPBlock::kAllDifferent , { 0 , 2 } , {} ) );
 assert( bad( CSPBlock::kAllDifferent , { -1 } , {} ) );
 assert( bad( CSPBlock::kAllDifferent , { 0 } , { 1 } ) );
 assert( bad( CSPBlock::kNumberKinds , { 0 } , {} ) );
 assert( bad( CSPBlock::kExtension , {} , {} ) );
 assert( bad( CSPBlock::kExtension , { 0 , 1 } , { 1 , 0 , 1 } ) );
 assert( bad( CSPBlock::kLinear , { 0 , 1 } , { 1 , 1 , CSPBlock::kLE } ) );
 assert( bad( CSPBlock::kLinear , { 0 } , { 1 , 7 , 0 } ) );
 assert( bad( CSPBlock::kElement , { 0 , 1 } , { 0 } ) );
 assert( bad( CSPBlock::kCumulative , { 0 } , { 1 , -1 , 1 } ) );
 assert( bad( CSPBlock::kCumulative , { 0 } , { 1 , 1 } ) );
 assert( bad( CSPBlock::kAllDifferent , { 0 } , {} , -1 ) );
 assert( bad( CSPBlock::kAllDifferent , { 0 } , {} ,
	      std::numeric_limits< double >::quiet_NaN() ) );
 assert( ( b.get_number_variables() == 2 ) && ( b.get_weights()[ 0 ] == 5 ) );

 // the empty scopes allowed: linear, alldifferent and cumulative
 b.load( CSPBlock::v_v_int( u ) ,
	 { CSPBlock::kLinear , CSPBlock::kAllDifferent ,
	   CSPBlock::kCumulative } , { {} , {} , {} } ,
	 { { CSPBlock::kLE , -1 } , {} , { 0 } } );
 const v_int a = { 0 , 0 };
 assert( ( ! b.is_satisfied( 0 , a ) ) && b.is_satisfied( 1 , a ) &&
	 b.is_satisfied( 2 , a ) && b.all_hard() );
 round_trips( b );

 // no variable at all
 b.load( {} , {} , {} , {} );
 assert( b.is_feasible( v_int{} ) );
 round_trips( b );
 }

/*--------------------------------------------------------------------------*/

static void test_flatzinc( void )
{
 CSPBlock b;
 load_string( b ,
  "% a comment\n"
  "predicate my_pred(array [int] of var int: x);\n"
  "int: n = 3;\n"
  "bool: t = true;\n"
  "array [1..3] of int: c = [2, -1, 1];\n"
  "set of int: S = {2, 4};\n"
  "var 1..3: a :: output_var;\n"
  "var {1, 5, 7}: b2 :: var_is_introduced;\n"
  "var bool: f;\n"
  "var 0..9: k = 4;\n"
  "var 1..2: al = a;\n"
  "array [1..3] of var int: xs :: output_array([1..3]) = [a, b2, 5];\n"
  "constraint int_lin_le(c, [a, b2, k], 10);\n"
  "constraint int_lin_eq([1, 1], [a, 3], n) :: domain;\n"
  "constraint int_lt(a, 2);\n"
  "constraint int_ge(a, b2);\n"
  "constraint bool2int(f, a);\n"
  "constraint fzn_all_different_int(xs) :: weight(2.5);\n"
  "constraint array_var_int_element(a, xs, b2) :: inactive;\n"
  "constraint array_int_element(a, [10, 20, 30], 20);\n"
  "constraint set_in(k, S);\n"
  "constraint table_int([a, f], [1, 0, 2, 1]);\n"
  "constraint fzn_cumulative([a, b2], [2, 1], [1, 1], 1);\n"
  "solve :: int_search(xs, input_order, indomain_min, complete) "
  "maximize b2;\n" );

 // a, b2, f, k, and the variables of the values 5 and 20 (al is a, and
 // the values in the linear constraints go in their right-hand sides)
 assert( b.get_number_variables() == 6 );
 assert( ( b.get_domain( 0 ) == v_int{ 1 , 2 } ) &&
	 ( b.get_universe( 0 ).size() == 3 ) );
 assert( ( b.get_domain( 3 ) == v_int{ 4 } ) &&
	 ( b.get_universe( 3 ).size() == 10 ) );
 assert( ( b.get_domain( 4 ) == v_int{ 5 } ) );
 assert( b.get_number_constraints() == 11 );

 const auto & k = b.get_kinds();
 assert( ( k == v_int{ CSPBlock::kLinear , CSPBlock::kLinear ,
		       CSPBlock::kLinear , CSPBlock::kLinear ,
		       CSPBlock::kLinear , CSPBlock::kAllDifferent ,
		       CSPBlock::kElement , CSPBlock::kExtension ,
		       CSPBlock::kExtension , CSPBlock::kExtension ,
		       CSPBlock::kCumulative } ) );
 // the value 3 among the terms goes in the right-hand side
 const auto p1 = b.get_parameters()[ 1 ];
 assert( ( b.get_scopes()[ 1 ].size() == 1 ) &&
	 ( v_int( p1.begin() , p1.end() ) ==
	   v_int{ 1 , CSPBlock::kEQ , 0 } ) );
 // a < 2 is a - 2 <= -1, i.e., a <= 1
 const auto p2 = b.get_parameters()[ 2 ];
 assert( ( v_int( p2.begin() , p2.end() ) ==
	   v_int{ 1 , CSPBlock::kLE , 1 } ) );
 assert( ( b.get_weights()[ 5 ] == 2.5 ) && ( ! b.is_active( 6 ) ) &&
	 b.is_active( 5 ) );

 // maximize b2: the value d costs - d
 assert( ( b.get_costs()[ b.value_index( 1 , 7 ) ] == -7 ) &&
	 ( b.get_costs()[ b.value_index( 0 , 1 ) ] == 0 ) );

 // a = 1, b2 = 1, f = 1: a + 3 = 3 fails, as do the alldifferent and
 // the table, the others hold
 v_int x = { 1 , 1 , 1 , 4 , 5 , 20 };
 assert( ( ! b.is_satisfied( 1 , x ) ) && b.is_satisfied( 0 , x ) &&
	 b.is_satisfied( 2 , x ) && b.is_satisfied( 3 , x ) &&
	 b.is_satisfied( 4 , x ) );
 assert( ( ! b.is_satisfied( 5 , x ) ) && b.is_satisfied( 8 , x ) &&
	 ( ! b.is_satisfied( 9 , x ) ) );
 assert( ! b.is_feasible( x ) );
 round_trips( b );

 // what is not supported, or wrong
 const auto bad = [ & ]( const std::string & t ) {
  return( throws( [ & ]() { load_string( b , t ); } ) );
  };
 assert( bad( "var int: x;\nsolve satisfy;\n" ) );
 assert( bad( "var float: x;\nsolve satisfy;\n" ) );
 assert( bad( "var 1..2: x;\nconstraint int_times(x, x, x);\n" ) );
 assert( bad( "var 1..2: x;\nconstraint int_le(x, y);\n" ) );
 assert( bad( "var 1..2: x;\nvar 1..2: x;\n" ) );
 assert( bad( "var 1..2: x\n" ) );
 assert( bad( "var 1..2: x :: cost([1]);\n" ) );
 assert( bad( "var 1..2: x :: universe({1});\n" ) );
 assert( bad( "var 1..2: x;\nconstraint fzn_table_int([x, x], [1]);\n" ) );
 assert( bad( "var 1..2: x;\nconstraint int_le(x, 1) :: weight(-1);\n" ) );
 assert( bad( "var 1..2: x;\nconstraint int_le(x, 9999999999);\n" ) );
 assert( bad( "var 1..99999999: x;\n" ) );
 assert( bad( "var 1..2: x;\nsolve optimize x;\n" ) );
 assert( bad( "var 1..2: x @ 1;\n" ) );
 assert( b.get_number_variables() == 6 );  // left as it was
 assert( throws( [ & ]() {
  std::istringstream in( "" );
  b.load( in , 'D' ); } ) );
 }

/*--------------------------------------------------------------------------*/

static void test_constraints( void )
{
 CSPBlock b;
 // five variables over 0..9
 b.load( CSPBlock::v_v_int( 5 , { 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 } ) ,
	 { CSPBlock::kExtension , CSPBlock::kConflicts , CSPBlock::kLinear ,
	   CSPBlock::kLinear , CSPBlock::kLinear , CSPBlock::kAllDifferent ,
	   CSPBlock::kElement , CSPBlock::kElement , CSPBlock::kCumulative ,
	   CSPBlock::kCumulative , CSPBlock::kExtension } ,
	 { { 0 , 1 } , { 0 , 1 } , { 0 , 1 , 0 } , { 0 , 1 } , { 0 , 1 } ,
	   { 0 , 1 , 0 } , { 0 , 1 , 2 , 3 } , { 0 , 1 , 2 , 3 } ,
	   { 0 , 1 , 2 } , { 0 , 1 } , { 0 } } ,
	 { { 1 , 2 , 3 , 4 } , { 1 , 2 } ,
	   { 2 , -1 , 1 , CSPBlock::kLE , 6 } ,  // 3 a_0 - a_1 <= 6
	   { 1 , 1 , CSPBlock::kEQ , 5 } , { 1 , 1 , CSPBlock::kNE , 5 } ,
	   {} , { 1 } , { -2 } ,
	   { 2 , 2 , 0 , 1 , 1 , 5 , 1 } ,   // the third task uses nothing
	   { 0 , 3 , 3 , 3 , 3 } , {} } );    // an empty table
 // the values not given are 0
 const auto sat = [ & ]( Block::Index k , v_int a ) {
  a.resize( 5 , 0 );
  return( b.is_satisfied( k , a ) );
  };
 // the tables
 assert( sat( 0 , { 3 , 4 } ) && ( ! sat( 0 , { 3 , 2 } ) ) );
 assert( ( ! sat( 1 , { 1 , 2 } ) ) && sat( 1 , { 2 , 1 } ) );
 assert( ! sat( 10 , { 0 } ) );
 // linear, a variable appearing twice, at the limit and beyond
 assert( sat( 2 , { 2 } ) && ( ! sat( 2 , { 3 , 2 } ) ) );
 assert( sat( 3 , { 2 , 3 } ) && ( ! sat( 3 , { 2 , 2 } ) ) );
 assert( ( ! sat( 4 , { 2 , 3 } ) ) && sat( 4 , { 2 , 2 } ) );
 // alldifferent with a variable twice is never satisfied
 assert( ! sat( 5 , { 1 , 2 } ) );
 // element of offset 1: x_{i - 1} = z, i in 1..2
 assert( sat( 6 , { 1 , 7 , 7 , 3 } ) && sat( 6 , { 2 , 3 , 7 , 3 } ) );
 assert( ( ! sat( 6 , { 2 , 7 , 7 , 3 } ) ) &&
	 ( ! sat( 6 , { 0 , 7 , 7 , 3 } ) ) &&
	 ( ! sat( 6 , { 3 , 3 , 7 , 3 } ) ) );
 // element of offset -2: i - o is at least 2 for i in 0..9, out of range
 assert( ! sat( 7 , { 0 , 3 , 7 , 3 } ) );
 // cumulative: tasks of 2 units, usage 1 and 1, capacity 1; one task ends
 // where the next starts, which is allowed; the third uses nothing
 assert( sat( 8 , { 0 , 2 , 1 } ) && sat( 8 , { 2 } ) );
 assert( ! sat( 8 , { 0 , 1 } ) );
 // a task of duration 0 uses nothing, the other alone fits exactly
 assert( sat( 9 , { 4 , 4 } ) );
 // the assignment must have one value per variable
 assert( throws( [ & ]() { (void) b.is_satisfied( 0 , v_int{ 1 } ); } ) );
 assert( throws( [ & ]() { (void) b.is_satisfied( 11 , v_int( 5 , 0 ) ); } ) );
 }

/*--------------------------------------------------------------------------*/

static void test_feasible( void )
{
 CSPBlock b;
 load_string( b ,
  "var 1..3: x :: universe(0..3) :: cost([0.5, 1, 2, 3]);\n"
  "var 1..3: y;\n"
  "constraint int_ne(x, y);\n"
  "constraint int_le(x, y) :: weight(4);\n"
  "constraint int_eq(x, 2) :: inactive;\n"
  "solve satisfy;\n" );
 assert( b.has_costs() && ( ! b.all_hard() ) );
 assert( b.is_feasible( v_int{ 3 , 1 } ) );
 assert( ( b.get_violated_weight( v_int{ 3 , 1 } ) == 4 ) &&
	 ( b.get_objective_value( v_int{ 3 , 1 } ) == 7 ) );
 assert( ( b.get_objective_value( v_int{ 1 , 3 } ) == 1 ) );
 assert( ! b.is_feasible( v_int{ 2 , 2 } ) );     // a hard one
 assert( ! b.is_feasible( v_int{ 0 , 1 } ) );     // 0 out of the domain
 assert( b.get_objective_value( v_int{ 0 , 1 } ) == 0.5 );
 assert( throws( [ & ]() {
  (void) b.get_objective_value( v_int{ 9 , 1 } ); } ) );
 assert( throws( [ & ]() { (void) b.is_feasible( v_int{ 1 } ); } ) );
 round_trips( b );
 }

/*--------------------------------------------------------------------------*/

static void test_round_trips( void )
{
 // every kind, both operators, an offset other than 1, groups, a domain
 // emptied, soft and inactive constraints, costs of no simple form
 CSPBlock b;
 b.load( { { -3 , 0 , 7 } , { 1 , 2 } , { 0 , 1 , 2 } , { 5 } } ,
	 { CSPBlock::kExtension , CSPBlock::kConflicts , CSPBlock::kLinear ,
	   CSPBlock::kLinear , CSPBlock::kAllDifferent , CSPBlock::kElement ,
	   CSPBlock::kElement , CSPBlock::kCumulative , CSPBlock::kLinear } ,
	 { { 0 , 1 } , { 1 , 2 } , { 0 , 1 , 1 } , {} , { 1 , 2 , 3 } ,
	   { 2 , 0 , 1 , 3 } , { 2 , 0 , 1 , 3 } , { 0 , 2 } , { 3 } } ,
	 { { -3 , 1 , 7 , 2 } , { 1 , 1 } ,
	   { 1 , -2 , 3 , CSPBlock::kEQ , 4 } ,
	   { CSPBlock::kLE , 0 } , {} , { 1 } , { 0 } , { 1 , 2 , 1 , 1 , 1 } ,
	   { 1 , CSPBlock::kNE , 5 } } ,
	 { Inf< double >() , 3 , 0 , Inf< double >() , 0.1 , Inf< double >() ,
	   7 , Inf< double >() , 1e-300 } );
 b.set_variable_groups( { 0 , 1 , 1 , 0 } );
 b.chg_domains( Block::Subset{ 0 , 3 } , { { 7 , -3 } , {} } );
 std::vector< double > c = { 1.5 , 0 , -1 , 0 , 0 , 0 , 0 , 0 , 1.0 / 3 };
 b.chg_costs( c , Block::Range( 0 , 9 ) );
 b.chg_active( false , Block::Subset{ 2 , 5 } );
 assert( ( b.get_domain( 3 ).empty() ) &&
	 ( b.get_domain( 0 ) == v_int{ -3 , 7 } ) );
 round_trips( b );

 // the costs of minimize and of maximize of a single variable are written
 // as the solve item, and read back the same
 for( double s : { 1.0 , -1.0 } ) {
  CSPBlock m;
  m.load( { { 1 , 2 } , { -4 , 0 , 6 } } , {} , {} , {} );
  std::vector< double > mc = { 0 , 0 , -4 * s , 0 , 6 * s };
  m.chg_costs( mc , Block::Range( 0 , 5 ) );
  std::ostringstream o;
  m.print( o , 'C' );
  assert( o.str().find( s > 0 ? "solve minimize x1;" : "solve maximize x1;" )
	  != std::string::npos );
  assert( o.str().find( "cost(" ) == std::string::npos );
  round_trips( m );
  }

 // the summary line
 std::ostringstream o;
 b.print( o );
 assert( o.str() == "CSPBlock with 4 variables and 9 constraints\n" );
 }

/*--------------------------------------------------------------------------*/

static void test_modifications( void )
{
 using M = Modification;
 const auto inf = Inf< double >();
 CSPBlock b;
 load_string( b ,
  "var 1..3: x;\n"
  "var 1..3: y;\n"
  "constraint int_ne(x, y);\n"
  "constraint int_le(x, y) :: weight(4);\n"
  "solve satisfy;\n" );
 FakeSolver fs;
 b.register_Solver( & fs );

 // the call issues exactly n CSPBlockMod, of type t, the first saying it
 // changes c, and none if n is 0
 const auto mods = [ & ]( const std::function< void( void ) > & call ,
			  unsigned n , int t = 0 , M::ModConcern c = 0 ) {
  fs.lock_Modification_list();
  fs.get_Modification_list().clear();
  fs.unlock_Modification_list();
  call();
  std::vector< std::shared_ptr< CSPBlockMod > > got;
  fs.lock_Modification_list();
  for( auto & mod : fs.get_Modification_list() )
   if( auto sm = std::dynamic_pointer_cast< CSPBlockMod >( mod ) )
    got.push_back( sm );
  fs.get_Modification_list().clear();
  fs.unlock_Modification_list();
  assert( got.size() == n );
  if( n ) {
   assert( got.front()->type() == t );
   assert( got.front()->changes() == c );
   }
  return( got );
  };

 // domains: shrink, grow, both, nothing
 mods( [ & ]() { b.chg_domains( Block::Subset{ 1 } , { { 1 , 2 } } ); } , 1 ,
       CSPBlockMod::eChgDomain , M::eModPhys | M::eRegnShrink );
 mods( [ & ]() { b.chg_domains( CSPBlock::v_v_int{ { 1 , 2 , 3 } } ,
					 Block::Range( 1 , 2 ) ); } ,
       1 , CSPBlockMod::eChgDomain , M::eModPhys | M::eRegnGrow );
 auto got = mods( [ & ]() {
  b.chg_domains( Block::Subset{ 1 , 0 } , { { 3 } , { 1 , 2 , 3 } } ); } , 1 ,
  CSPBlockMod::eChgDomain , M::eModPhys | M::eRegnShrink );
 assert( ( std::static_pointer_cast< CSPBlockSbstMod >( got.front() )->nms()
	   == Block::Subset{ 1 } ) );
 mods( [ & ]() { b.chg_domains( Block::Subset{ 1 } , { { 3 , 3 } } ); } , 0 );
 mods( [ & ]() {
  b.chg_domains( Block::Subset{ 0 , 1 } , { { 1 } , { 1 , 2 , 3 } } ); } , 1 ,
  CSPBlockMod::eChgDomain , M::eModPhys | M::eRegnShrink | M::eRegnGrow );
 // wrong ones change nothing
 assert( throws( [ & ]() {
  b.chg_domains( Block::Subset{ 0 } , { { 4 } } ); } ) );
 assert( throws( [ & ]() {
  b.chg_domains( Block::Subset{ 0 , 0 } , { {} , {} } ); } ) );
 assert( throws( [ & ]() {
  b.chg_domains( Block::Subset{ 2 } , { {} } ); } ) );
 assert( b.get_domain( 0 ) == v_int{ 1 } );

 // values: new ones join the universe and the domain, and the slots after
 // them move with their costs
 std::vector< double > c = { 9 };
 b.chg_costs( c , Block::Subset{ 3 } );  // y = 1
 mods( [ & ]() { b.add_values( { 0 } , { { 0 , 2 , 5 } } ); } , 1 ,
       CSPBlockMod::eAddValues , M::eModPhys | M::eRegnGrow );
 assert( ( b.get_domain( 0 ) == v_int{ 0 , 1 , 5 } ) &&
	 ( b.get_universe( 0 ).size() == 5 ) &&
	 ( b.get_number_values() == 8 ) );
 assert( b.get_costs()[ b.value_index( 1 , 1 ) ] == 9 );
 mods( [ & ]() { b.add_values( { 0 } , { { 1 , 5 } } ); } , 0 );

 // costs, by Range and by Subset, and nothing if they do not change
 c = { 1 , 2 };
 mods( [ & ]() { b.chg_costs( c , Block::Range( 0 , 2 ) ); } , 1 ,
       CSPBlockMod::eChgCost , M::ePhysObj );
 mods( [ & ]() { b.chg_costs( c , Block::Range( 0 , 2 ) ); } , 0 );
 c = { 3 , 4 };
 mods( [ & ]() { b.chg_costs( c , Block::Subset{ 7 , 6 } ); } , 1 ,
       CSPBlockMod::eChgCost , M::ePhysObj );
 assert( ( b.get_costs()[ 6 ] == 4 ) && ( b.get_costs()[ 7 ] == 3 ) );
 c = { inf };
 assert( throws( [ & ]() { b.chg_costs( c , Block::Range( 0 , 1 ) ); } ) );

 // weights
 const auto obj = M::ePhysObj;
 std::vector< double > w = { 2 };
 mods( [ & ]() { b.chg_weights( w , Block::Range( 1 , 2 ) ); } , 1 ,
       CSPBlockMod::eChgWeight , obj );
 w = { inf };
 mods( [ & ]() { b.chg_weights( w , Block::Subset{ 1 } ); } , 1 ,
       CSPBlockMod::eChgWeight , obj | M::eModCnsSide | M::eRegnShrink );
 w = { 1 , inf };
 mods( [ & ]() { b.chg_weights( w , Block::Range( 0 , 2 ) ); } , 1 ,
       CSPBlockMod::eChgWeight , obj | M::eModCnsSide | M::eRegnGrow );
 w = { -1 };
 assert( throws( [ & ]() { b.chg_weights( w , Block::Subset{ 0 } ); } ) );

 // activity: constraint 0 soft, 1 hard
 mods( [ & ]() { b.chg_active( false , Block::Range( 0 , 2 ) ); } , 1 ,
       CSPBlockMod::eChgActive ,
       M::eModPhys | M::eModCnsSide | M::eRegnGrow | M::eObjDown );
 mods( [ & ]() { b.chg_active( false , Block::Subset{ 1 } ); } , 0 );
 mods( [ & ]() { b.chg_active( true , Block::Subset{ 1 } ); } , 1 ,
       CSPBlockMod::eChgActive ,
       M::eModPhys | M::eModCnsSide | M::eRegnShrink );

 // constraints and variables added
 mods( [ & ]() {
  b.add_constraints( { CSPBlock::kAllDifferent } , { { 0 , 1 } } , { {} } );
  } , 1 , CSPBlockMod::eAddConstraints ,
  M::eModPhys | M::eModCnsSet | M::eRegnShrink | M::eObjUp );
 assert( ( b.get_number_constraints() == 3 ) && b.is_active( 2 ) );
 assert( throws( [ & ]() {
  b.add_constraints( { CSPBlock::kAllDifferent } , { { 0 , 2 } } , { {} } );
  } ) );
 mods( [ & ]() { b.add_variables( { { 4 , 2 } } ); } , 1 ,
       CSPBlockMod::eAddVariables , M::eModPhys | M::eEffAny );
 assert( ( b.get_number_variables() == 3 ) &&
	 ( b.get_domain( 2 ) == v_int{ 2 , 4 } ) &&
	 ( b.get_costs().back() == 0 ) );
 assert( throws( [ & ]() { b.add_variables( { {} } ); } ) );
 assert( throws( [ & ]() { b.add_variables( { { 1 } } , { 0 } ); } ) );
 b.set_variable_groups( { 0 , 0 , 1 } );
 assert( throws( [ & ]() { b.add_variables( { { 1 } } ); } ) );
 b.add_variables( { { 1 } } , { 2 } );
 assert( b.get_variable_groups() == v_int( { 0 , 0 , 1 , 2 } ) );
 assert( throws( [ & ]() { b.set_variable_groups( { 0 } ); } ) );
 assert( throws( [ & ]() { b.set_variable_groups( { 0 , 0 , -1 , 0 } ); } ) );
 round_trips( b );

 b.unregister_Solver( & fs );
 }

/*--------------------------------------------------------------------------*/

static void test_changes( void )
{
 CSPBlock b;
 b.load( { { 1 , 2 , 3 , 4 } , { 1 , 2 , 3 , 4 } } , {} , {} , {} );

 // a change, then its undo, brings the domains back
 const auto check = [ & ]( CSPBlockChange & ch , v_int d0 , v_int d1 ) {
  const auto o0 = b.get_domain( 0 ) , o1 = b.get_domain( 1 );
  auto undo = ch.apply( & b , true );
  assert( ( b.get_domain( 0 ) == d0 ) && ( b.get_domain( 1 ) == d1 ) );
  assert( undo && static_cast< CSPBlockChange * >( undo )->type() ==
	  CSPBlockChange::eSetDomain );
  // the undo through netCDF
  const char * file = "CSPBlock_unit_test_change.nc4";
  {
   netCDF::NcFile f( file , netCDF::NcFile::replace );
   auto g = f.addGroup( "Change" );
   undo->serialize( g );
   }
  delete undo;
  {
   netCDF::NcFile f( file , netCDF::NcFile::read );
   std::unique_ptr< Change > u( Change::new_Change( f.getGroup( "Change" ) ) );
   assert( u );
   u->apply( & b );
   }
  std::remove( file );
  assert( ( b.get_domain( 0 ) == o0 ) && ( b.get_domain( 1 ) == o1 ) );
  ch.apply( & b );
  };

 CSPBlockChange a( CSPBlockChange::eSplitLE , { 0 , 1 } , { 3 , 1 } );
 check( a , { 1 , 2 , 3 } , { 1 } );
 CSPBlockChange r( CSPBlockChange::eRemove , { 0 , 1 } , { 2 , 3 } );
 check( r , { 1 , 3 } , { 1 } );
 CSPBlockChange s( CSPBlockChange::eSplitGT , { 0 } , { 1 } );
 check( s , { 3 } , { 1 } );
 CSPBlockChange e( CSPBlockChange::eAssign , { 1 } , { 4 } );  // not there
 check( e , { 3 } , {} );
 CSPBlockChange d( CSPBlockChange::eSetDomain , { 1 , 0 } , { 2 , 4 } ,
		   { 2 , 0 } );
 check( d , {} , { 2 , 4 } );

 // what is wrong throws, and changes nothing
 CSPBlockChange w1( CSPBlockChange::eAssign , { 0 , 0 } , { 1 , 1 } );
 CSPBlockChange w2( CSPBlockChange::eAssign , { 2 } , { 1 } );
 CSPBlockChange w3( CSPBlockChange::eSetDomain , { 0 } , { 1 } , { 2 } );
 CSPBlockChange w4;
 for( auto ch : { & w1 , & w2 , & w3 , & w4 } )
  assert( throws( [ & ]() { delete ch->apply( & b , true ); } ) );
 assert( ( b.get_domain( 0 ).empty() ) &&
	 ( b.get_domain( 1 ) == v_int{ 2 , 4 } ) );
 }

/*--------------------------------------------------------------------------*/
/// generates the abstract representation and checks it against the physical
/** Builds the abstract Variable and Constraint of \p b, then, for each
 * assignment in \p cases, sets the values of the ColVariable and checks that
 * the abstract feasibility is_feasible( true ) agrees with the physical
 * is_feasible( std::span ); the domains of \p b must be contiguous and free
 * of holes, so that the box of each variable is exactly its domain and the
 * two representations coincide. */

static void check_abstract( CSPBlock & b , const std::vector< v_int > & cases )
{
 b.generate_abstract_variables();
 b.generate_abstract_constraints();
 assert( b.get_abstract_variables().size() == b.get_number_variables() );

 auto & var = b.get_abstract_variables();
 for( const auto & a : cases ) {
  assert( a.size() == var.size() );
  for( std::size_t i = 0 ; i < a.size() ; ++i )
   var[ i ].set_value( a[ i ] );
  const bool phys = b.is_feasible( std::span< const int >( a ) );
  const bool abst = b.is_feasible( true );
  assert( phys == abst );
  }
 }

/*--------------------------------------------------------------------------*/

static void test_abstract( void )
{
 // an alldifferent instance: 3 variables over 0..3, all pairwise different
 {
  CSPBlock b;
  b.load( CSPBlock::v_v_int( 3 , { 0 , 1 , 2 , 3 } ) ,
	  { CSPBlock::kAllDifferent } , { { 0 , 1 , 2 } } , { {} } );
  check_abstract( b , { { 0 , 1 , 2 } ,    // feasible
			{ 3 , 2 , 1 } ,    // feasible
			{ 0 , 0 , 1 } ,    // alldifferent violated
			{ 2 , 2 , 2 } ,    // alldifferent violated
			{ 0 , 1 , 3 } } ); // feasible
 }

 // a linear instance: 2 variables over 0..5, a <=, an == and a constraint
 // with a repeated variable (2 x - x = x <= 3), to exercise the aggregation
 {
  CSPBlock b;
  b.load( CSPBlock::v_v_int( 2 , { 0 , 1 , 2 , 3 , 4 , 5 } ) ,
	  { CSPBlock::kLinear , CSPBlock::kLinear , CSPBlock::kLinear } ,
	  { { 0 , 1 } , { 0 , 1 } , { 0 , 0 } } ,
	  { { 1 , 1 , CSPBlock::kLE , 5 } ,    // x + y <= 5
	    { 1 , -1 , CSPBlock::kEQ , 1 } ,   // x - y == 1
	    { 2 , -1 , CSPBlock::kLE , 3 } } );  // x <= 3
  check_abstract( b , { { 3 , 2 } ,    // feasible: 5<=5, 1==1, 3<=3
			{ 1 , 0 } ,    // feasible
			{ 4 , 2 } ,    // x + y = 6 > 5
			{ 2 , 2 } ,    // x - y = 0 != 1
			{ 4 , 3 } ,    // x = 4 > 3 (the repeated-variable one)
			{ 5 , 4 } } ); // x + y = 9 > 5
 }

 // a not-equal linear constraint, which is not a single row: x != 3
 {
  CSPBlock b;
  b.load( { { 0 , 1 , 2 , 3 , 4 } } , { CSPBlock::kLinear } , { { 0 } } ,
	  { { 1 , CSPBlock::kNE , 3 } } );
  check_abstract( b , { { 0 } , { 2 } , { 3 } , { 4 } } );
 }

 // an inactive hard constraint and a soft one are not in the abstract
 // representation, so an assignment violating only them is abstract-feasible
 // while it is still physically feasible (is_feasible ignores both)
 {
  CSPBlock b;
  b.load( CSPBlock::v_v_int( 2 , { 0 , 1 , 2 , 3 } ) ,
	  { CSPBlock::kAllDifferent , CSPBlock::kLinear } ,
	  { { 0 , 1 } , { 0 , 1 } } ,
	  { {} , { 1 , 1 , CSPBlock::kLE , 2 } } ,
	  { Inf< double >() , 4 } );         // the linear one is soft
  b.chg_active( false , Block::Subset{ 0 } );  // the alldifferent is inactive
  check_abstract( b , { { 3 , 3 } , { 0 , 0 } , { 2 , 3 } } );
 }

 // the linear relaxation of alldifferent: the "sum" bound, built on its own
 {
  CSPBlock b;
  b.load( CSPBlock::v_v_int( 3 , { 0 , 1 , 2 , 3 } ) ,
	  { CSPBlock::kAllDifferent } , { { 0 , 1 , 2 } } , { {} } );
  b.generate_abstract_variables();
  b.generate_abstract_constraints();
  auto & var = b.get_abstract_variables();

  AllDifferentConstraint ad( & b );
  ad.set_data( & b , 0 , { & var[ 0 ] , & var[ 1 ] , & var[ 2 ] } );
  auto rel = ad.lin_relaxation();
  assert( rel.size() == 1 );
  auto * fr = rel.front();
  // 3 distinct values out of {0,1,2,3}: sum in [0+1+2, 3+2+1] = [3, 6]
  assert( ( fr->get_lhs() == 3 ) && ( fr->get_rhs() == 6 ) );
  // the relaxation is valid: an alldifferent-feasible point lies in it, and it
  // rejects a point whose sum is out of the bound
  var[ 0 ].set_value( 0 ); var[ 1 ].set_value( 1 ); var[ 2 ].set_value( 2 );
  fr->compute();
  assert( fr->feasible() );
  var[ 0 ].set_value( 0 ); var[ 1 ].set_value( 0 ); var[ 2 ].set_value( 1 );
  fr->compute();
  assert( ! fr->feasible() );  // sum 1 < 3
  delete fr;
 }

 // a table, an element and a cumulative give no linear view in the prototype
 {
  CSPBlock b;
  b.load( CSPBlock::v_v_int( 3 , { 0 , 1 , 2 , 3 } ) , {} , {} , {} );
  b.generate_abstract_variables();
  auto & var = b.get_abstract_variables();
  TableConstraint tb( & b );
  tb.set_data( & b , 0 , { & var[ 0 ] , & var[ 1 ] } );
  ElementConstraint el( & b );
  el.set_data( & b , 0 , { & var[ 0 ] , & var[ 1 ] , & var[ 2 ] } );
  CumulativeConstraint cu( & b );
  cu.set_data( & b , 0 , { & var[ 0 ] , & var[ 1 ] } );
  assert( tb.lin_relaxation().empty() && el.lin_relaxation().empty() &&
	  cu.lin_relaxation().empty() );
 }
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 test_factory();
 test_load();
 test_flatzinc();
 test_constraints();
 test_feasible();
 test_round_trips();
 test_modifications();
 test_changes();
 test_abstract();

 std::cout << "CSPBlock: all tests passed" << std::endl;
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- End File test.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/

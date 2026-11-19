/*--------------------------------------------------------------------------*/
/*----------------------------- File battery.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Battery for CSPBlock on small instances of known families, each written
 * in FlatZinc as MiniZinc would flatten it: the n queens, the colouring of
 * the map of Australia with k colours, a sudoku, and the scheduling of
 * tasks on a resource of given capacity, minimizing the makespan. For each
 * instance, the CSPBlock read from FlatZinc is written and read back
 * through netCDF and through FlatZinc, and its solutions are enumerated by
 * backtracking on the check of the constraints [see
 * CSPBlock::is_satisfied()]: each must be feasible [see
 * CSPBlock::is_feasible()], and their number must be the one of the batch
 * file, or, for the scheduling, the one counted by a check of the resource
 * time by time that does not use CSPBlock, which also gives the optimal
 * makespan.
 *
 *     CSPBlock_battery_test <family> <size> <solutions | ->
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

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CSPBlock.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using v_int = CSPBlock::v_int;

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// a FlatZinc array of the names prefix<i> for the given i

static std::string names( const std::string & prefix , const v_int & i )
{
 std::string s = "[";
 for( std::size_t k = 0 ; k < i.size() ; ++k )
  s += ( k ? ", " : "" ) + prefix + std::to_string( i[ k ] );
 return( s + "]" );
 }

/*--------------------------------------------------------------------------*/
/// the n queens: q_i is the row of the queen of column i

static std::string queens( int n )
{
 std::ostringstream o;
 v_int all( n );
 for( int i = 0 ; i < n ; ++i ) {
  all[ i ] = i;
  o << "var 1.." << n << ": q" << i << " :: output_var;\n";
  }
 o << "array [1.." << n << "] of var int: q :: output_array([1.." << n
   << "]) = " << names( "q" , all ) << ";\n";
 o << "constraint fzn_all_different_int(q);\n";
 for( int i = 0 ; i < n ; ++i )
  for( int j = i + 1 ; j < n ; ++j ) {
   const auto x = names( "q" , { i , j } );
   o << "constraint int_lin_ne([1, -1], " << x << ", " << j - i << ");\n"
     << "constraint int_lin_ne([1, -1], " << x << ", " << i - j << ");\n";
   }
 o << "solve satisfy;\n";
 return( o.str() );
 }

/*--------------------------------------------------------------------------*/
/// the map of Australia with k colours

static std::string map_colouring( int k )
{
 static const char * region[] = { "wa" , "nt" , "sa" , "q" , "nsw" , "v" ,
				  "t" };
 static const int edge[][ 2 ] = { { 0 , 1 } , { 0 , 2 } , { 1 , 2 } ,
				  { 1 , 3 } , { 2 , 3 } , { 2 , 4 } ,
				  { 2 , 5 } , { 3 , 4 } , { 4 , 5 } };
 std::ostringstream o;
 for( auto r : region )
  o << "var 1.." << k << ": " << r << ";\n";
 for( auto & e : edge )
  o << "constraint int_ne(" << region[ e[ 0 ] ] << ", " << region[ e[ 1 ] ]
    << ");\n";
 o << "solve satisfy;\n";
 return( o.str() );
 }

/*--------------------------------------------------------------------------*/
/// a sudoku of side n (4 or 9), or, for n = 0, a 4 x 4 one with no solution

static std::string sudoku( int n )
{
 static const std::string p4 = "1...3..2...3..2.";
 static const std::string p0 = "1.1.............";
 static const std::string p9 =
  "53..7....6..195....98....6.8...6...34..8.3..17...2...6"
  ".6....28....419..5....8..79";
 const auto & p = ( n == 9 ) ? p9 : ( n == 4 ? p4 : p0 );
 if( n == 0 )
  n = 4;
 const int b = ( n == 9 ) ? 3 : 2;
 if( int( p.size() ) != n * n )
  throw( std::logic_error( "sudoku: wrong puzzle" ) );

 std::ostringstream o;
 for( int c = 0 ; c < n * n ; ++c ) {
  o << "var 1.." << n << ": c" << c;
  if( p[ c ] != '.' )
   o << " = " << p[ c ];  // a given: its value only in the domain
  o << ";\n";
  }
 for( int r = 0 ; r < n ; ++r ) {
  v_int row , col , box;
  for( int k = 0 ; k < n ; ++k ) {
   row.push_back( r * n + k );
   col.push_back( k * n + r );
   box.push_back( ( r / b * b + k / b ) * n + r % b * b + k % b );
   }
  for( const auto & x : { row , col , box } )
   o << "constraint fzn_all_different_int(" << names( "c" , x ) << ");\n";
  }
 o << "solve satisfy;\n";
 return( o.str() );
 }

/*--------------------------------------------------------------------------*/
/// the durations, usages and capacity of the scheduling of n tasks

static const int duration[] = { 2 , 2 , 1 , 3 , 1 };
static const int usage[] = { 1 , 1 , 2 , 1 , 2 };
static const int capacity = 2;
static const int horizon = 5;   ///< the latest start

/*--------------------------------------------------------------------------*/
/// n tasks on a resource, starts in 0..horizon, makespan m to minimize

static std::string scheduling( int n )
{
 std::ostringstream o;
 v_int all( n );
 int maxd = 0;
 for( int j = 0 ; j < n ; ++j ) {
  all[ j ] = j;
  maxd = std::max( maxd , duration[ j ] );
  o << "var 0.." << horizon << ": s" << j << ";\n";
  }
 o << "var 0.." << horizon + maxd << ": m;\n";
 o << "array [1.." << n << "] of int: d = [";
 for( int j = 0 ; j < n ; ++j )
  o << ( j ? ", " : "" ) << duration[ j ];
 o << "];\n";
 o << "constraint fzn_cumulative(" << names( "s" , all ) << ", d, [";
 for( int j = 0 ; j < n ; ++j )
  o << ( j ? ", " : "" ) << usage[ j ];
 o << "], " << capacity << ");\n";
 for( int j = 0 ; j < n ; ++j )  // s_j + d_j <= m
  o << "constraint int_lin_le([1, -1], [s" << j << ", m], "
    << - duration[ j ] << ");\n";
 o << "solve minimize m;\n";
 return( o.str() );
 }

/*--------------------------------------------------------------------------*/
/// the solutions of the scheduling of n tasks counted without CSPBlock,
/// checking the resource at each time, and the optimal makespan

static std::pair< long , int > scheduling_count( int n )
{
 int maxd = 0;
 for( int j = 0 ; j < n ; ++j )
  maxd = std::max( maxd , duration[ j ] );
 const int mmax = horizon + maxd;
 long count = 0;
 int best = mmax + 1;
 v_int s( n , 0 );
 for( ; ; ) {
  int end = 0;
  bool ok = true;
  for( int t = 0 ; ok && ( t < mmax ) ; ++t ) {
   int use = 0;
   for( int j = 0 ; j < n ; ++j )
    if( ( s[ j ] <= t ) && ( t < s[ j ] + duration[ j ] ) )
     use += usage[ j ];
   ok = ( use <= capacity );
   }
  for( int j = 0 ; j < n ; ++j )
   end = std::max( end , s[ j ] + duration[ j ] );
  if( ok && ( end <= mmax ) ) {
   count += mmax - end + 1;  // each m from end to mmax
   best = std::min( best , end );
   }
  int j = 0;
  while( ( j < n ) && ( ++s[ j ] > horizon ) )
   s[ j++ ] = 0;
  if( j == n )
   break;
  }
 return( std::make_pair( count , best ) );
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
/// the enumeration of the solutions of a CSPBlock by backtracking

class Enumerator {
 public:
 explicit Enumerator( const CSPBlock & b ) : f_b( b ) ,
  f_a( b.get_number_variables() , 0 ) ,
  f_check( b.get_number_variables() ) ,
  f_alldiff( b.get_number_variables() ) {
  // each active hard constraint is checked once its last variable has a
  // value, and an alldifferent one also at each of its variables
  for( Block::Index k = 0 ; k < b.get_number_constraints() ; ++k ) {
   if( ( ! b.is_active( k ) ) || ( ! b.is_hard( k ) ) )
    continue;
   const auto s = b.get_scopes()[ k ];
   if( s.empty() ) {
    f_empty.push_back( k );
    continue;
    }
   f_check[ *std::max_element( s.begin() , s.end() ) ].push_back( k );
   if( b.get_kinds()[ k ] == CSPBlock::kAllDifferent )
    for( auto v : s )
     f_alldiff[ v ].push_back( k );
   }
  }

 /// enumerates all the solutions
 void run( void ) {
  for( auto k : f_empty )
   if( ! f_b.is_satisfied( k , f_a ) )
    return;
  rec( 0 );
  }

 long count = 0;        ///< the number of solutions
 double best = Inf< double >();  ///< the best objective value
 bool wrong = false;    ///< true if a solution found is not feasible

 private:
 void rec( std::size_t i ) {
  if( i == f_a.size() ) {
   ++count;
   if( ! f_b.is_feasible( f_a ) )
    wrong = true;
   best = std::min( best , f_b.get_objective_value( f_a ) );
   return;
   }
  for( auto d : f_b.get_domain( i ) ) {
   f_a[ i ] = d;
   if( ok( i ) )
    rec( i + 1 );
   }
  }

 bool ok( std::size_t i ) const {
  for( auto k : f_alldiff[ i ] )
   for( auto v : f_b.get_scopes()[ k ] )
    if( ( std::size_t( v ) < i ) && ( f_a[ v ] == f_a[ i ] ) )
     return( false );
  for( auto k : f_check[ i ] )
   if( ! f_b.is_satisfied( k , f_a ) )
    return( false );
  return( true );
  }

 const CSPBlock & f_b;
 v_int f_a;
 std::vector< std::vector< Block::Index > > f_check;
 std::vector< std::vector< Block::Index > > f_alldiff;
 std::vector< Block::Index > f_empty;
 };

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 if( argc < 4 ) {
  std::cerr << "usage: " << argv[ 0 ] << " <family> <size> <solutions | ->"
	    << std::endl;
  return( 1 );
  }
 const std::string family = argv[ 1 ];
 const int size = std::stoi( argv[ 2 ] );
 long expected = ( std::string( argv[ 3 ] ) == "-" ) ? -1
						      : std::stol( argv[ 3 ] );

 std::string text;
 double optimum = 0;
 if( family == "queens" )
  text = queens( size );
 else
  if( family == "mapcolouring" )
   text = map_colouring( size );
  else
   if( family == "sudoku" )
    text = sudoku( size );
   else
    if( family == "scheduling" ) {
     text = scheduling( size );
     const auto c = scheduling_count( size );
     if( expected < 0 )
      expected = c.first;
     optimum = c.second;
     }
    else {
     std::cerr << "unknown family " << family << std::endl;
     return( 1 );
     }

 CSPBlock b;
 std::istringstream in( text );
 b.load( in );

 // the round trips through netCDF and through FlatZinc
 const std::string file = "CSPBlock_battery_" + family + "_" +
                          std::to_string( size ) + ".nc4";
 bool trips = true;
 {
  netCDF::NcFile f( file , netCDF::NcFile::replace );
  auto g = f.addGroup( "Block" );
  b.serialize( g );
  }
 {
  netCDF::NcFile f( file , netCDF::NcFile::read );
  auto d = Block::new_Block( f.getGroup( "Block" ) );
  auto c = dynamic_cast< CSPBlock * >( d );
  trips = c && same( *c , b );
  delete d;
  }
 std::remove( file.c_str() );
 {
  std::ostringstream full;
  b.print( full , 'C' );
  CSPBlock p;
  std::istringstream pin( full.str() );
  p.load( pin );
  trips = trips && same( p , b );
  }

 Enumerator e( b );
 e.run();

 const bool ok = trips && ( ! e.wrong ) && ( e.count == expected ) &&
                 ( ( family != "scheduling" ) || ( e.best == optimum ) );
 std::cout << family << " " << size << ": " << b.get_number_variables()
	   << " variables, " << b.get_number_constraints()
	   << " constraints, " << e.count << " solutions (expected "
	   << expected << ")";
 if( family == "scheduling" )
  std::cout << ", optimal makespan " << e.best << " (expected " << optimum
	    << ")";
 if( ! trips )
  std::cout << ", round trip FAILED";
 if( e.wrong )
  std::cout << ", a solution not feasible";
 std::cout << ( ok ? ": OK" : ": FAILED" ) << std::endl;
 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- End File battery.cpp -------------------------*/
/*--------------------------------------------------------------------------*/

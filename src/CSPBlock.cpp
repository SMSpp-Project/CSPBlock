/*--------------------------------------------------------------------------*/
/*--------------------------- File CSPBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the CSPBlock class, of its FlatZinc reader and writer,
 * and of the CSPBlockChange class.
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

#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>

#include "CSPBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register CSPBlock in the Block factory

SMSpp_insert_in_factory_cpp_1( CSPBlock );

// register CSPBlockChange in the Change factory

SMSpp_insert_in_factory_cpp_1( CSPBlockChange );

/*--------------------------------------------------------------------------*/
/*---------------------------- LOCAL FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

/// sorts the values and keeps each once

void sort_unique( CSPBlock::v_int & v )
{
 std::sort( v.begin() , v.end() );
 v.erase( std::unique( v.begin() , v.end() ) , v.end() );
 }

/*--------------------------------------------------------------------------*/
/// the ragged rows as a vector of entries and a vector of starts

std::pair< CSPBlock::v_int , std::vector< std::size_t > > flatten(
					      CSPBlock::v_v_int && rows )
{
 std::pair< CSPBlock::v_int , std::vector< std::size_t > > f;
 std::size_t tot = 0;
 for( const auto & r : rows )
  tot += r.size();
 f.first.reserve( tot );
 f.second.reserve( rows.size() + 1 );
 f.second.push_back( 0 );
 for( auto & r : rows ) {
  f.first.insert( f.first.end() , r.begin() , r.end() );
  f.second.push_back( f.first.size() );
  CSPBlock::v_int().swap( r );  // its memory goes at once
  }
 return( f );
 }

/*--------------------------------------------------------------------------*/
/// writes rows of int as entries over a dimension and starts over \p rows

void put_ragged( netCDF::NcGroup & group , const CSPBlock::RaggedArray & r ,
		 const std::string & name , const std::string & dim ,
		 const netCDF::NcDim & rows )
{
 if( r.values().empty() )  // all the rows are empty
  return;
 auto d = group.addDim( dim , r.values().size() );
 group.addVar( name , netCDF::NcInt() , d ).putVar( r.values().data() );
 std::vector< unsigned int > start( r.starts().begin() ,
				    r.starts().end() - 1 );
 group.addVar( name + "Start" , netCDF::NcUint() , rows ).putVar(
							     start.data() );
 }

/*--------------------------------------------------------------------------*/
/// reads \p n rows of int written by put_ragged()

void get_ragged( const netCDF::NcGroup & group , CSPBlock::RaggedArray & r ,
		 const std::string & name , std::size_t n )
{
 auto v = group.getVar( name );
 auto s = group.getVar( name + "Start" );
 if( v.isNull() != s.isNull() )
  throw( std::invalid_argument( "CSPBlock::deserialize: " + name + " and " +
				name + "Start must be both there or neither"
				) );
 if( v.isNull() ) {  // all the rows are empty
  r.assign( {} , std::vector< std::size_t >( n + 1 , 0 ) );
  return;
  }
 if( ( v.getDimCount() != 1 ) || ( s.getDimCount() != 1 ) ||
     ( s.getDim( 0 ).getSize() != n ) )
  throw( std::invalid_argument( "CSPBlock::deserialize: " + name +
				"Start must have one element per row" ) );
 CSPBlock::v_int val( v.getDim( 0 ).getSize() );
 if( ! val.empty() )
  v.getVar( val.data() );
 std::vector< unsigned int > st( n );
 if( n )
  s.getVar( st.data() );
 std::vector< std::size_t > start( st.begin() , st.end() );
 start.push_back( val.size() );
 r.assign( std::move( val ) , std::move( start ) );
 }

/*--------------------------------------------------------------------------*/
/// reads a 0/1 ubyte variable over a dimension of size \p n, all 1 if absent

void get_flags( const netCDF::NcGroup & group , CSPBlock::v_flag & f ,
		const std::string & name , std::size_t n )
{
 f.assign( n , 1 );
 auto v = group.getVar( name );
 if( v.isNull() || ( ! n ) )
  return;
 if( ( v.getDimCount() != 1 ) || ( v.getDim( 0 ).getSize() != n ) )
  throw( std::invalid_argument( "CSPBlock::deserialize: " + name +
				" has the wrong size" ) );
 v.getVar( f.data() );
 for( auto & x : f )
  x = ( x != 0 );
 }

/*--------------------------------------------------------------------------*/
/// reads a double variable over a dimension of size \p n, \p dflt if absent

void get_doubles( const netCDF::NcGroup & group , CSPBlock::v_dbl & f ,
		  const std::string & name , std::size_t n , double dflt )
{
 f.assign( n , dflt );
 auto v = group.getVar( name );
 if( v.isNull() || ( ! n ) )
  return;
 if( ( v.getDimCount() != 1 ) || ( v.getDim( 0 ).getSize() != n ) )
  throw( std::invalid_argument( "CSPBlock::deserialize: " + name +
				" has the wrong size" ) );
 v.getVar( f.data() );
 }

/*--------------------------------------------------------------------------*/
/// what a change of the weight from \p old to \p w does

Modification::ModConcern weight_changes( double old , double w )
{
 const double inf = Inf< double >();
 if( ( old == inf ) == ( w == inf ) )
  return( Modification::ePhysObj );
 return( Modification::ePhysObj | Modification::eModCnsSide |
	 ( w == inf ? Modification::eRegnShrink : Modification::eRegnGrow ) );
 }

/*--------------------------------------------------------------------------*/
/// sorts the Subset \p nms, and the elements of \p v along with it

template< class T >
void sort_along( Block::Subset & nms , std::vector< T > & v )
{
 std::vector< std::size_t > p( nms.size() );
 std::iota( p.begin() , p.end() , 0 );
 std::sort( p.begin() , p.end() , [ & ]( auto a , auto b ) {
  return( nms[ a ] < nms[ b ] ); } );
 Block::Subset sn( nms.size() );
 std::vector< T > sv;
 sv.reserve( nms.size() );
 for( std::size_t k = 0 ; k < p.size() ; ++k ) {
  sn[ k ] = nms[ p[ k ] ];
  sv.push_back( std::move( v[ p[ k ] ] ) );
  }
 nms = std::move( sn );
 v = std::move( sv );
 }

/*--------------------------------------------------------------------------*/
/// throws if the ordered Subset \p nms is repeated or is not below \p n

void check_subset( const Block::Subset & nms , std::size_t n ,
		   const std::string & who , const char * what )
{
 for( std::size_t k = 0 ; k < nms.size() ; ++k ) {
  if( nms[ k ] >= n )
   throw( std::invalid_argument( who + ": " + what + " " +
				 std::to_string( nms[ k ] ) +
				 " does not exist" ) );
  if( k && ( nms[ k ] == nms[ k - 1 ] ) )
   throw( std::invalid_argument( who + ": " + what + " " +
				 std::to_string( nms[ k ] ) + " repeated" ) );
  }
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- FlatZinc reader ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a token of FlatZinc: an int, a float, a name, a string or a symbol

struct Token {
 enum { kInt , kFloat , kName , kString , kSym } t;
 std::string s;     ///< the name, the string or the symbol
 long long i = 0;   ///< the int
 double d = 0;      ///< the float (also set for an int)
 };

/*--------------------------------------------------------------------------*/
/// an expression of FlatZinc

struct Expr {
 enum { kInt , kFloat , kName , kString , kArray , kSet , kRange , kCall } t;
 long long i = 0 , hi = 0;   ///< the int, or the range i .. hi
 double d = 0;               ///< the float (also set for an int)
 std::string s;              ///< the name, of a call too, or the string
 std::vector< Expr > e;      ///< the elements, or the arguments of a call
 };

/*--------------------------------------------------------------------------*/
/// the reader of the FlatZinc statements into a CSPBlock

class FznReader {
 public:
 using v_int = CSPBlock::v_int;
 using v_v_int = CSPBlock::v_v_int;

 /// reads the whole input
 explicit FznReader( std::istream & in ) {
  std::string text( ( std::istreambuf_iterator< char >( in ) ) ,
		    std::istreambuf_iterator< char >() );
  tokenize( text );
  while( f_p < f_tok.size() )
   statement();
  }

 v_v_int universes;   ///< the universe of each variable
 v_v_int domains;     ///< the current domain of each variable
 std::vector< CSPBlock::v_dbl > costs;  ///< the costs, empty if none
 v_int groups;        ///< the group of each variable, -1 if none
 v_int kinds;         ///< the kind of each constraint
 v_v_int scopes;      ///< the scope of each constraint
 v_v_int params;      ///< the parameters of each constraint
 CSPBlock::v_dbl weights;  ///< the weight of each constraint
 CSPBlock::v_flag active;  ///< the activity of each constraint
 int obj_var = -1;    ///< the variable of the solve item, if any
 int obj_sign = 0;    ///< 1 to minimize it, -1 to maximize it

 private:
 static constexpr long long kMaxValues = 1 << 24;  ///< values of a range

 [[noreturn]] void fail( const std::string & what ) const {
  throw( std::invalid_argument( "CSPBlock::load: " + what ) );
  }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// splits the text into tokens, the comments left out
 void tokenize( const std::string & s ) {
  std::size_t i = 0;
  while( i < s.size() ) {
   const char c = s[ i ];
   if( std::isspace( static_cast< unsigned char >( c ) ) ) { ++i; continue; }
   if( c == '%' ) {
    while( ( i < s.size() ) && ( s[ i ] != '\n' ) )
     ++i;
    continue;
    }
   Token t;
   if( c == '"' ) {
    std::size_t j = i + 1;
    while( ( j < s.size() ) && ( s[ j ] != '"' ) )
     j += ( s[ j ] == '\\' ) ? 2 : 1;
    if( j >= s.size() )
     fail( "unterminated string" );
    t.t = Token::kString;
    t.s = s.substr( i + 1 , j - i - 1 );
    i = j + 1;
    }
   else
    if( std::isdigit( static_cast< unsigned char >( c ) ) ||
	( ( c == '-' ) && ( i + 1 < s.size() ) &&
	  std::isdigit( static_cast< unsigned char >( s[ i + 1 ] ) ) ) ) {
     std::size_t j = i + 1;
     while( ( j < s.size() ) &&
	    std::isdigit( static_cast< unsigned char >( s[ j ] ) ) )
      ++j;
     bool flt = false;
     if( ( j + 1 < s.size() ) && ( s[ j ] == '.' ) &&
	 std::isdigit( static_cast< unsigned char >( s[ j + 1 ] ) ) ) {
      flt = true;  // a fraction, not the .. of a range
      for( ++j ; ( j < s.size() ) &&
	     std::isdigit( static_cast< unsigned char >( s[ j ] ) ) ; )
       ++j;
      }
     if( ( j < s.size() ) && ( ( s[ j ] == 'e' ) || ( s[ j ] == 'E' ) ) ) {
      std::size_t k = j + 1;
      if( ( k < s.size() ) && ( ( s[ k ] == '+' ) || ( s[ k ] == '-' ) ) )
       ++k;
      if( ( k < s.size() ) &&
	  std::isdigit( static_cast< unsigned char >( s[ k ] ) ) ) {
       flt = true;
       for( j = k ; ( j < s.size() ) &&
	      std::isdigit( static_cast< unsigned char >( s[ j ] ) ) ; )
	++j;
       }
      }
     const std::string num = s.substr( i , j - i );
     try {
      if( flt ) {
       t.t = Token::kFloat;
       t.d = std::stod( num );
       }
      else {
       t.t = Token::kInt;
       t.i = std::stoll( num );
       t.d = double( t.i );
       }
      }
     catch( std::exception & ) {
      fail( "number " + num + " out of range" );
      }
     i = j;
     }
    else
     if( std::isalpha( static_cast< unsigned char >( c ) ) || ( c == '_' ) ) {
      std::size_t j = i + 1;
      while( ( j < s.size() ) &&
	     ( std::isalnum( static_cast< unsigned char >( s[ j ] ) ) ||
	       ( s[ j ] == '_' ) ) )
       ++j;
      t.t = Token::kName;
      t.s = s.substr( i , j - i );
      i = j;
      }
     else {
      t.t = Token::kSym;
      if( ( ( c == ':' ) || ( c == '.' ) ) && ( i + 1 < s.size() ) &&
	  ( s[ i + 1 ] == c ) ) {
       t.s = std::string( 2 , c );
       i += 2;
       }
      else {
       if( std::string( "[](){},:=;" ).find( c ) == std::string::npos )
	fail( std::string( "unexpected character '" ) + c + "'" );
       t.s = std::string( 1 , c );
       ++i;
       }
      }
   f_tok.push_back( std::move( t ) );
   }
  }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// true if the current token is the symbol \p sym
 bool is_sym( const char * sym ) const {
  return( ( f_p < f_tok.size() ) && ( f_tok[ f_p ].t == Token::kSym ) &&
	  ( f_tok[ f_p ].s == sym ) );
  }

 /// true if the current token is the name \p nm
 bool is_name( const char * nm ) const {
  return( ( f_p < f_tok.size() ) && ( f_tok[ f_p ].t == Token::kName ) &&
	  ( f_tok[ f_p ].s == nm ) );
  }

 /// consumes the symbol \p sym, which must be there
 void expect( const char * sym ) {
  if( ! is_sym( sym ) )
   fail( std::string( "'" ) + sym + "' expected" +
	 ( f_p < f_tok.size() ? " before '" + f_tok[ f_p ].s + "'" : "" ) );
  ++f_p;
  }

 /// consumes a name, which must be there, and returns it
 std::string name( void ) {
  if( ( f_p >= f_tok.size() ) || ( f_tok[ f_p ].t != Token::kName ) )
   fail( "a name expected" );
  return( f_tok[ f_p++ ].s );
  }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// reads an expression
 Expr expr( void ) {
  if( f_p >= f_tok.size() )
   fail( "unexpected end of the input" );
  Expr x;
  const Token & t = f_tok[ f_p++ ];
  switch( t.t ) {
   case( Token::kInt ):
    x.t = Expr::kInt;
    x.i = t.i;
    x.d = t.d;
    if( is_sym( ".." ) ) {
     ++f_p;
     if( ( f_p >= f_tok.size() ) || ( f_tok[ f_p ].t != Token::kInt ) )
      fail( "an int expected after '..'" );
     x.t = Expr::kRange;
     x.hi = f_tok[ f_p++ ].i;
     }
    return( x );
   case( Token::kFloat ):
    x.t = Expr::kFloat;
    x.d = t.d;
    return( x );
   case( Token::kString ):
    x.t = Expr::kString;
    x.s = t.s;
    return( x );
   case( Token::kName ):
    if( ( t.s == "true" ) || ( t.s == "false" ) ) {
     x.t = Expr::kInt;
     x.i = ( t.s == "true" );
     x.d = double( x.i );
     return( x );
     }
    x.s = t.s;
    if( is_sym( "(" ) ) {
     x.t = Expr::kCall;
     ++f_p;
     x.e = list( ")" );
     }
    else
     x.t = Expr::kName;
    return( x );
   default:
    if( t.s == "[" ) {
     x.t = Expr::kArray;
     x.e = list( "]" );
     return( x );
     }
    if( t.s == "{" ) {
     x.t = Expr::kSet;
     x.e = list( "}" );
     return( x );
     }
    fail( "unexpected '" + t.s + "'" );
   }
  }

 /// reads the expressions separated by commas up to the symbol \p close
 std::vector< Expr > list( const char * close ) {
  std::vector< Expr > v;
  if( is_sym( close ) ) {
   ++f_p;
   return( v );
   }
  for( ; ; ) {
   v.push_back( expr() );
   if( is_sym( close ) ) {
    ++f_p;
    return( v );
    }
   expect( "," );
   }
  }

 /// reads the annotations, if any
 std::vector< Expr > annotations( void ) {
  std::vector< Expr > v;
  while( is_sym( "::" ) ) {
   ++f_p;
   v.push_back( expr() );
   }
  return( v );
  }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// the int of an expression: an int or a parameter
 long long int1( const Expr & x ) const {
  if( x.t == Expr::kInt )
   return( x.i );
  if( x.t == Expr::kName ) {
   auto it = m_par.find( x.s );
   if( it != m_par.end() )
    return( it->second );
   }
  fail( "an int expected" );
  }

 /// an int that fits an int
 int to_int( long long v ) const {
  if( ( v < INT_MIN ) || ( v > INT_MAX ) )
   fail( "value " + std::to_string( v ) + " out of the range of int" );
  return( int( v ) );
  }

 /// the ints of an expression: an array of ints or of parameters, or an
 /// array parameter
 v_int ints( const Expr & x ) const {
  if( x.t == Expr::kName ) {
   auto it = m_par_array.find( x.s );
   if( it == m_par_array.end() )
    fail( "array of int " + x.s + " not declared" );
   return( it->second );
   }
  if( x.t != Expr::kArray )
   fail( "an array of int expected" );
  v_int v;
  v.reserve( x.e.size() );
  for( const auto & y : x.e )
   v.push_back( to_int( int1( y ) ) );
  return( v );
  }

 /// the values of a set: a range or a set of ints
 v_int set( const Expr & x ) const {
  v_int v;
  if( x.t == Expr::kRange ) {
   if( x.hi - x.i >= kMaxValues )
    fail( "range " + std::to_string( x.i ) + ".." + std::to_string( x.hi ) +
	  " too large" );
   for( long long d = x.i ; d <= x.hi ; ++d )
    v.push_back( to_int( d ) );
   return( v );
   }
  if( x.t == Expr::kSet )
   for( const auto & y : x.e )
    if( y.t == Expr::kRange ) {
     auto r = set( y );
     v.insert( v.end() , r.begin() , r.end() );
     }
    else
     v.push_back( to_int( int1( y ) ) );
  else
   if( x.t == Expr::kName ) {
    auto it = m_par_set.find( x.s );
    if( it == m_par_set.end() )
     fail( "set " + x.s + " not declared" );
    return( it->second );
    }
   else
    fail( "a set of int expected" );
  sort_unique( v );
  return( v );
  }

 /// a term of an array of variables: a variable (>= 0), or a value
 struct Term { int var; long long val; };

 /// the term of an expression where a variable is expected
 Term term( const Expr & x ) const {
  if( x.t == Expr::kName ) {
   auto it = m_var.find( x.s );
   if( it != m_var.end() )
    return( Term{ it->second , 0 } );
   auto p = m_par.find( x.s );
   if( p != m_par.end() )
    return( Term{ -1 , p->second } );
   fail( "variable " + x.s + " not declared" );
   }
  return( Term{ -1 , int1( x ) } );
  }

 /// the terms of an array of variables, or of an array variable
 std::vector< Term > terms( const Expr & x ) const {
  if( x.t == Expr::kName ) {
   auto it = m_var_array.find( x.s );
   if( it != m_var_array.end() )
    return( it->second );
   auto p = m_par_array.find( x.s );
   if( p != m_par_array.end() ) {
    std::vector< Term > v;
    for( auto d : p->second )
     v.push_back( { -1 , d } );
    return( v );
    }
   fail( "array " + x.s + " not declared" );
   }
  if( x.t != Expr::kArray )
   fail( "an array of variables expected" );
  std::vector< Term > v;
  v.reserve( x.e.size() );
  for( const auto & y : x.e )
   v.push_back( term( y ) );
  return( v );
  }

 /// the variable of a term, a value being the variable of that value only
 int var( const Term & t ) {
  if( t.var >= 0 )
   return( t.var );
  const int d = to_int( t.val );
  auto it = m_fixed.find( d );
  if( it != m_fixed.end() )
   return( it->second );
  const int v = new_var( { d } );
  m_fixed[ d ] = v;
  return( v );
  }

 /// the variables of the terms
 v_int vars( const std::vector< Term > & t ) {
  v_int v;
  v.reserve( t.size() );
  for( const auto & x : t )
   v.push_back( var( x ) );
  return( v );
  }

 /// adds a variable with the given universe
 int new_var( v_int && u ) {
  universes.push_back( u );
  domains.push_back( std::move( u ) );
  costs.emplace_back();
  groups.push_back( -1 );
  return( int( universes.size() - 1 ) );
  }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// reads a statement, up to its ';'
 void statement( void ) {
  if( is_name( "predicate" ) )
   while( ( f_p < f_tok.size() ) && ( ! is_sym( ";" ) ) )
    ++f_p;
  else
   if( is_name( "var" ) )
    var_decl();
   else
    if( is_name( "array" ) )
     array_decl();
    else
     if( is_name( "constraint" ) )
      constraint();
     else
      if( is_name( "solve" ) )
       solve();
      else
       par_decl();
  expect( ";" );
  }

 /// reads the type of a variable after "var": its values
 v_int var_type( void ) {
  if( is_name( "bool" ) ) {
   ++f_p;
   return( v_int{ 0 , 1 } );
   }
  if( is_name( "int" ) || is_name( "float" ) )
   fail( "variable of type " + f_tok[ f_p ].s + " without a finite set of "
	 "values" );
  if( is_name( "set" ) )
   fail( "set variables are not supported" );
  return( set( expr() ) );
  }

 /// reads the declaration of a variable
 void var_decl( void ) {
  ++f_p;  // "var"
  auto dom = var_type();
  expect( ":" );
  const auto nm = name();
  const auto ann = annotations();
  if( m_var.count( nm ) || m_par.count( nm ) )
   fail( nm + " declared twice" );
  if( is_sym( "=" ) ) {
   ++f_p;
   const auto t = term( expr() );
   if( t.var >= 0 ) {  // the same variable
    m_var[ nm ] = t.var;
    auto & d = domains[ t.var ];
    v_int both;
    std::set_intersection( d.begin() , d.end() , dom.begin() , dom.end() ,
			   std::back_inserter( both ) );
    d = std::move( both );
    return;
    }
   const int v = new_var( std::move( dom ) );
   const int d = to_int( t.val );
   auto & cd = domains[ v ];
   const bool in = std::binary_search( cd.begin() , cd.end() , d );
   cd.clear();
   if( in )
    cd.push_back( d );
   m_var[ nm ] = v;
   var_annotations( v , ann );
   return;
   }
  const int v = new_var( std::move( dom ) );
  m_var[ nm ] = v;
  var_annotations( v , ann );
  }

 /// applies the annotations universe(), cost() and group() of a variable
 void var_annotations( int v , const std::vector< Expr > & ann ) {
  for( const auto & a : ann ) {
   if( ( a.t != Expr::kCall ) || ( a.e.size() != 1 ) )
    continue;
   if( a.s == "universe" ) {
    auto u = set( a.e[ 0 ] );
    if( ! std::includes( u.begin() , u.end() , domains[ v ].begin() ,
			 domains[ v ].end() ) )
     fail( "the universe of a variable must hold its domain" );
    universes[ v ] = std::move( u );
    }
   }
  for( const auto & a : ann ) {
   if( ( a.t != Expr::kCall ) || ( a.e.size() != 1 ) )
    continue;
   if( a.s == "cost" ) {
    if( a.e[ 0 ].t != Expr::kArray )
     fail( "cost() wants an array" );
    CSPBlock::v_dbl c;
    for( const auto & y : a.e[ 0 ].e ) {
     if( ( y.t != Expr::kInt ) && ( y.t != Expr::kFloat ) )
      fail( "cost() wants numbers" );
     c.push_back( y.d );
     }
    if( c.size() != universes[ v ].size() )
     fail( "cost() wants one cost per value of the universe" );
    costs[ v ] = std::move( c );
    }
   else
    if( a.s == "group" ) {
     const auto g = int1( a.e[ 0 ] );
     if( g < 0 )
      fail( "a group cannot be negative" );
     groups[ v ] = to_int( g );
     }
   }
  }

 /// reads the declaration of an array
 void array_decl( void ) {
  ++f_p;  // "array"
  expect( "[" );
  const auto ix = expr();
  if( ix.t != Expr::kRange )
   fail( "the index set of an array must be a range" );
  expect( "]" );
  if( ! is_name( "of" ) )
   fail( "'of' expected" );
  ++f_p;
  if( is_name( "var" ) ) {
   ++f_p;
   if( is_name( "int" ) || is_name( "bool" ) )
    ++f_p;  // the values are those of the variables in the array
   else
    var_type();
   expect( ":" );
   const auto nm = name();
   annotations();
   expect( "=" );
   m_var_array[ nm ] = terms( expr() );
   return;
   }
  if( is_name( "set" ) ) {  // arrays of sets are not used
   while( ! is_sym( ";" ) && ( f_p < f_tok.size() ) )
    ++f_p;
   return;
   }
  name();  // int or bool
  expect( ":" );
  const auto nm = name();
  annotations();
  expect( "=" );
  m_par_array[ nm ] = ints( expr() );
  }

 /// reads the declaration of a parameter
 void par_decl( void ) {
  if( is_name( "set" ) ) {
   ++f_p;
   if( ! is_name( "of" ) )
    fail( "'of' expected" );
   ++f_p;
   name();
   expect( ":" );
   const auto nm = name();
   annotations();
   expect( "=" );
   m_par_set[ nm ] = set( expr() );
   return;
   }
  const auto ty = name();
  if( ( ty != "int" ) && ( ty != "bool" ) )
   fail( "unsupported statement beginning with " + ty );
  expect( ":" );
  const auto nm = name();
  annotations();
  expect( "=" );
  m_par[ nm ] = int1( expr() );
  }

 /// adds a constraint
 void add( int kind , v_int && scope , v_int && param ) {
  kinds.push_back( kind );
  scopes.push_back( std::move( scope ) );
  params.push_back( std::move( param ) );
  weights.push_back( Inf< double >() );
  active.push_back( 1 );
  }

 /// adds a linear constraint, the values among the terms in the right-hand
 /// side
 void linear( const v_int & a , const std::vector< Term > & x , int op ,
	      long long b ) {
  if( a.size() != x.size() )
   fail( "a linear constraint wants as many coefficients as variables" );
  v_int scope , param;
  for( std::size_t i = 0 ; i < a.size() ; ++i )
   if( x[ i ].var >= 0 ) {
    scope.push_back( x[ i ].var );
    param.push_back( a[ i ] );
    }
   else
    b -= ( long long ) a[ i ] * x[ i ].val;
  param.push_back( op );
  param.push_back( to_int( b ) );
  add( CSPBlock::kLinear , std::move( scope ) , std::move( param ) );
  }

 /// reads a constraint
 void constraint( void ) {
  ++f_p;  // "constraint"
  const auto c = expr();
  if( c.t != Expr::kCall )
   fail( "a constraint expected" );
  const auto ann = annotations();
  const auto & a = c.e;
  const auto & n = c.s;
  const auto argc = [ & ]( std::size_t k ) {
   if( a.size() != k )
    fail( n + " wants " + std::to_string( k ) + " arguments" );
   };

  int op = -1;
  if( ( n == "int_lin_le" ) || ( n == "csp_lin_le" ) ) op = CSPBlock::kLE;
  if( ( n == "int_lin_eq" ) || ( n == "csp_lin_eq" ) ) op = CSPBlock::kEQ;
  if( ( n == "int_lin_ne" ) || ( n == "csp_lin_ne" ) ) op = CSPBlock::kNE;
  if( op >= 0 ) {
   argc( 3 );
   linear( ints( a[ 0 ] ) , terms( a[ 1 ] ) , op , int1( a[ 2 ] ) );
   }
  else {
   // the binary ones: a op b, as a - b op 0 or b - a op 0
   static const std::map< std::string , std::pair< int , int > > bin = {
    { "int_le" , { CSPBlock::kLE , 0 } } ,
    { "int_lt" , { CSPBlock::kLE , 1 } } ,
    { "int_ge" , { CSPBlock::kLE , 2 } } ,
    { "int_gt" , { CSPBlock::kLE , 3 } } ,
    { "int_eq" , { CSPBlock::kEQ , 0 } } ,
    { "int_ne" , { CSPBlock::kNE , 0 } } ,
    { "bool2int" , { CSPBlock::kEQ , 0 } } };
   auto it = bin.find( n );
   if( it != bin.end() ) {
    argc( 2 );
    const int s = it->second.second;
    // a - b <= 0, a - b <= -1, b - a <= 0, b - a <= -1
    const v_int coef = ( s < 2 ) ? v_int{ 1 , -1 } : v_int{ -1 , 1 };
    linear( coef , { term( a[ 0 ] ) , term( a[ 1 ] ) } , it->second.first ,
	    ( s % 2 ) ? -1 : 0 );
    }
   else
    if( ( n == "fzn_all_different_int" ) || ( n == "all_different_int" ) ) {
     argc( 1 );
     add( CSPBlock::kAllDifferent , vars( terms( a[ 0 ] ) ) , {} );
     }
    else
     if( ( n == "fzn_table_int" ) || ( n == "table_int" ) ||
	 ( n == "csp_table_int_neg" ) ) {
      argc( 2 );
      auto x = vars( terms( a[ 0 ] ) );
      auto t = ints( a[ 1 ] );
      if( x.empty() || ( t.size() % x.size() ) )
       fail( n + " wants tuples of as many values as the variables" );
      add( n == "csp_table_int_neg" ? CSPBlock::kConflicts
				    : CSPBlock::kExtension ,
	   std::move( x ) , std::move( t ) );
      }
     else
      if( n == "array_var_int_element" ) {
       argc( 3 );
       v_int s = { var( term( a[ 0 ] ) ) , var( term( a[ 2 ] ) ) };
       auto x = vars( terms( a[ 1 ] ) );
       s.insert( s.end() , x.begin() , x.end() );
       add( CSPBlock::kElement , std::move( s ) , { 1 } );
       }
      else
       if( n == "csp_element_offset" ) {
	argc( 4 );
	v_int s = { var( term( a[ 0 ] ) ) , var( term( a[ 3 ] ) ) };
	auto x = vars( terms( a[ 2 ] ) );
	s.insert( s.end() , x.begin() , x.end() );
	add( CSPBlock::kElement , std::move( s ) ,
	     { to_int( int1( a[ 1 ] ) ) } );
	}
       else
	if( n == "array_int_element" ) {
	 argc( 3 );
	 auto as = ints( a[ 1 ] );
	 v_int t;
	 for( std::size_t j = 0 ; j < as.size() ; ++j ) {
	  t.push_back( int( j + 1 ) );
	  t.push_back( as[ j ] );
	  }
	 add( CSPBlock::kExtension ,
	      { var( term( a[ 0 ] ) ) , var( term( a[ 2 ] ) ) } ,
	      std::move( t ) );
	 }
	else
	 if( n == "set_in" ) {
	  argc( 2 );
	  add( CSPBlock::kExtension , { var( term( a[ 0 ] ) ) } ,
	       set( a[ 1 ] ) );
	  }
	 else
	  if( ( n == "fzn_cumulative" ) || ( n == "cumulative" ) ) {
	   argc( 4 );
	   auto s = vars( terms( a[ 0 ] ) );
	   auto d = ints( a[ 1 ] );
	   auto r = ints( a[ 2 ] );
	   if( ( d.size() != s.size() ) || ( r.size() != s.size() ) )
	    fail( n + " wants one duration and one usage per task" );
	   d.insert( d.end() , r.begin() , r.end() );
	   d.push_back( to_int( int1( a[ 3 ] ) ) );
	   add( CSPBlock::kCumulative , std::move( s ) , std::move( d ) );
	   }
	  else
	   fail( "constraint " + n + " is not supported" );
   }

  for( const auto & x : ann )
   if( ( x.t == Expr::kName ) && ( x.s == "inactive" ) )
    active.back() = 0;
   else
    if( ( x.t == Expr::kCall ) && ( x.s == "weight" ) &&
	( x.e.size() == 1 ) ) {
     if( ( x.e[ 0 ].t != Expr::kInt ) && ( x.e[ 0 ].t != Expr::kFloat ) )
      fail( "weight() wants a number" );
     weights.back() = x.e[ 0 ].d;
     }
  }

 /// reads the solve item
 void solve( void ) {
  ++f_p;  // "solve"
  annotations();
  if( is_name( "satisfy" ) ) {
   ++f_p;
   return;
   }
  if( ! ( is_name( "minimize" ) || is_name( "maximize" ) ) )
   fail( "satisfy, minimize or maximize expected" );
  obj_sign = is_name( "minimize" ) ? 1 : -1;
  ++f_p;
  const auto t = term( expr() );
  if( t.var < 0 )  // a constant: no cost
   obj_sign = 0;
  else
   obj_var = t.var;
  }

 std::vector< Token > f_tok;  ///< the tokens
 std::size_t f_p = 0;         ///< the current one

 std::map< std::string , int > m_var;               ///< the variables
 std::map< std::string , std::vector< Term > > m_var_array;  ///< arrays
 std::map< std::string , long long > m_par;         ///< the int parameters
 std::map< std::string , v_int > m_par_array;       ///< the int arrays
 std::map< std::string , v_int > m_par_set;         ///< the set parameters
 std::map< int , int > m_fixed;  ///< the variable of each value
 };

/*--------------------------------------------------------------------------*/
/// writes the values as a FlatZinc set, a range if they are contiguous

void print_set( std::ostream & out , std::span< const int > v )
{
 if( ( v.size() > 1 ) &&
     ( ( long long ) v.back() - v.front() + 1 == ( long long ) v.size() ) ) {
  out << v.front() << ".." << v.back();
  return;
  }
 out << "{";
 for( std::size_t i = 0 ; i < v.size() ; ++i )
  out << ( i ? "," : "" ) << v[ i ];
 out << "}";
 }

/*--------------------------------------------------------------------------*/
/// writes the ints as a FlatZinc array, each with the given prefix

void print_array( std::ostream & out , std::span< const int > v ,
		  const char * prefix = "" )
{
 out << "[";
 for( std::size_t i = 0 ; i < v.size() ; ++i )
  out << ( i ? "," : "" ) << prefix << v[ i ];
 out << "]";
 }

/*--------------------------------------------------------------------------*/

 }  // end( anonymous namespace )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS OF RaggedArray --------------------------*/
/*--------------------------------------------------------------------------*/

void CSPBlock::RaggedArray::assign( v_int && val ,
				    std::vector< std::size_t > && start )
{
 if( start.empty() || start.front() || ( start.back() != val.size() ) ||
     ( ! std::is_sorted( start.begin() , start.end() ) ) )
  throw( std::invalid_argument( "CSPBlock::RaggedArray::assign: wrong "
				"starts of the rows" ) );
 v_val = std::move( val );
 v_start = std::move( start );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/

void CSPBlock::check_constraint( int kind , std::span< const int > scope ,
				 std::span< const int > param ) const
{
 const auto n = get_number_variables();
 const auto k = scope.size();
 const auto fail = [ & ]( const std::string & what ) {
  throw( std::invalid_argument( "CSPBlock::check_constraint: " + what ) );
  };

 for( auto v : scope )
  if( ( v < 0 ) || ( unsigned( v ) >= n ) )
   fail( "variable " + std::to_string( v ) + " out of range" );

 switch( kind ) {
  case( kExtension ):
  case( kConflicts ):
   if( ( ! k ) || ( param.size() % k ) )
    fail( "a table wants a nonempty scope and tuples of its size" );
   break;
  case( kLinear ):
   if( param.size() != k + 2 )
    fail( "a linear constraint wants k coefficients, op and b" );
   if( ( param[ k ] < kLE ) || ( param[ k ] > kNE ) )
    fail( "wrong operator " + std::to_string( param[ k ] ) );
   break;
  case( kAllDifferent ):
   if( ! param.empty() )
    fail( "an alldifferent constraint has no parameters" );
   break;
  case( kElement ):
   if( ( k < 3 ) || ( param.size() != 1 ) )
    fail( "an element constraint wants a scope i, z, x_0, ... and an "
	  "offset" );
   break;
  case( kCumulative ):
   if( param.size() != 2 * k + 1 )
    fail( "a cumulative constraint wants k durations, k usages and the "
	  "capacity" );
   if( std::any_of( param.begin() , param.end() ,
		    []( int x ) { return( x < 0 ); } ) )
    fail( "durations, usages and capacity cannot be negative" );
   break;
  default:
   fail( "wrong kind " + std::to_string( kind ) );
  }
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::check_weights( const v_dbl & weights , Index first ) const
{
 for( Index i = first ; i < weights.size() ; ++i )
  if( std::isnan( weights[ i ] ) || ( weights[ i ] < 0 ) )
   throw( std::invalid_argument( "CSPBlock::check_weights: weight " +
				 std::to_string( weights[ i ] ) +
				 " of constraint " + std::to_string( i ) ) );
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::guts_of_destructor( void )
{
 v_universe.clear();
 v_flag().swap( v_domain );
 v_dbl().swap( v_costs );
 v_int().swap( v_group );
 v_int().swap( v_kind );
 v_scope.clear();
 v_param.clear();
 v_dbl().swap( v_weights );
 v_flag().swap( v_active );
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::load( v_v_int && universes , v_int && kinds ,
		     v_v_int && scopes , v_v_int && params , v_dbl && weights )
{
 for( auto & u : universes ) {
  sort_unique( u );
  if( u.empty() )
   throw( std::invalid_argument( "CSPBlock::load: empty universe" ) );
  }
 if( ( scopes.size() != kinds.size() ) || ( params.size() != kinds.size() ) )
  throw( std::invalid_argument( "CSPBlock::load: kinds, scopes and "
				"parameters must have the same size" ) );
 if( weights.empty() )
  weights.assign( kinds.size() , Inf< double >() );
 if( weights.size() != kinds.size() )
  throw( std::invalid_argument( "CSPBlock::load: one weight per constraint "
				"is needed" ) );

 // the variables first, which the constraints are checked against
 CSPBlock nb;
 auto u = flatten( std::move( universes ) );
 nb.v_universe.assign( std::move( u.first ) , std::move( u.second ) );
 for( Index c = 0 ; c < kinds.size() ; ++c )
  nb.check_constraint( kinds[ c ] , scopes[ c ] , params[ c ] );
 check_weights( weights );

 guts_of_destructor();
 v_universe = std::move( nb.v_universe );
 v_domain.assign( get_number_values() , 1 );
 v_costs.assign( get_number_values() , 0 );
 v_kind = std::move( kinds );
 auto s = flatten( std::move( scopes ) );
 v_scope.assign( std::move( s.first ) , std::move( s.second ) );
 auto p = flatten( std::move( params ) );
 v_param.assign( std::move( p.first ) , std::move( p.second ) );
 v_weights = std::move( weights );
 v_active.assign( v_kind.size() , 1 );
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::load( std::istream & input , char frmt )
{
 if( frmt && ( frmt != 'F' ) )
  throw( std::invalid_argument( std::string( "CSPBlock::load: unknown "
					     "format " ) + frmt ) );

 FznReader r( input );

 const auto dom = std::move( r.domains );
 auto costs = std::move( r.costs );
 auto groups = std::move( r.groups );
 auto active = std::move( r.active );
 const int ov = r.obj_var;
 const int os = r.obj_sign;

 load( std::move( r.universes ) , std::move( r.kinds ) ,
       std::move( r.scopes ) , std::move( r.params ) ,
       std::move( r.weights ) );

 // the current domains, the costs, the activity and the groups
 for( Index i = 0 ; i < dom.size() ; ++i ) {
  const auto u = get_universe( i );
  const auto first = v_universe.starts()[ i ];
  for( Index j = 0 ; j < u.size() ; ++j )
   v_domain[ first + j ] = std::binary_search( dom[ i ].begin() ,
					       dom[ i ].end() , u[ j ] );
  if( ! costs[ i ].empty() )
   std::copy( costs[ i ].begin() , costs[ i ].end() ,
	      v_costs.begin() + first );
  if( int( i ) == ov )
   for( Index j = 0 ; j < u.size() ; ++j )
    v_costs[ first + j ] += os * double( u[ j ] );
  }
 if( ! std::all_of( v_costs.begin() , v_costs.end() ,
		    []( double c ) { return( std::isfinite( c ) ); } ) )
  throw( std::invalid_argument( "CSPBlock::load: a cost is not finite" ) );
 v_active = std::move( active );
 if( std::any_of( groups.begin() , groups.end() ,
		  []( int g ) { return( g >= 0 ); } ) ) {
  if( std::any_of( groups.begin() , groups.end() ,
		   []( int g ) { return( g < 0 ); } ) )
   throw( std::invalid_argument( "CSPBlock::load: either all the variables "
				 "have a group or none has" ) );
  set_variable_groups( std::move( groups ) );
  }
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::deserialize( const netCDF::NcGroup & group )
{
 auto nv = group.getDim( "NumberVariables" );
 if( nv.isNull() )
  throw( std::invalid_argument( "CSPBlock::deserialize: dimension "
				"NumberVariables is required" ) );
 const auto n = nv.getSize();

 CSPBlock nb;
 get_ragged( group , nb.v_universe , "Universe" , n );
 for( Index i = 0 ; i < n ; ++i ) {
  const auto u = nb.get_universe( i );
  if( u.empty() || ( std::adjacent_find( u.begin() , u.end() ,
					 std::greater_equal< int >() ) !=
		     u.end() ) )
   throw( std::invalid_argument( "CSPBlock::deserialize: the universe of "
				 "variable " + std::to_string( i ) +
				 " is empty or not sorted" ) );
  }
 const auto nval = nb.get_number_values();
 get_flags( group , nb.v_domain , "Domain" , nval );
 get_doubles( group , nb.v_costs , "Costs" , nval , 0 );
 if( ! std::all_of( nb.v_costs.begin() , nb.v_costs.end() ,
		    []( double c ) { return( std::isfinite( c ) ); } ) )
  throw( std::invalid_argument( "CSPBlock::deserialize: a cost is not "
				"finite" ) );

 auto g = group.getVar( "VariableGroups" );
 if( ( ! g.isNull() ) && n ) {
  if( ( g.getDimCount() != 1 ) || ( g.getDim( 0 ).getSize() != n ) )
   throw( std::invalid_argument( "CSPBlock::deserialize: VariableGroups "
				 "must have one element per variable" ) );
  v_int groups( n );
  g.getVar( groups.data() );
  nb.set_variable_groups( std::move( groups ) );
  }

 auto nc = group.getDim( "NumberConstraints" );
 const std::size_t m = nc.isNull() ? 0 : nc.getSize();
 if( m ) {
  nb.v_kind.resize( m );
  auto k = group.getVar( "ConstraintKind" );
  if( k.isNull() || ( k.getDimCount() != 1 ) ||
      ( k.getDim( 0 ).getSize() != m ) )
   throw( std::invalid_argument( "CSPBlock::deserialize: ConstraintKind "
				 "must have one element per constraint" ) );
  k.getVar( nb.v_kind.data() );
  }
 get_ragged( group , nb.v_scope , "Scope" , m );
 get_ragged( group , nb.v_param , "Parameters" , m );
 for( Index c = 0 ; c < m ; ++c )
  nb.check_constraint( nb.v_kind[ c ] , nb.v_scope[ c ] , nb.v_param[ c ] );
 get_doubles( group , nb.v_weights , "Weights" , m , Inf< double >() );
 check_weights( nb.v_weights );
 get_flags( group , nb.v_active , "Active" , m );

 guts_of_destructor();
 v_universe = std::move( nb.v_universe );
 v_domain = std::move( nb.v_domain );
 v_costs = std::move( nb.v_costs );
 v_group = std::move( nb.v_group );
 v_kind = std::move( nb.v_kind );
 v_scope = std::move( nb.v_scope );
 v_param = std::move( nb.v_param );
 v_weights = std::move( nb.v_weights );
 v_active = std::move( nb.v_active );

 Block::deserialize( group );
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::set_variable_groups( v_int && groups )
{
 if( ( ! groups.empty() ) && ( groups.size() != get_number_variables() ) )
  throw( std::invalid_argument( "CSPBlock::set_variable_groups: one group "
				"per variable is needed" ) );
 if( std::any_of( groups.begin() , groups.end() ,
		  []( int g ) { return( g < 0 ); } ) )
  throw( std::invalid_argument( "CSPBlock::set_variable_groups: a group "
				"cannot be negative" ) );
 v_group = std::move( groups );
 }

/*--------------------------------------------------------------------------*/
/*----------------- Methods for reading the data of the CSPBlock -----------*/
/*--------------------------------------------------------------------------*/

Block::Index CSPBlock::value_index( Index i , int d ) const
{
 if( i >= get_number_variables() )
  return( Inf< Index >() );
 const auto u = get_universe( i );
 auto it = std::lower_bound( u.begin() , u.end() , d );
 if( ( it == u.end() ) || ( *it != d ) )
  return( Inf< Index >() );
 return( v_universe.starts()[ i ] + Index( it - u.begin() ) );
 }

/*--------------------------------------------------------------------------*/

CSPBlock::v_int CSPBlock::get_domain( Index i ) const
{
 v_int d;
 const auto u = get_universe( i );
 const auto first = v_universe.starts()[ i ];
 for( Index j = 0 ; j < u.size() ; ++j )
  if( v_domain[ first + j ] )
   d.push_back( u[ j ] );
 return( d );
 }

/*--------------------------------------------------------------------------*/

bool CSPBlock::has_costs( void ) const
{
 return( std::any_of( v_costs.begin() , v_costs.end() ,
		      []( double c ) { return( c != 0 ); } ) );
 }

/*--------------------------------------------------------------------------*/

bool CSPBlock::all_hard( void ) const
{
 return( std::all_of( v_weights.begin() , v_weights.end() ,
		      []( double w ) { return( w == Inf< double >() ); } ) );
 }

/*--------------------------------------------------------------------------*/
/*------------- Methods for checking the state of the CSPBlock -------------*/
/*--------------------------------------------------------------------------*/

void CSPBlock::check_assignment( std::span< const int > a ,
				 const char * method ) const
{
 if( a.size() != get_number_variables() )
  throw( std::invalid_argument( std::string( "CSPBlock::" ) + method +
				": one value per variable is needed" ) );
 }

/*--------------------------------------------------------------------------*/

bool CSPBlock::is_satisfied( Index k , std::span< const int > a ) const
{
 check_assignment( a , "is_satisfied" );
 if( k >= get_number_constraints() )
  throw( std::invalid_argument( "CSPBlock::is_satisfied: constraint " +
				std::to_string( k ) + " does not exist" ) );
 const auto s = v_scope[ k ];
 const auto p = v_param[ k ];
 const auto n = s.size();

 switch( v_kind[ k ] ) {
  case( kExtension ):
  case( kConflicts ): {
   bool found = false;
   for( std::size_t t = 0 ; ( ! found ) && ( t < p.size() ) ; t += n ) {
    std::size_t i = 0;
    while( ( i < n ) && ( a[ s[ i ] ] == p[ t + i ] ) )
     ++i;
    found = ( i == n );
    }
   return( found == ( v_kind[ k ] == kExtension ) );
   }
  case( kLinear ): {
   long long lhs = 0;
   for( std::size_t i = 0 ; i < n ; ++i )
    lhs += ( long long ) p[ i ] * a[ s[ i ] ];
   const long long b = p[ n + 1 ];
   switch( p[ n ] ) {
    case( kLE ): return( lhs <= b );
    case( kEQ ): return( lhs == b );
    default:     return( lhs != b );
    }
   }
  case( kAllDifferent ): {
   v_int v;
   v.reserve( n );
   for( auto i : s )
    v.push_back( a[ i ] );
   std::sort( v.begin() , v.end() );
   return( std::adjacent_find( v.begin() , v.end() ) == v.end() );
   }
  case( kElement ): {
   const long long idx = ( long long ) a[ s[ 0 ] ] - p[ 0 ];
   if( ( idx < 0 ) || ( idx >= ( long long ) n - 2 ) )
    return( false );
   return( a[ s[ 2 + idx ] ] == a[ s[ 1 ] ] );
   }
  case( kCumulative ): {
   // the usage changes at the starts and at the ends of the tasks, the
   // ends coming first at the same time, as a task is over at s_j + d_j
   std::vector< std::pair< long long , long long > > ev;
   for( std::size_t j = 0 ; j < n ; ++j )
    if( p[ j ] && p[ n + j ] ) {
     ev.emplace_back( a[ s[ j ] ] , p[ n + j ] );
     ev.emplace_back( ( long long ) a[ s[ j ] ] + p[ j ] , - p[ n + j ] );
     }
   std::sort( ev.begin() , ev.end() );
   long long use = 0;
   for( const auto & e : ev )
    if( ( use += e.second ) > p[ 2 * n ] )
     return( false );
   return( true );
   }
  default:
   return( false );
  }
 }

/*--------------------------------------------------------------------------*/

bool CSPBlock::is_feasible( std::span< const int > a ) const
{
 check_assignment( a , "is_feasible" );
 for( Index i = 0 ; i < a.size() ; ++i )
  if( ! in_domain( i , a[ i ] ) )
   return( false );
 for( Index k = 0 ; k < get_number_constraints() ; ++k )
  if( v_active[ k ] && is_hard( k ) && ( ! is_satisfied( k , a ) ) )
   return( false );
 return( true );
 }

/*--------------------------------------------------------------------------*/

double CSPBlock::get_violated_weight( std::span< const int > a ) const
{
 check_assignment( a , "get_violated_weight" );
 double w = 0;
 for( Index k = 0 ; k < get_number_constraints() ; ++k )
  if( v_active[ k ] && ( ! is_hard( k ) ) && ( ! is_satisfied( k , a ) ) )
   w += v_weights[ k ];
 return( w );
 }

/*--------------------------------------------------------------------------*/

double CSPBlock::get_objective_value( std::span< const int > a ) const
{
 double c = get_violated_weight( a );
 for( Index i = 0 ; i < a.size() ; ++i ) {
  const auto s = value_index( i , a[ i ] );
  if( s == Inf< Index >() )
   throw( std::invalid_argument( "CSPBlock::get_objective_value: value " +
				 std::to_string( a[ i ] ) + " not in the "
				 "universe of variable " +
				 std::to_string( i ) ) );
  c += v_costs[ s ];
  }
 return( c );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- Methods for serializing ------------------------*/
/*--------------------------------------------------------------------------*/

void CSPBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );

 auto nv = group.addDim( "NumberVariables" , get_number_variables() );
 if( ! get_number_variables() )
  return;
 auto nval = group.addDim( "NumberValues" , get_number_values() );
 group.addVar( "Universe" , netCDF::NcInt() , nval ).putVar(
					     v_universe.values().data() );
 std::vector< unsigned int > start( v_universe.starts().begin() ,
				    v_universe.starts().end() - 1 );
 group.addVar( "UniverseStart" , netCDF::NcUint() , nv ).putVar(
							     start.data() );
 if( std::find( v_domain.begin() , v_domain.end() , 0 ) != v_domain.end() )
  group.addVar( "Domain" , netCDF::NcUbyte() , nval ).putVar(
							   v_domain.data() );
 if( has_costs() )
  group.addVar( "Costs" , netCDF::NcDouble() , nval ).putVar(
							   v_costs.data() );
 if( ! v_group.empty() )
  group.addVar( "VariableGroups" , netCDF::NcInt() , nv ).putVar(
							    v_group.data() );

 if( v_kind.empty() )
  return;
 auto nc = group.addDim( "NumberConstraints" , v_kind.size() );
 group.addVar( "ConstraintKind" , netCDF::NcInt() , nc ).putVar(
							     v_kind.data() );
 put_ragged( group , v_scope , "Scope" , "NumberScopeEntries" , nc );
 put_ragged( group , v_param , "Parameters" , "NumberParameterEntries" , nc );
 if( ! all_hard() )
  group.addVar( "Weights" , netCDF::NcDouble() , nc ).putVar(
							  v_weights.data() );
 if( std::find( v_active.begin() , v_active.end() , 0 ) != v_active.end() )
  group.addVar( "Active" , netCDF::NcUbyte() , nc ).putVar(
							   v_active.data() );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- Methods for printing --------------------------*/
/*--------------------------------------------------------------------------*/

void CSPBlock::print( std::ostream & output , char vlvl ) const
{
 const auto n = get_number_variables();
 if( vlvl != 'C' ) {
  output << "CSPBlock with " << n << " variables and "
	 << get_number_constraints() << " constraints" << std::endl;
  return;
  }

 // the costs are those of "minimize x_v" (sign 1) or "maximize x_v" (sign
 // -1) if they are 0 but for the variable v, whose value d costs sign * d
 int ov = -1 , os = 0;
 bool annotate = false;
 for( Index i = 0 ; ( ! annotate ) && ( i < n ) ; ++i ) {
  const auto u = get_universe( i );
  const auto first = v_universe.starts()[ i ];
  if( std::all_of( v_costs.begin() + first ,
		   v_costs.begin() + first + u.size() ,
		   []( double c ) { return( c == 0 ); } ) )
   continue;
  int sign = 0;
  for( int s : { 1 , -1 } ) {
   Index j = 0;
   while( ( j < u.size() ) &&
	  ( v_costs[ first + j ] == s * double( u[ j ] ) ) )
    ++j;
   if( j == u.size() )
    sign = s;
   }
  if( ( ov >= 0 ) || ( ! sign ) )
   annotate = true;  // more than one variable, or not a linear cost
  else {
   ov = int( i );
   os = sign;
   }
  }

 const auto prec = output.precision(
			     std::numeric_limits< double >::max_digits10 );
 output << "% CSPBlock: " << n << " variables, " << get_number_constraints()
	<< " constraints" << std::endl;
 for( Index i = 0 ; i < n ; ++i ) {
  const auto u = get_universe( i );
  const auto d = get_domain( i );
  output << "var ";
  print_set( output , d );
  output << ": x" << i;
  if( d.size() != u.size() ) {
   output << " :: universe(";
   print_set( output , u );
   output << ")";
   }
  if( ! v_group.empty() )
   output << " :: group(" << v_group[ i ] << ")";
  const auto first = v_universe.starts()[ i ];
  if( annotate &&
      std::any_of( v_costs.begin() + first ,
		   v_costs.begin() + first + u.size() ,
		   []( double c ) { return( c != 0 ); } ) ) {
   output << " :: cost([";
   for( Index j = 0 ; j < u.size() ; ++j )
    output << ( j ? "," : "" ) << v_costs[ first + j ];
   output << "])";
   }
  output << ";" << std::endl;
  }

 static const char * lin[] = { "int_lin_le" , "int_lin_eq" , "int_lin_ne" };
 for( Index k = 0 ; k < get_number_constraints() ; ++k ) {
  const auto s = v_scope[ k ];
  const auto p = v_param[ k ];
  const auto m = s.size();
  output << "constraint ";
  switch( v_kind[ k ] ) {
   case( kExtension ):
   case( kConflicts ):
    output << ( v_kind[ k ] == kExtension ? "fzn_table_int("
					  : "csp_table_int_neg(" );
    print_array( output , s , "x" );
    output << ",";
    print_array( output , p );
    break;
   case( kLinear ):
    output << lin[ p[ m ] ] << "(";
    print_array( output , p.first( m ) );
    output << ",";
    print_array( output , s , "x" );
    output << "," << p[ m + 1 ];
    break;
   case( kAllDifferent ):
    output << "fzn_all_different_int(";
    print_array( output , s , "x" );
    break;
   case( kElement ):
    if( p[ 0 ] == 1 )
     output << "array_var_int_element(x" << s[ 0 ] << ",";
    else
     output << "csp_element_offset(x" << s[ 0 ] << "," << p[ 0 ] << ",";
    print_array( output , s.subspan( 2 ) , "x" );
    output << ",x" << s[ 1 ];
    break;
   case( kCumulative ):
    output << "fzn_cumulative(";
    print_array( output , s , "x" );
    output << ",";
    print_array( output , p.first( m ) );
    output << ",";
    print_array( output , p.subspan( m , m ) );
    output << "," << p[ 2 * m ];
    break;
   }
  output << ")";
  if( ! is_hard( k ) )
   output << " :: weight(" << v_weights[ k ] << ")";
  if( ! v_active[ k ] )
   output << " :: inactive";
  output << ";" << std::endl;
  }

 if( ( ov >= 0 ) && ( ! annotate ) )
  output << "solve " << ( os > 0 ? "minimize" : "maximize" ) << " x" << ov
	 << ";" << std::endl;
 else
  output << "solve satisfy;" << std::endl;
 output.precision( prec );
 }

/*--------------------------------------------------------------------------*/
/*----------- METHODS FOR MODIFYING THE PHYSICAL REPRESENTATION ------------*/
/*--------------------------------------------------------------------------*/

void CSPBlock::chg_domains( Subset && nms , v_v_int && domains , bool ordered ,
			    ModParam issueMod )
{
 if( nms.empty() )  // nothing to change
  return;
 if( domains.size() < nms.size() )
  throw( std::invalid_argument( "CSPBlock::chg_domains: fewer domains than "
				"variables" ) );
 domains.resize( nms.size() );
 if( ! ordered )
  sort_along( nms , domains );
 check_subset( nms , get_number_variables() , "CSPBlock::chg_domains" ,
	       "variable" );

 // the new flags of each variable, all checked before anything changes;
 // in a dry run they have been changed already, and the type says all
 // that may have happened to all the variables
 const bool dry = ! not_dry_run( issueMod );
 Modification::ModConcern c = 0;
 Subset changed;
 std::vector< v_flag > nf;
 for( Index k = 0 ; k < nms.size() ; ++k ) {
  const auto i = nms[ k ];
  const auto first = v_universe.starts()[ i ];
  v_flag f( get_universe( i ).size() , 0 );
  for( auto d : domains[ k ] ) {
   const auto s = value_index( i , d );
   if( s == Inf< Index >() )
    throw( std::invalid_argument( "CSPBlock::chg_domains: value " +
				  std::to_string( d ) + " not in the "
				  "universe of variable " +
				  std::to_string( i ) ) );
   f[ s - first ] = 1;
   }
  bool any = false;
  for( Index j = 0 ; j < f.size() ; ++j )
   if( f[ j ] != v_domain[ first + j ] ) {
    any = true;
    c |= f[ j ] ? Modification::eRegnGrow : Modification::eRegnShrink;
    }
  if( any || dry ) {
   changed.push_back( i );
   nf.push_back( std::move( f ) );
   }
  }
 if( changed.empty() )
  return;  // nothing changes, avoid issuing the Modification

 if( dry )
  c = CSPBlockMod::changes_of( CSPBlockMod::eChgDomain );
 else
  for( Index k = 0 ; k < changed.size() ; ++k )
   std::copy( nf[ k ].begin() , nf[ k ].end() ,
	      v_domain.begin() + v_universe.starts()[ changed[ k ] ] );

 c |= Modification::eModPhys;
 if( issue_pmod( issueMod , c ) ) {
  auto mod = std::make_shared< CSPBlockSbstMod >( this ,
			     CSPBlockMod::eChgDomain , std::move( changed ) );
  mod->set_changes( c );
  add_Modification( mod , Observer::par2chnl( issueMod ) );
  }
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::chg_domains( v_v_int && domains , Range rng ,
			    ModParam issueMod )
{
 rng.second = std::min( rng.second , Index( get_number_variables() ) );
 if( rng.second <= rng.first )  // nothing to change
  return;
 Subset nms( rng.second - rng.first );
 std::iota( nms.begin() , nms.end() , rng.first );
 chg_domains( std::move( nms ) , std::move( domains ) , true , issueMod );
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::add_values( Subset && nms , v_v_int && values ,
			   ModParam issueMod )
{
 if( nms.empty() )  // nothing to change
  return;
 if( values.size() < nms.size() )
  throw( std::invalid_argument( "CSPBlock::add_values: fewer sets of values "
				"than variables" ) );
 values.resize( nms.size() );
 sort_along( nms , values );
 check_subset( nms , get_number_variables() , "CSPBlock::add_values" ,
	       "variable" );

 // the values really new of each variable; in a dry run they have been
 // added already, and all the variables are said to have changed
 const bool dry = ! not_dry_run( issueMod );
 Subset changed;
 std::vector< v_int > nv( get_number_variables() );
 for( Index k = 0 ; k < nms.size() ; ++k ) {
  const auto i = nms[ k ];
  auto & v = values[ k ];
  sort_unique( v );
  const auto u = get_universe( i );
  std::set_difference( v.begin() , v.end() , u.begin() , u.end() ,
		       std::back_inserter( nv[ i ] ) );
  if( dry || ( ! nv[ i ].empty() ) )
   changed.push_back( i );
  }
 if( changed.empty() )
  return;  // nothing changes, avoid issuing the Modification

 if( ! dry ) {
  // the universes merged, and the domains and the costs along with them
  v_int val;
  std::vector< std::size_t > start = { 0 };
  v_flag dom;
  v_dbl cost;
  for( Index i = 0 ; i < get_number_variables() ; ++i ) {
   const auto u = get_universe( i );
   auto first = v_universe.starts()[ i ];
   std::size_t a = 0 , b = 0;
   while( ( a < u.size() ) || ( b < nv[ i ].size() ) )
    if( ( b == nv[ i ].size() ) ||
	( ( a < u.size() ) && ( u[ a ] < nv[ i ][ b ] ) ) ) {
     val.push_back( u[ a ] );
     dom.push_back( v_domain[ first + a ] );
     cost.push_back( v_costs[ first + a ] );
     ++a;
     }
    else {
     val.push_back( nv[ i ][ b++ ] );
     dom.push_back( 1 );
     cost.push_back( 0 );
     }
   start.push_back( val.size() );
   }
  v_universe.assign( std::move( val ) , std::move( start ) );
  v_domain = std::move( dom );
  v_costs = std::move( cost );
  }

 if( issue_pmod( issueMod , CSPBlockMod::changes_of(
				  CSPBlockMod::eAddValues ) ) )
  add_Modification( std::make_shared< CSPBlockSbstMod >( this ,
			     CSPBlockMod::eAddValues , std::move( changed ) ) ,
		    Observer::par2chnl( issueMod ) );
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::chg_costs( MF_dbl_sp NCost , Range rng , ModParam issueMod )
{
 rng.second = std::min( rng.second , Index( get_number_values() ) );
 if( rng.second <= rng.first )  // nothing to change
  return;
 if( NCost.size() < rng.second - rng.first )
  throw( std::invalid_argument( "CSPBlock::chg_costs: the span is shorter "
				"than the Range" ) );
 for( Index i = rng.first ; i < rng.second ; ++i )
  if( ! std::isfinite( NCost[ i - rng.first ] ) )
   throw( std::invalid_argument( "CSPBlock::chg_costs: cost " +
				 std::to_string( NCost[ i - rng.first ] ) +
				 " of slot " + std::to_string( i ) ) );
 if( not_dry_run( issueMod ) ) {
  if( std::equal( NCost.begin() , NCost.begin() + ( rng.second - rng.first ) ,
		  v_costs.begin() + rng.first ) )
   return;  // nothing changes, avoid issuing the Modification
  std::copy( NCost.begin() , NCost.begin() + ( rng.second - rng.first ) ,
	     v_costs.begin() + rng.first );
  }

 if( issue_pmod( issueMod , CSPBlockMod::changes_of(
				  CSPBlockMod::eChgCost ) ) )
  add_Modification( std::make_shared< CSPBlockRngdMod >( this ,
					   CSPBlockMod::eChgCost , rng ) ,
		    Observer::par2chnl( issueMod ) );
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::chg_costs( MF_dbl_sp NCost , Subset && nms , bool ordered ,
			  ModParam issueMod )
{
 if( nms.empty() )  // nothing to change
  return;
 if( NCost.size() < nms.size() )
  throw( std::invalid_argument( "CSPBlock::chg_costs: the span is shorter "
				"than the Subset" ) );
 v_dbl nc( NCost.begin() , NCost.begin() + nms.size() );
 if( ! ordered )
  sort_along( nms , nc );
 check_subset( nms , get_number_values() , "CSPBlock::chg_costs" , "slot" );
 for( Index k = 0 ; k < nms.size() ; ++k )
  if( ! std::isfinite( nc[ k ] ) )
   throw( std::invalid_argument( "CSPBlock::chg_costs: cost " +
				 std::to_string( nc[ k ] ) + " of slot " +
				 std::to_string( nms[ k ] ) ) );
 if( not_dry_run( issueMod ) ) {
  bool any = false;
  for( Index k = 0 ; k < nms.size() ; ++k )
   if( v_costs[ nms[ k ] ] != nc[ k ] ) {
    v_costs[ nms[ k ] ] = nc[ k ];
    any = true;
    }
  if( ! any )
   return;  // nothing changes, avoid issuing the Modification
  }

 if( issue_pmod( issueMod , CSPBlockMod::changes_of(
				  CSPBlockMod::eChgCost ) ) )
  add_Modification( std::make_shared< CSPBlockSbstMod >( this ,
				   CSPBlockMod::eChgCost , std::move( nms ) ) ,
		    Observer::par2chnl( issueMod ) );
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::chg_weights( MF_dbl_sp NWeight , Range rng ,
			    ModParam issueMod )
{
 rng.second = std::min( rng.second , Index( get_number_constraints() ) );
 if( rng.second <= rng.first )  // nothing to change
  return;
 Subset nms( rng.second - rng.first );
 std::iota( nms.begin() , nms.end() , rng.first );
 if( NWeight.size() < nms.size() )
  throw( std::invalid_argument( "CSPBlock::chg_weights: the span is "
				"shorter than the Range" ) );
 v_dbl nw( NWeight.begin() , NWeight.begin() + nms.size() );
 check_weights( nw );

 auto c = CSPBlockMod::changes_of( CSPBlockMod::eChgWeight );
 if( not_dry_run( issueMod ) ) {
  if( std::equal( nw.begin() , nw.end() , v_weights.begin() + rng.first ) )
   return;  // nothing changes, avoid issuing the Modification
  c = 0;
  for( Index k = 0 ; k < nms.size() ; ++k ) {
   c |= weight_changes( v_weights[ nms[ k ] ] , nw[ k ] );
   v_weights[ nms[ k ] ] = nw[ k ];
   }
  }

 if( issue_pmod( issueMod , c ) ) {
  auto mod = std::make_shared< CSPBlockRngdMod >( this ,
					   CSPBlockMod::eChgWeight , rng );
  mod->set_changes( c );
  add_Modification( mod , Observer::par2chnl( issueMod ) );
  }
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::chg_weights( MF_dbl_sp NWeight , Subset && nms , bool ordered ,
			    ModParam issueMod )
{
 if( nms.empty() )  // nothing to change
  return;
 if( NWeight.size() < nms.size() )
  throw( std::invalid_argument( "CSPBlock::chg_weights: the span is "
				"shorter than the Subset" ) );
 v_dbl nw( NWeight.begin() , NWeight.begin() + nms.size() );
 if( ! ordered )
  sort_along( nms , nw );
 check_subset( nms , get_number_constraints() , "CSPBlock::chg_weights" ,
	       "constraint" );
 check_weights( nw );

 auto c = CSPBlockMod::changes_of( CSPBlockMod::eChgWeight );
 if( not_dry_run( issueMod ) ) {
  c = 0;
  bool any = false;
  for( Index k = 0 ; k < nms.size() ; ++k )
   if( v_weights[ nms[ k ] ] != nw[ k ] ) {
    c |= weight_changes( v_weights[ nms[ k ] ] , nw[ k ] );
    v_weights[ nms[ k ] ] = nw[ k ];
    any = true;
    }
  if( ! any )
   return;  // nothing changes, avoid issuing the Modification
  }

 if( issue_pmod( issueMod , c ) ) {
  auto mod = std::make_shared< CSPBlockSbstMod >( this ,
				  CSPBlockMod::eChgWeight , std::move( nms ) );
  mod->set_changes( c );
  add_Modification( mod , Observer::par2chnl( issueMod ) );
  }
 }

/*--------------------------------------------------------------------------*/

Modification::ModConcern CSPBlock::set_active( Index k , bool active )
{
 if( bool( v_active[ k ] ) == active )
  return( 0 );
 v_active[ k ] = active;
 if( is_hard( k ) )
  return( Modification::eModPhys | Modification::eModCnsSide |
	  ( active ? Modification::eRegnShrink : Modification::eRegnGrow ) );
 return( Modification::eModPhys | Modification::eModCnsSide |
	 ( active ? Modification::eObjUp : Modification::eObjDown ) );
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::chg_active( bool active , Range rng , ModParam issueMod )
{
 rng.second = std::min( rng.second , Index( get_number_constraints() ) );
 if( rng.second <= rng.first )  // nothing to change
  return;

 auto c = CSPBlockMod::changes_of( CSPBlockMod::eChgActive );
 if( not_dry_run( issueMod ) ) {
  c = 0;
  for( Index k = rng.first ; k < rng.second ; ++k )
   c |= set_active( k , active );
  if( ! c )
   return;  // nothing changes, avoid issuing the Modification
  }

 if( issue_pmod( issueMod , c ) ) {
  auto mod = std::make_shared< CSPBlockRngdMod >( this ,
					   CSPBlockMod::eChgActive , rng );
  mod->set_changes( c );
  add_Modification( mod , Observer::par2chnl( issueMod ) );
  }
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::chg_active( bool active , Subset && nms , bool ordered ,
			   ModParam issueMod )
{
 if( nms.empty() )  // nothing to change
  return;
 if( ! ordered )
  std::sort( nms.begin() , nms.end() );
 check_subset( nms , get_number_constraints() , "CSPBlock::chg_active" ,
	       "constraint" );

 auto c = CSPBlockMod::changes_of( CSPBlockMod::eChgActive );
 if( not_dry_run( issueMod ) ) {
  c = 0;
  for( auto k : nms )
   c |= set_active( k , active );
  if( ! c )
   return;  // nothing changes, avoid issuing the Modification
  }

 if( issue_pmod( issueMod , c ) ) {
  auto mod = std::make_shared< CSPBlockSbstMod >( this ,
				  CSPBlockMod::eChgActive , std::move( nms ) );
  mod->set_changes( c );
  add_Modification( mod , Observer::par2chnl( issueMod ) );
  }
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::add_constraints( v_int && kinds , v_v_int && scopes ,
				v_v_int && params , v_dbl && weights ,
				ModParam issueMod )
{
 if( kinds.empty() )  // nothing to add
  return;
 if( ( scopes.size() != kinds.size() ) || ( params.size() != kinds.size() ) )
  throw( std::invalid_argument( "CSPBlock::add_constraints: kinds, scopes "
				"and parameters must have the same size" ) );
 if( weights.empty() )
  weights.assign( kinds.size() , Inf< double >() );
 if( weights.size() != kinds.size() )
  throw( std::invalid_argument( "CSPBlock::add_constraints: one weight per "
				"constraint is needed" ) );
 for( Index c = 0 ; c < kinds.size() ; ++c )
  check_constraint( kinds[ c ] , scopes[ c ] , params[ c ] );
 check_weights( weights );

 const Range rng( get_number_constraints() ,
		  get_number_constraints() + kinds.size() );
 if( not_dry_run( issueMod ) ) {
  for( Index c = 0 ; c < kinds.size() ; ++c ) {
   v_scope.push_back( scopes[ c ] );
   v_param.push_back( params[ c ] );
   }
  v_kind.insert( v_kind.end() , kinds.begin() , kinds.end() );
  v_weights.insert( v_weights.end() , weights.begin() , weights.end() );
  v_active.resize( v_kind.size() , 1 );
  }

 if( issue_pmod( issueMod , CSPBlockMod::changes_of(
				  CSPBlockMod::eAddConstraints ) ) )
  add_Modification( std::make_shared< CSPBlockRngdMod >( this ,
				   CSPBlockMod::eAddConstraints , rng ) ,
		    Observer::par2chnl( issueMod ) );
 }

/*--------------------------------------------------------------------------*/

void CSPBlock::add_variables( v_v_int && universes , v_int && groups ,
			      ModParam issueMod )
{
 if( universes.empty() )  // nothing to add
  return;
 for( auto & u : universes ) {
  sort_unique( u );
  if( u.empty() )
   throw( std::invalid_argument( "CSPBlock::add_variables: empty "
				 "universe" ) );
  }
 if( v_group.empty() != groups.empty() )
  throw( std::invalid_argument( "CSPBlock::add_variables: groups are "
				"needed if and only if the variables have "
				"them" ) );
 if( ( ! groups.empty() ) &&
     ( ( groups.size() != universes.size() ) ||
       std::any_of( groups.begin() , groups.end() ,
		    []( int g ) { return( g < 0 ); } ) ) )
  throw( std::invalid_argument( "CSPBlock::add_variables: one non-negative "
				"group per variable is needed" ) );

 const Range rng( get_number_variables() ,
		  get_number_variables() + universes.size() );
 if( not_dry_run( issueMod ) ) {
  for( const auto & u : universes ) {
   v_universe.push_back( u );
   v_domain.resize( v_domain.size() + u.size() , 1 );
   v_costs.resize( v_costs.size() + u.size() , 0 );
   }
  v_group.insert( v_group.end() , groups.begin() , groups.end() );
  }

 if( issue_pmod( issueMod , CSPBlockMod::changes_of(
				  CSPBlockMod::eAddVariables ) ) )
  add_Modification( std::make_shared< CSPBlockRngdMod >( this ,
				    CSPBlockMod::eAddVariables , rng ) ,
		    Observer::par2chnl( issueMod ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF CSPBlockChange ------------------------*/
/*--------------------------------------------------------------------------*/

Change * CSPBlockChange::apply( Block * block , bool doUndo ,
				ModParam issueMod , ModParam )
{
 auto csp = dynamic_cast< CSPBlock * >( block );
 if( ! csp )
  throw( std::invalid_argument( "CSPBlockChange::apply: the Block is not a "
				"CSPBlock" ) );
 if( ( f_type < eAssign ) || ( f_type >= eEmpty ) )
  throw( std::invalid_argument( "CSPBlockChange::apply: empty "
				"CSPBlockChange" ) );
 const bool set = ( f_type == eSetDomain );
 if( set ? ( ( f_sizes.size() != f_nms.size() ) ||
	     ( std::accumulate( f_sizes.begin() , f_sizes.end() ,
				std::size_t( 0 ) ) != f_values.size() ) )
	 : ( f_values.size() != f_nms.size() ) )
  throw( std::invalid_argument( "CSPBlockChange::apply: the values do not "
				"match the variables" ) );
 for( auto i : f_nms )
  if( i >= csp->get_number_variables() )
   throw( std::invalid_argument( "CSPBlockChange::apply: variable " +
				 std::to_string( i ) + " out of range" ) );

 // the new domains, and the old ones for the undo
 CSPBlock::v_v_int nd , od;
 std::size_t pos = 0;
 for( std::size_t k = 0 ; k < f_nms.size() ; ++k ) {
  auto d = csp->get_domain( f_nms[ k ] );
  CSPBlock::v_int n;
  if( set ) {
   n.assign( f_values.begin() + pos , f_values.begin() + pos + f_sizes[ k ] );
   pos += f_sizes[ k ];
   }
  else {
   const int v = f_values[ k ];
   for( auto x : d )
    if( ( f_type == eAssign ) ? ( x == v ) :
	( f_type == eRemove ) ? ( x != v ) :
	( f_type == eSplitLE ) ? ( x <= v ) : ( x > v ) )
     n.push_back( x );
   }
  nd.push_back( std::move( n ) );
  if( doUndo )
   od.push_back( std::move( d ) );
  }

 Change * undo = nullptr;
 if( doUndo ) {
  std::vector< int > val;
  std::vector< unsigned int > sz;
  for( const auto & d : od ) {
   val.insert( val.end() , d.begin() , d.end() );
   sz.push_back( d.size() );
   }
  undo = new CSPBlockChange( eSetDomain , Block::Subset( f_nms ) ,
			     std::move( val ) , std::move( sz ) );
  }

 try {
  csp->chg_domains( Block::Subset( f_nms ) , std::move( nd ) , false ,
		    issueMod );
  }
 catch( ... ) {
  delete undo;
  throw;
  }
 return( undo );
 }

/*--------------------------------------------------------------------------*/

void CSPBlockChange::deserialize( const netCDF::NcGroup & group )
{
 auto ftype = group.getAtt( "CSPBlockChange_type" );
 if( ftype.isNull() )
  throw( std::invalid_argument( "CSPBlockChange::deserialize: no type" ) );
 ftype.getValues( & f_type );

 f_nms.clear();
 f_values.clear();
 f_sizes.clear();
 const auto nv = group.getDim( "NumVar" );
 if( nv.isNull() || ( ! nv.getSize() ) )
  return;
 f_nms.resize( nv.getSize() );
 group.getVar( "Index" ).getVar( f_nms.data() );
 const auto val = group.getVar( "Value" );
 if( ! val.isNull() ) {
  f_values.resize( val.getDim( 0 ).getSize() );
  if( ! f_values.empty() )
   val.getVar( f_values.data() );
  }
 if( f_type == eSetDomain ) {
  std::vector< unsigned int > start( f_nms.size() , 0 );
  const auto st = group.getVar( "ValueStart" );
  if( ! st.isNull() )
   st.getVar( start.data() );
  f_sizes.resize( f_nms.size() );
  for( std::size_t k = 0 ; k < f_nms.size() ; ++k ) {
   const auto end = ( k + 1 < f_nms.size() ) ? start[ k + 1 ]
					      : f_values.size();
   if( ( end < start[ k ] ) || ( end > f_values.size() ) )
    throw( std::invalid_argument( "CSPBlockChange::deserialize: wrong "
				  "ValueStart" ) );
   f_sizes[ k ] = end - start[ k ];
   }
  }
 }

/*--------------------------------------------------------------------------*/

void CSPBlockChange::serialize( netCDF::NcGroup & group ) const
{
 Change::serialize( group );
 group.putAtt( "CSPBlockChange_type" , netCDF::NcInt() , f_type );
 if( f_nms.empty() )
  return;
 const auto nv = group.addDim( "NumVar" , f_nms.size() );
 group.addVar( "Index" , netCDF::NcUint() , nv ).putVar( f_nms.data() );
 if( ! f_values.empty() ) {
  const auto nval = group.addDim( "NumValues" , f_values.size() );
  group.addVar( "Value" , netCDF::NcInt() , nval ).putVar( f_values.data() );
  }
 if( f_type == eSetDomain ) {
  std::vector< unsigned int > start( f_nms.size() , 0 );
  for( std::size_t k = 1 ; k < f_nms.size() ; ++k )
   start[ k ] = start[ k - 1 ] + f_sizes[ k - 1 ];
  group.addVar( "ValueStart" , netCDF::NcUint() , nv ).putVar(
							     start.data() );
  }
 }

/*--------------------------------------------------------------------------*/
/*----------- Methods for the abstract representation (phase 1) ------------*/
/*--------------------------------------------------------------------------*/

void CSPBlock::generate_abstract_variables( Configuration * stvv )
{
 if( f_AR & HasVar )  // the variables are there already
  return;             // nothing to do

 const auto n = get_number_variables();
 v_var.resize( n );   // resized once, so the ColVariable do not move
 for( Index i = 0 ; i < n ; ++i )
  v_var[ i ].is_integer( true , eNoBlck );  // integer; bounds are the box

 add_static_variable( v_var , "x" );

 f_AR |= HasVar;

 }  // end( CSPBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void CSPBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( ! ( f_AR & HasVar ) )  // the Variable are built first, as the protocol
  generate_abstract_variables();  // asks [see Block.h]

 if( f_AR & HasCon )  // the constraints are there already
  return;             // nothing to do

 const auto n = get_number_variables();
 const auto m = get_number_constraints();

 // the box of each variable, i.e., the min..max of its current domain; a
 // domain with holes gives the same box, its missing values being left to
 // the physical current domain (TODO: enforce them in the abstract form too,
 // e.g., with a small table), and an empty domain gives an infeasible box
 v_box.resize( n );
 for( Index i = 0 ; i < n ; ++i ) {
  const auto d = get_domain( i );
  v_box[ i ].set_variable( & v_var[ i ] , eNoBlck );
  if( d.empty() ) {           // no value is left: an infeasible box
   v_box[ i ].set_lhs( 1 , eNoBlck );
   v_box[ i ].set_rhs( 0 , eNoBlck );
   }
  else {
   v_box[ i ].set_lhs( d.front() , eNoBlck );
   v_box[ i ].set_rhs( d.back() , eNoBlck );
   }
  }
 add_static_constraint( v_box , "bnd" );

 // the distinct variables of the scope of the k-th constraint, in the order
 // of their first appearance
 const auto distinct_scope = [ this ]( Index k ) {
  CSPConstraint::v_var_ptr vars;
  for( auto s : get_scopes()[ k ] )
   if( std::find( vars.begin() , vars.end() , & v_var[ s ] ) == vars.end() )
    vars.push_back( & v_var[ s ] );
  return( vars );
  };

 // an active and hard constraint is in exactly one of the groups below, by
 // its kind, and, for a kLinear, by its operator; a soft or inactive one is
 // not in the abstract representation [see the class comment]. Count first,
 // so each vector is resized once and its :Constraint never move
 Index n_lin = 0 , n_ne = 0 , n_ad = 0 , n_tb = 0 , n_el = 0 , n_cu = 0;
 for( Index k = 0 ; k < m ; ++k ) {
  if( ( ! v_active[ k ] ) || ( ! is_hard( k ) ) )
   continue;
  switch( v_kind[ k ] ) {
   case( kLinear ): {
    const auto p = get_parameters()[ k ];
    ( p[ p.size() - 2 ] == kNE ? n_ne : n_lin )++;
    break;
    }
   case( kAllDifferent ): ++n_ad; break;
   case( kExtension ):
   case( kConflicts ):    ++n_tb; break;
   case( kElement ):      ++n_el; break;
   case( kCumulative ):   ++n_cu; break;
   default: break;
   }
  }

 v_lin.resize( n_lin );
 v_ne.resize( n_ne );
 v_alldiff.resize( n_ad );
 v_table.resize( n_tb );
 v_element.resize( n_el );
 v_cumul.resize( n_cu );

 // fill the groups; the counters now index into the resized vectors
 Index i_lin = 0 , i_ne = 0 , i_ad = 0 , i_tb = 0 , i_el = 0 , i_cu = 0;
 for( Index k = 0 ; k < m ; ++k ) {
  if( ( ! v_active[ k ] ) || ( ! is_hard( k ) ) )
   continue;

  switch( v_kind[ k ] ) {
   case( kLinear ): {
    const auto s = get_scopes()[ k ];
    const auto p = get_parameters()[ k ];
    const int op = p[ p.size() - 2 ];
    const long long b = p[ p.size() - 1 ];

    if( op == kNE ) {       // not a single row: an abstract verifier only
     v_ne[ i_ne ].set_data( this , k , distinct_scope( k ) );
     ++i_ne;
     break;
     }

    // kLE or kEQ: an FRowConstraint over a LinearFunction, the coefficients
    // of a repeated variable summed so that each variable appears once
    std::map< Index , double > agg;
    for( Index j = 0 ; j < s.size() ; ++j )
     agg[ s[ j ] ] += double( p[ j ] );
    LinearFunction::v_coeff_pair pairs;
    pairs.reserve( agg.size() );
    for( const auto & [ idx , coef ] : agg )
     if( coef != 0 )
      pairs.emplace_back( & v_var[ idx ] , coef );

    auto & fr = v_lin[ i_lin ];
    if( op == kEQ )
     fr.set_both( double( b ) , eNoBlck );
    else {                  // kLE: - INF <= sum <= b
     fr.set_lhs( - RowConstraint::RHSINF , eNoBlck );
     fr.set_rhs( double( b ) , eNoBlck );
     }
    fr.set_function( new LinearFunction( std::move( pairs ) ) , eNoBlck );
    ++i_lin;
    break;
    }
   case( kAllDifferent ):
    v_alldiff[ i_ad ].set_data( this , k , distinct_scope( k ) );
    ++i_ad;
    break;
   case( kExtension ):
   case( kConflicts ):
    v_table[ i_tb ].set_data( this , k , distinct_scope( k ) );
    ++i_tb;
    break;
   case( kElement ):
    v_element[ i_el ].set_data( this , k , distinct_scope( k ) );
    ++i_el;
    break;
   case( kCumulative ):
    v_cumul[ i_cu ].set_data( this , k , distinct_scope( k ) );
    ++i_cu;
    break;
   default: break;
   }
  }

 if( ! v_lin.empty() )     add_static_constraint( v_lin , "lin" );
 if( ! v_ne.empty() )      add_static_constraint( v_ne , "ne" );
 if( ! v_alldiff.empty() ) add_static_constraint( v_alldiff , "alldiff" );
 if( ! v_table.empty() )   add_static_constraint( v_table , "table" );
 if( ! v_element.empty() ) add_static_constraint( v_element , "element" );
 if( ! v_cumul.empty() )   add_static_constraint( v_cumul , "cumulative" );

 f_AR |= HasCon;

 }  // end( CSPBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

bool CSPBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 if( ! useabstract )
  return( Block::is_feasible( false , fsbc ) );

 if( ! ( f_AR & HasCon ) )
  throw( std::logic_error( "CSPBlock::is_feasible: the abstract constraints "
			   "are not generated" ) );

 const auto chk = []( auto & vec ) {
  for( auto & c : vec ) {
   c.compute();
   if( ! c.feasible() )
    return( false );
   }
  return( true );
  };

 return( chk( v_box ) && chk( v_lin ) && chk( v_ne ) && chk( v_alldiff ) &&
	 chk( v_table ) && chk( v_element ) && chk( v_cumul ) );
 }

/*--------------------------------------------------------------------------*/
/*------ Methods of the abstract global constraints [see the header] -------*/
/*--------------------------------------------------------------------------*/

int CSPConstraint::compute( bool changedvars )
{
 const auto n = f_csp->get_number_variables();
 auto & var = f_csp->get_abstract_variables();
 CSPBlock::v_int a( n );
 for( Block::Index i = 0 ; i < n ; ++i )
  a[ i ] = int( std::lround( var[ i ].get_value() ) );
 f_feasible = f_csp->is_satisfied( f_k , a );
 return( kOK );
 }

/*--------------------------------------------------------------------------*/

std::list< FRowConstraint * > AllDifferentConstraint::lin_relaxation( void )
 const
{
 const auto k = f_var.size();
 if( k == 0 )
  return( std::list< FRowConstraint * >() );

 // the valid "sum" bound (Williams and Yan): the k distinct variables take
 // distinct values, so their sum is at least the sum of the k smallest and
 // at most the sum of the k largest distinct values any of them may take,
 // read off the union of the current domains of the scope
 CSPBlock::v_int vals;
 CSPBlock::v_int seen;
 for( auto s : f_csp->get_scopes()[ f_k ] ) {
  if( std::find( seen.begin() , seen.end() , s ) != seen.end() )
   continue;
  seen.push_back( s );
  const auto d = f_csp->get_domain( s );
  vals.insert( vals.end() , d.begin() , d.end() );
  }
 std::sort( vals.begin() , vals.end() );
 vals.erase( std::unique( vals.begin() , vals.end() ) , vals.end() );
 if( vals.size() < k ) {
  // fewer distinct values than variables: by the pigeonhole principle no
  // assignment of pairwise distinct values exists, so the constraint, and
  // thus its relaxation, is infeasible; witness it with a row that no point
  // satisfies (1 <= sum x_i <= 0) rather than dropping the infeasibility
  LinearFunction::v_coeff_pair pairs;
  pairs.reserve( k );
  for( auto v : f_var )
   pairs.emplace_back( v , double( 1 ) );
  auto * fr = new FRowConstraint();
  fr->set_lhs( double( 1 ) , eNoBlck );
  fr->set_rhs( double( 0 ) , eNoBlck );
  fr->set_function( new LinearFunction( std::move( pairs ) ) , eNoBlck );
  return( std::list< FRowConstraint * >{ fr } );
  }

 long long lo = 0 , hi = 0;
 for( std::size_t j = 0 ; j < k ; ++j ) {
  lo += vals[ j ];
  hi += vals[ vals.size() - 1 - j ];
  }

 LinearFunction::v_coeff_pair pairs;
 pairs.reserve( k );
 for( auto v : f_var )
  pairs.emplace_back( v , double( 1 ) );

 auto * fr = new FRowConstraint();
 fr->set_lhs( double( lo ) , eNoBlck );
 fr->set_rhs( double( hi ) , eNoBlck );
 fr->set_function( new LinearFunction( std::move( pairs ) ) , eNoBlck );
 return( std::list< FRowConstraint * >{ fr } );
 }

/*--------------------------------------------------------------------------*/

std::list< FRowConstraint * > TableConstraint::lin_relaxation( void ) const
{
 std::list< FRowConstraint * > rel;

 // only a positive table (kExtension) admits the cheap per-variable box; a
 // negative one (kConflicts) has a non-convex feasible region [see header].
 // The kind is read from the CSPBlock, since the stored index alone does not
 // tell the two apart
 if( ( f_k >= f_csp->get_number_constraints() ) ||
     ( f_csp->get_kinds()[ f_k ] != CSPBlock::kExtension ) )
  return( rel );

 const auto s = f_csp->get_scopes()[ f_k ];
 const auto p = f_csp->get_parameters()[ f_k ];
 const auto n = s.size();
 if( ( n == 0 ) || p.empty() )  // no column or no tuple: nothing to bound
  return( rel );

 // for each distinct variable of the scope, the min and the max of the values
 // it takes across the tuples, which are p read n values at a time [see
 // CSPBlock.h]; a variable at more than one position is bounded by the
 // min..max over all of them, which is a valid necessary condition
 auto & var = f_csp->get_abstract_variables();
 CSPBlock::v_int seen;
 for( std::size_t j = 0 ; j < n ; ++j ) {
  if( std::find( seen.begin() , seen.end() , s[ j ] ) != seen.end() )
   continue;
  seen.push_back( s[ j ] );
  long long lo = p[ j ] , hi = p[ j ];
  for( std::size_t jj = j ; jj < n ; ++jj )
   if( s[ jj ] == s[ j ] )
    for( std::size_t t = 0 ; t < p.size() ; t += n ) {
     lo = std::min< long long >( lo , p[ t + jj ] );
     hi = std::max< long long >( hi , p[ t + jj ] );
     }
  LinearFunction::v_coeff_pair pairs;
  pairs.emplace_back( & var[ s[ j ] ] , double( 1 ) );
  auto * fr = new FRowConstraint();
  fr->set_lhs( double( lo ) , eNoBlck );
  fr->set_rhs( double( hi ) , eNoBlck );
  fr->set_function( new LinearFunction( std::move( pairs ) ) , eNoBlck );
  rel.push_back( fr );
  }
 return( rel );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File CSPBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/

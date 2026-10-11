/*--------------------------------------------------------------------------*/
/*------------------------ File BranchingRule.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the components that a RelaxationSolver can use to
 * choose the variable to branch upon: BranchingFeatures, StrongBranching and
 * LearnedBranching.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <filesystem>

#include <fstream>

#include <iomanip>

#include <limits>

#include <stdexcept>

#include "BranchingRule.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF BranchingFeatures ---------------------*/
/*--------------------------------------------------------------------------*/

double BranchingFeatures::fractionality( double value )
{
 return( std::abs( value - std::round( value ) ) );
 }

/*--------------------------------------------------------------------------*/

void BranchingFeatures::sort_by_fractionality(
				      std::vector< BranchCandidate > & cands )
{
 std::stable_sort( cands.begin() , cands.end() ,
		   []( const BranchCandidate & a , const BranchCandidate & b ) {
		    return( fractionality( a.value ) >
			    fractionality( b.value ) ); } );
 }

/*--------------------------------------------------------------------------*/

std::vector< BranchingFeatures::Features > BranchingFeatures::features(
				const std::vector< BranchCandidate > & cands ,
				double fixed )
{
 double maxc = 0;
 for( const auto & cd : cands )
  maxc = std::max( maxc , std::abs( cd.cost ) );
 if( maxc <= 0 )
  maxc = 1;
 const std::size_t n = cands.size();

 std::vector< Features > x( n );
 for( std::size_t c = 0 ; c < n ; ++c ) {
  const auto & cd = cands[ c ];
  x[ c ] = { fractionality( cd.value ) , cd.value , cd.second ,
	     std::abs( cd.value - cd.second ) , cd.cost / maxc ,
	     cd.position , fixed , n > 1 ? double( c ) / ( n - 1 ) : 0.0 };
  }
 return( x );

 }  // end( BranchingFeatures::features )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS OF StrongBranching ----------------------*/
/*--------------------------------------------------------------------------*/

std::vector< StrongBranching::Eval > StrongBranching::evaluate(
				const std::vector< BranchCandidate > & cands ,
				const Evaluator & children , double node ,
				bool maximize , std::ostream * log )
{
 // the improvement of a child bound over that of the node, on the side of
 // the relaxation; at least eps, so that a candidate with one child not
 // improving still ranks by the other one
 const double eps = 1e-6 * std::max( std::abs( node ) , 1.0 );
 auto gain = [ & ]( double child ) {
  const double g = maximize ? node - child : child - node;
  return( std::max( g , eps ) );
  };

 std::vector< Eval > evals( cands.size() );
 for( std::size_t c = 0 ; c < cands.size() ; ++c ) {
  const auto bounds = children( cands[ c ] );
  evals[ c ] = { bounds[ 0 ] , bounds[ 1 ] ,
		 gain( bounds[ 0 ] ) * gain( bounds[ 1 ] ) };
  if( log )
   *log << "StrongBranching::evaluate: candidate " << c << " value "
	<< cands[ c ].value << " down " << bounds[ 0 ] << " up "
	<< bounds[ 1 ] << " score " << evals[ c ][ 2 ] << std::endl;
  }
 return( evals );

 }  // end( StrongBranching::evaluate )

/*--------------------------------------------------------------------------*/

std::size_t StrongBranching::best( const std::vector< Eval > & evals )
{
 std::size_t best_c = 0;
 for( std::size_t c = 1 ; c < evals.size() ; ++c )
  if( evals[ c ][ 2 ] > evals[ best_c ][ 2 ] )
   best_c = c;
 return( best_c );
 }

/*--------------------------------------------------------------------------*/

void StrongBranching::write_log( const std::string & file ,
				 unsigned long call , double node ,
				 const std::vector< BranchCandidate > & cands ,
				 const std::vector< BranchingFeatures::Features >
				 & x , const std::vector< Eval > & evals ,
				 std::size_t chosen )
{
 const bool is_new = ! std::filesystem::exists( file );
 std::ofstream out( file , std::ios::app );
 if( ! out )
  throw( std::runtime_error( "StrongBranching::write_log: cannot open " +
			     file ) );
 if( is_new )
  out << "call,node,fixed,rank,sb,pos,value,lagr,frac,cost,down,up,"
      << "score,chosen" << std::endl;
 out << std::setprecision( 10 );
 for( std::size_t c = 0 ; c < cands.size() ; ++c ) {
  const auto & cd = cands[ c ];
  out << call << "," << node << "," << x[ c ][ 6 ] << "," << c << ","
      << cd.block << "," << x[ c ][ 5 ] << "," << cd.value << ","
      << x[ c ][ 2 ] << "," << x[ c ][ 0 ] << "," << cd.cost << ","
      << evals[ c ][ 0 ] << "," << evals[ c ][ 1 ] << "," << evals[ c ][ 2 ]
      << "," << ( c == chosen ? 1 : 0 ) << std::endl;
  }

 }  // end( StrongBranching::write_log )

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF LearnedBranching ----------------------*/
/*--------------------------------------------------------------------------*/

void LearnedBranching::load( const std::string & file )
{
 static const std::string _prfx = "LearnedBranching::load: ";
 constexpr auto NF = BranchingFeatures::NFeatures;

 f_w1.clear();
 if( file.empty() )
  throw( std::invalid_argument( _prfx + "no file given" ) );
 std::ifstream in( file );
 if( ! in )
  throw( std::runtime_error( _prfx + "cannot open " + file ) );

 std::string tag;
 std::size_t nin , nhid;
 if( ! ( in >> tag >> nin >> nhid ) || ( tag != "LBRModel" ) ||
     ( nin != NF ) || ( nhid == 0 ) )
  throw( std::invalid_argument( _prfx + file + " is not a model of " +
				std::to_string( NF ) + " features" ) );

 std::vector< std::vector< double > > w1( nhid ,
					  std::vector< double >( nin ) );
 std::vector< double > b1( nhid ) , w2( nhid );
 double b2;
 for( auto & row : w1 )
  for( auto & w : row )
   in >> w;
 for( auto & b : b1 )
  in >> b;
 for( auto & w : w2 )
  in >> w;
 in >> b2;
 if( ! in )
  throw( std::invalid_argument( _prfx + file + " is truncated" ) );

 f_w1 = std::move( w1 );
 f_b1 = std::move( b1 );
 f_w2 = std::move( w2 );
 f_b2 = b2;

 }  // end( LearnedBranching::load )

/*--------------------------------------------------------------------------*/

double LearnedBranching::score( const BranchingFeatures::Features & x ) const
{
 double score = f_b2;
 for( std::size_t h = 0 ; h < f_w1.size() ; ++h ) {
  double a = f_b1[ h ];
  for( std::size_t i = 0 ; i < BranchingFeatures::NFeatures ; ++i )
   a += f_w1[ h ][ i ] * x[ i ];
  score += f_w2[ h ] * std::tanh( a );
  }
 return( score );

 }  // end( LearnedBranching::score )

/*--------------------------------------------------------------------------*/

std::size_t LearnedBranching::best(
		     const std::vector< BranchingFeatures::Features > & x ) const
{
 std::size_t best_c = 0;
 double best_score = - std::numeric_limits< double >::infinity();
 for( std::size_t c = 0 ; c < x.size() ; ++c ) {
  const double s = score( x[ c ] );
  if( s > best_score ) {
   best_score = s;
   best_c = c;
   }
  }
 return( best_c );

 }  // end( LearnedBranching::best )

/*--------------------------------------------------------------------------*/
/*---------------------- End File BranchingRule.cpp ------------------------*/
/*--------------------------------------------------------------------------*/

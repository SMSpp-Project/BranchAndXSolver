/*--------------------------------------------------------------------------*/
/*---------------------------- File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Unit test of BranchAndXSolver on Blocks of the core alone.
 *
 * The Block is an AbstractBlock with N ColVariable x, a BoxConstraint on
 * each of them and a separable quadratic DQuadFunction as Objective,
 *
 *     min ( or max ) sum_i q_i x_i^2 + b_i x_i   s.t.   l_i <= x_i <= u_i ,
 *                                                        x_i integer ,
 *
 * with q_i > 0 for a minimization and q_i < 0 for a maximization, so that
 * the integer optimum is known in closed form: in each coordinate the
 * nearest integer to the unconstrained optimum - b_i / ( 2 q_i ) inside
 * [ ceil( l_i ) , floor( u_i ) ]. The closed form is checked against the
 * enumeration of all the integer points of the box.
 *
 * The core has no RelaxationSolver, so the test defines the two inner Solver
 * that BranchAndXSolver drives, both BoxSolver underneath:
 *
 * - BoxRelaxation, a BoxSolver and a RelaxationSolver: the ColVariable are
 *   continuous in the Block and their integrality is the knowledge of this
 *   Solver, which solves the continuous relaxation with BoxSolver, reports
 *   it as a true solution when it is integer, and branches on its most
 *   fractional coordinate by tightening the BoxConstraint (a BoundChange);
 *
 * - BoxRounding, a BoxSolver and a (heuristic) ChangeSolver, which rounds
 *   the solution of the relaxation to the nearest integer inside the box;
 *
 * - BoxEstimate, the BoxRelaxation that also estimates the best value under
 *   each child of a branching [see BranchEstimator], as the value of the
 *   rounding in the box of the child, and BoxWrongEstimate, which swaps the
 *   estimates of the two children and at times gives none.
 *
 * The inner Solver are given to BranchAndXSolver by the BlockSolverConfig
 * in InnerBSCfg.txt (the relaxation alone), InnerBSCfg-rounding.txt (the
 * relaxation and the heuristic), InnerBSCfg-estimate.txt and
 * InnerBSCfg-wrong-estimate.txt (the relaxations with the estimates), read
 * from the working directory.
 *
 * The test needs nothing but the core, so that the CI of BranchAndXSolver
 * builds this module alone. The relaxation acts on the Block rather than on
 * its own state, hence the parallel depth-first exploration, whose workers
 * move their own Solver independently on the one Block, is not exercised.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "AbstractBlock.h"
#include "BoxSolver.h"
#include "BranchAndXSolver.h"
#include "BranchingRule.h"
#include "ChangeSolver.h"
#include "DQuadFunction.h"
#include "FRealObjective.h"
#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BoundChange ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the Change setting both the bounds of a BoxConstraint

class BoundChange : public Change {

 public:

 explicit BoundChange( BoxConstraint * c = nullptr , double l = 0 ,
                       double u = 0 ) : Change() , f_c( c ) , f_l( l ) ,
                                        f_u( u ) {}

 Change * apply( Block * block , bool doUndo = false ,
                 ModParam issueMod = eNoBlck ,
                 ModParam issueAMod = eNoBlck ) override {
  auto undo = doUndo ? new BoundChange( f_c , f_c->get_lhs() ,
                                        f_c->get_rhs() ) : nullptr;
  f_c->set_lhs( f_l , issueMod );
  f_c->set_rhs( f_u , issueMod );
  return( undo );
  }

 void deserialize( const netCDF::NcGroup & group ) override {
  throw( std::logic_error( "BoundChange::deserialize: not supported" ) );
  }

 private:

 BoxConstraint * f_c;  ///< the BoxConstraint
 double f_l;           ///< its new lower bound
 double f_u;           ///< its new upper bound

 SMSpp_insert_in_factory_h;

 };  // end( class( BoundChange ) )

SMSpp_insert_in_factory_cpp_0( BoundChange );

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS BoxBase ------------------------------*/
/*--------------------------------------------------------------------------*/
/// a BoxSolver that re-solves at each compute() and keeps its solution
/** The Block must have the static groups "x" of ColVariable and "box" of
 * BoxConstraint, one per ColVariable. Each compute() solves the problem
 * from scratch, since the BoundChange are applied to the Block without any
 * Modification, and keeps the values of the ColVariable it has found. */

class BoxBase : public BoxSolver {

 public:

 BoxBase( void ) : BoxSolver() { set_sol( 1 ); }

 void set_Block( Block * block ) override {
  BoxSolver::set_Block( block );
  f_x = nullptr;
  f_box = nullptr;
  if( ! block )
   return;
  f_x = block->get_static_variable_v< ColVariable >( "x" );
  f_box = block->get_static_constraint_v< BoxConstraint >( "box" );
  if( ( ! f_x ) || ( ! f_box ) || ( f_x->size() != f_box->size() ) ||
      f_x->empty() )
   throw( std::invalid_argument( "BoxBase::set_Block: unsupported Block" ) );
  }

 int compute( bool changedvars = true ) override {
  // a ConstraintMod on a BoxConstraint makes the BoxSolver forget its
  // solution, i.e., solve again
  sp_Mod mod = std::make_shared< ConstraintMod >( & f_box->front() );
  BoxSolver::add_Modification( mod );
  const int status = BoxSolver::compute( changedvars );
  ++n_compute;
  f_val.clear();
  if( BoxSolver::has_var_solution() )
   for( const auto & v : *f_x )
    f_val.push_back( v.get_value() );
  return( status );
  }

 protected:

 /// apply a BoundChange to the Block, issuing no Modification
 Change * apply_bounds( Change * chg , bool doUndo ) {
  if( ! dynamic_cast< BoundChange * >( chg ) )
   throw( std::invalid_argument( "BoxBase::apply_bounds: not a BoundChange"
                                 ) );
  return( chg->apply( f_Block , doUndo , eNoMod , eNoMod ) );
  }

 /// write the given values in the ColVariable
 void write( const std::vector< double > & val ) {
  for( std::size_t i = 0 ; i < val.size() ; ++i )
   ( *f_x )[ i ].set_value( val[ i ] );
  }

 /// the value of the Objective at the point x
 double objective_value( const std::vector< double > & x ) const {
  auto obj = static_cast< FRealObjective * >( f_Block->get_objective() );
  auto qf = static_cast< DQuadFunction * >( obj->get_function() );
  double v = qf->get_constant_term();
  for( const auto & el : qf->get_v_var() ) {
   const double r = x[ std::get< 0 >( el ) - & f_x->front() ];
   v += std::get< 2 >( el ) * r * r + std::get< 1 >( el ) * r;
   }
  return( v );
  }

 std::vector< ColVariable > * f_x = nullptr;     ///< the ColVariable
 std::vector< BoxConstraint > * f_box = nullptr;  ///< their bounds
 std::vector< double > f_val;  ///< the solution of the last compute()

 public:

 /// the compute() of all the BoxBase so far, which tells how many nodes a
 /// solve has evaluated
 static inline int n_compute = 0;

 };  // end( class( BoxBase ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS BoxRelaxation --------------------------*/
/*--------------------------------------------------------------------------*/
/// the continuous relaxation of the integer box problem, as RelaxationSolver

class BoxRelaxation : public BoxBase , public RelaxationSolver {

 public:

 Change * apply( Change * chg , bool doUndo = false ) override {
  return( apply_bounds( chg , doUndo ) );
  }

 /// branch on the most fractional coordinate, the nearest side first
 std::vector< Change * > branch( void ) override {
  std::size_t j = f_val.size();
  double dist = 0;
  for( std::size_t i = 0 ; i < f_val.size() ; ++i ) {
   const double fr = f_val[ i ] - std::floor( f_val[ i ] );
   const double d = std::min( fr , 1 - fr );
   if( d > dist + eps ) {
    dist = d;
    j = i;
    }
   }
  if( j == f_val.size() )
   throw( std::logic_error( "BoxRelaxation::branch: integer solution" ) );

  auto & c = ( *f_box )[ j ];
  const double v = f_val[ j ];
  Change * down = new BoundChange( & c , c.get_lhs() , std::floor( v ) );
  Change * up = new BoundChange( & c , std::ceil( v ) , c.get_rhs() );
  f_j = j;
  f_down_first = v - std::floor( v ) < 0.5;
  trace.push_back( double( j ) + v / 1000 );
  if( f_down_first )
   return( std::vector< Change * >{ down , up } );
  return( std::vector< Change * >{ up , down } );
  }

 bool has_true_var_solution( void ) override {
  if( f_val.empty() )
   return( false );
  for( auto v : f_val )
   if( std::abs( v - std::round( v ) ) > eps )
    return( false );
  return( true );
  }

 Solver::OFValue get_true_lb( void ) override {
  return( has_true_var_solution() ? get_var_value()
                                  : - Inf< Solver::OFValue >() );
  }

 Solver::OFValue get_true_ub( void ) override {
  return( has_true_var_solution() ? get_var_value()
                                  : Inf< Solver::OFValue >() );
  }

 void get_true_var_solution( Configuration * solc = nullptr ) override {
  write( f_val );
  }

 Solution * get_true_solution( Configuration * solc = nullptr ) override {
  get_true_var_solution( solc );
  return( f_Block->get_Solution( solc , false ) );
  }

 /// the coordinate and the value branched upon by every branch() so far,
 /// which tells the order in which the nodes are explored
 static inline std::vector< double > trace;

 protected:

 std::size_t f_j = 0;          ///< the coordinate of the last branch()
 bool f_down_first = true;     ///< whether its "down" child came first

 private:

 static constexpr double eps = 1e-9;  ///< integrality tolerance

 SMSpp_insert_in_factory_h;

 };  // end( class( BoxRelaxation ) )

SMSpp_insert_in_factory_cpp_0( BoxRelaxation );

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BoxRounding ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the rounding of the relaxation inside the box, as heuristic ChangeSolver

class BoxRounding : public BoxBase , public ChangeSolver {

 public:

 Change * apply( Change * chg , bool doUndo = false ) override {
  return( apply_bounds( chg , doUndo ) );
  }

 int compute( bool changedvars = true ) override {
  f_status = BoxBase::compute( changedvars );
  if( f_status != kOK )
   return( f_status );

  f_round.resize( f_val.size() );
  for( std::size_t i = 0 ; i < f_val.size() ; ++i ) {
   const double lo = std::ceil( ( *f_box )[ i ].get_lhs() );
   const double hi = std::floor( ( *f_box )[ i ].get_rhs() );
   if( lo > hi )           // no integer point in the box
    return( f_status = kInfeasible );
   f_round[ i ] = std::clamp( std::round( f_val[ i ] ) , lo , hi );
   }

  f_value = objective_value( f_round );
  return( f_status );
  }

 Solver::OFValue get_lb( void ) override { return( f_value ); }

 Solver::OFValue get_ub( void ) override { return( f_value ); }

 Solver::OFValue get_var_value( void ) override { return( f_value ); }

 bool has_var_solution( void ) override { return( f_status == kOK ); }

 void get_var_solution( Configuration * solc = nullptr ) override {
  write( f_round );
  }

 private:

 int f_status = kUnEval;       ///< the status of the last compute()
 std::vector< double > f_round;  ///< the rounded solution
 double f_value = 0;           ///< its value

 SMSpp_insert_in_factory_h;

 };  // end( class( BoxRounding ) )

SMSpp_insert_in_factory_cpp_0( BoxRounding );

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS BoxEstimate ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the relaxation that also estimates the best value under each child of a
/// branching [see BranchEstimator]: the value of the rounding of the
/// solution of the relaxation to the nearest integer in the box of the
/// child, NaN if that box has no integer point

class BoxEstimate : public BoxRelaxation , public BranchEstimator {

 public:

 std::vector< double > branch_estimates( void ) override {
  ++n_estimates;
  std::vector< double > e;
  for( bool down : { f_down_first , ! f_down_first } )
   e.push_back( rounded_value( down ) );
  return( e );
  }

 /// the calls of branch_estimates() of all the BoxEstimate so far
 static inline int n_estimates = 0;

 protected:

 // the value of the rounding of the solution of the relaxation in the box
 // of the child "down" (or "up") of the last branch(), NaN if none
 double rounded_value( bool down ) const {
  std::vector< double > r( f_val.size() );
  for( std::size_t i = 0 ; i < f_val.size() ; ++i ) {
   double lo = std::ceil( ( *f_box )[ i ].get_lhs() );
   double hi = std::floor( ( *f_box )[ i ].get_rhs() );
   if( i == f_j ) {
    if( down )
     hi = std::min( hi , std::floor( f_val[ i ] ) );
    else
     lo = std::max( lo , std::ceil( f_val[ i ] ) );
    }
   if( lo > hi )
    return( std::numeric_limits< double >::quiet_NaN() );
   r[ i ] = std::clamp( std::round( f_val[ i ] ) , lo , hi );
   }
  return( objective_value( r ) );
  }

 SMSpp_insert_in_factory_h;

 };  // end( class( BoxEstimate ) )

SMSpp_insert_in_factory_cpp_0( BoxEstimate );

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS BoxWrongEstimate -------------------------*/
/*--------------------------------------------------------------------------*/
/// the BoxEstimate with the estimates of the two children swapped, and none
/// at every third branching: the order of a best-estimate search is
/// misled, but not its optimum, since the nodes are still pruned by their
/// bounds

class BoxWrongEstimate : public BoxEstimate {

 public:

 std::vector< double > branch_estimates( void ) override {
  auto e = BoxEstimate::branch_estimates();
  if( ( ++f_calls % 3 ) == 0 )
   return( std::vector< double >() );
  std::swap( e[ 0 ] , e[ 1 ] );
  return( e );
  }

 private:

 int f_calls = 0;              ///< the calls so far

 SMSpp_insert_in_factory_h;

 };  // end( class( BoxWrongEstimate ) )

SMSpp_insert_in_factory_cpp_0( BoxWrongEstimate );


/*--------------------------------------------------------------------------*/
/*------------------------------- DATA -------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

const int N = 6;            // number of ColVariable
const double tol = 1e-8;    // absolute accuracy asked of the value
const double INF = Inf< double >();

int n_fail = 0;             // number of failed checks

/// the data of one ColVariable: bounds and objective coefficients
struct Coord {
 double l , u , q , b;
 };

/// an instance: the data, the sense and the Block built out of them
struct Instance {
 std::vector< Coord > data;
 bool max = false;
 AbstractBlock * block = nullptr;
 DQuadFunction * f = nullptr;
 // the groups of the Block, which does not own them
 std::unique_ptr< std::vector< ColVariable > > x;
 std::unique_ptr< std::vector< BoxConstraint > > box;
 };

/*--------------------------------------------------------------------------*/
/// record a check, printing it when it fails

void check( bool ok , const std::string & what )
{
 if( ok )
  return;
 ++n_fail;
 std::cout << "  KO: " << what << std::endl;
 }

/*--------------------------------------------------------------------------*/
/// random data, with fixed seed: fractional bounds and a fractional
/// unconstrained optimum in each coordinate

std::vector< Coord > random_data( unsigned seed , bool max )
{
 std::mt19937 rg( seed );
 std::uniform_real_distribution< double > ul( -4 , -0.5 ) , uu( 0.5 , 4 ) ,
  ua( 0.5 , 2 ) , ub( -6 , 6 );
 std::vector< Coord > data;
 for( int i = 0 ; i < N ; ++i ) {
  const double l = ul( rg ) , u = uu( rg ) , a = ua( rg ) , b = ub( rg );
  data.push_back( { l , u , max ? - a : a , b } );
  }
 return( data );
 }

/*--------------------------------------------------------------------------*/
/// the value of the integer point x

double value( const std::vector< Coord > & data ,
              const std::vector< double > & x )
{
 double v = 0;
 for( std::size_t i = 0 ; i < data.size() ; ++i )
  v += data[ i ].q * x[ i ] * x[ i ] + data[ i ].b * x[ i ];
 return( v );
 }

/*--------------------------------------------------------------------------*/
/// the integer optimum in closed form, INF ( -INF ) if there is none

double closed_form( const std::vector< Coord > & data , bool max )
{
 std::vector< double > x;
 for( const auto & d : data ) {
  const double lo = std::ceil( d.l ) , hi = std::floor( d.u );
  if( lo > hi )
   return( max ? - INF : INF );
  const double c = std::clamp( - d.b / ( 2 * d.q ) , lo , hi );
  const double f = std::floor( c ) , g = std::ceil( c );
  const double vf = d.q * f * f + d.b * f , vg = d.q * g * g + d.b * g;
  x.push_back( ( max ? vf >= vg : vf <= vg ) ? f : g );
  }
 return( value( data , x ) );
 }

/*--------------------------------------------------------------------------*/
/// the integer optimum by enumeration of all the integer points of the box

double enumeration( const std::vector< Coord > & data , bool max )
{
 double best = max ? - INF : INF;
 std::vector< double > x( data.size() );
 std::function< void( std::size_t ) > visit = [ & ]( std::size_t i ) {
  if( i == data.size() ) {
   const double v = value( data , x );
   if( max ? v > best : v < best )
    best = v;
   return;
   }
  for( double k = std::ceil( data[ i ].l ) ; k <= std::floor( data[ i ].u ) ;
       ++k ) {
   x[ i ] = k;
   visit( i + 1 );
   }
  };
 visit( 0 );
 return( best );
 }

/*--------------------------------------------------------------------------*/
/// build the Block of the instance

void build( Instance & inst )
{
 inst.block = new AbstractBlock();
 inst.x = std::make_unique< std::vector< ColVariable > >( inst.data.size() );
 auto & x = inst.x;
 inst.block->add_static_variable( *x , "x" );

 inst.box = std::make_unique< std::vector< BoxConstraint > >(
                                                          inst.data.size() );
 DQuadFunction::v_coeff_triple triples;
 for( std::size_t i = 0 ; i < inst.data.size() ; ++i ) {
  ( *inst.box )[ i ].set_variable( & ( *x )[ i ] );
  ( *inst.box )[ i ].set_lhs( inst.data[ i ].l );
  ( *inst.box )[ i ].set_rhs( inst.data[ i ].u );
  triples.emplace_back( & ( *x )[ i ] , inst.data[ i ].b ,
                        inst.data[ i ].q );
  }
 inst.block->add_static_constraint( *inst.box , "box" );

 inst.f = new DQuadFunction( std::move( triples ) );
 auto obj = new FRealObjective( inst.block , inst.f );
 obj->set_sense( inst.max ? Objective::eMax : Objective::eMin , eNoMod );
 inst.block->set_objective( obj );
 }

/// an instance with random data
Instance make( unsigned seed , bool max )
{
 Instance inst;
 inst.data = random_data( seed , max );
 inst.max = max;
 build( inst );
 return( inst );
 }

/*--------------------------------------------------------------------------*/
/// a BranchAndXSolver with the given inner BlockSolverConfig, parameters by
/// name

Solver * new_bnx( const std::string & bscfg = "InnerBSCfg.txt" )
{
 auto bnx = Solver::new_Solver( "BranchAndXSolver" );
 bnx->set_par( bnx->str_par_str2idx( "strNameOfBlockSolverConfigurationFile" )
               , bscfg );
 bnx->set_par( bnx->dbl_par_str2idx( "dblRelAcc" ) , 1e-12 );
 return( bnx );
 }

void set_int( Solver * s , const std::string & name , int value )
{
 s->set_par( s->int_par_str2idx( name ) , value );
 }

void set_str( Solver * s , const std::string & name ,
              const std::string & value )
{
 s->set_par( s->str_par_str2idx( name ) , value );
 }

int get_int( Solver * s , const std::string & name )
{
 return( s->get_int_par( s->int_par_str2idx( name ) ) );
 }

double get_dbl( Solver * s , const std::string & name )
{
 return( s->get_dbl_par( s->dbl_par_str2idx( name ) ) );
 }

std::string get_str( Solver * s , const std::string & name )
{
 return( s->get_str_par( s->str_par_str2idx( name ) ) );
 }

/*--------------------------------------------------------------------------*/
/// solve and check status, value, bounds and solution against the optimum

void solve_and_check( Solver * bnx , Instance & inst ,
                      const std::string & name )
{
 const double opt = closed_form( inst.data , inst.max );
 int status = Solver::kError;
 try {
  status = bnx->compute();
  }
 catch( std::exception & e ) {
  check( false , name + ": compute() threw " + e.what() );
  return;
  }
 const double v = bnx->get_var_value();
 const double lb = bnx->get_lb() , ub = bnx->get_ub();
 const bool ok = ( status == Solver::kOK ) && ( std::abs( v - opt ) <= tol );
 std::cout << std::left << std::setw( 40 ) << name << std::right
           << std::scientific << std::setprecision( 9 ) << " value " << v
           << " optimum " << opt << " status " << status
           << ( ok ? "  OK" : "  KO" ) << std::endl;
 check( status == Solver::kOK , name + ": status " +
        std::to_string( status ) + ", not kOK" );
 check( std::abs( v - opt ) <= tol , name + ": wrong value" );
 // the enumeration is complete, so both bounds are the optimum
 check( ( std::abs( lb - opt ) <= tol ) && ( std::abs( ub - opt ) <= tol ) ,
        name + ": lb " + std::to_string( lb ) + " ub " +
        std::to_string( ub ) + " not the optimum" );

 // the solution written in the Block is integer, in the box, optimal
 bnx->get_var_solution();
 auto x = inst.block->get_static_variable_v< ColVariable >( "x" );
 std::vector< double > xv;
 bool in_box = true;
 for( std::size_t i = 0 ; i < x->size() ; ++i ) {
  const double xi = ( *x )[ i ].get_value();
  xv.push_back( xi );
  in_box &= ( std::abs( xi - std::round( xi ) ) <= 1e-9 ) &&
            ( xi >= inst.data[ i ].l ) && ( xi <= inst.data[ i ].u );
  }
 check( in_box , name + ": the solution is not an integer point of the box" );
 check( std::abs( value( inst.data , xv ) - opt ) <= tol ,
        name + ": the solution is not optimal" );
 }

/*--------------------------------------------------------------------------*/
/// solve and check the status and that the bounds hold on both sides of the
/// optimum, as they have to after a stop

void solve_and_check_stop( Solver * bnx , Instance & inst , int expected ,
                           const std::string & name )
{
 const double opt = closed_form( inst.data , inst.max );
 int status = Solver::kError;
 try {
  status = bnx->compute();
  }
 catch( std::exception & e ) {
  check( false , name + ": compute() threw " + e.what() );
  return;
  }
 const double lb = bnx->get_lb() , ub = bnx->get_ub();
 std::cout << std::left << std::setw( 40 ) << name << std::right
           << std::scientific << std::setprecision( 9 ) << " lb " << lb
           << " ub " << ub << " status " << status
           << ( status == expected ? "  OK" : "  KO" ) << std::endl;
 check( status == expected , name + ": status " + std::to_string( status ) +
        ", not " + std::to_string( expected ) );
 check( ( lb <= opt + tol ) && ( ub >= opt - tol ) ,
        name + ": the bounds do not hold" );
 // the enumeration has not been completed: the bound on the side of the
 // nodes not yet explored is the best of their relaxation bounds, which is
 // finite once the root has been evaluated, as a budget of nodes ensures
 if( expected == Solver::kStopIter )
  check( inst.max ? ub < INF : lb > - INF ,
         name + ": no relaxation bound after the root has been evaluated" );
 }

/*--------------------------------------------------------------------------*/
/// delete the Block of the instance, the Solver registered to it and its
/// Objective, which the Block does not own

void destroy( Instance & inst )
{
 inst.block->unregister_Solvers( true );
 auto obj = inst.block->get_objective();
 delete inst.block;  // it clear()-s the Objective, which is then deleted
 delete obj;
 inst.block = nullptr;
 inst.box.reset();
 inst.x.reset();
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- THE CASES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// the closed form agrees with the enumeration, minimization and maximization

void test_closed_form( void )
{
 std::cout << "closed form against enumeration" << std::endl;
 for( unsigned seed : { 7u , 11u , 13u } )
  for( bool max : { false , true } ) {
   const auto data = random_data( seed , max );
   const double cf = closed_form( data , max );
   const double en = enumeration( data , max );
   check( std::abs( cf - en ) <= tol , "seed " + std::to_string( seed ) +
          ( max ? " max" : " min" ) + ": closed form " + std::to_string( cf )
          + " enumeration " + std::to_string( en ) );
   }
 }

/*--------------------------------------------------------------------------*/
/// every exploration strategy under both bounding protocols, minimization

void test_strategies( void )
{
 std::cout << "strategies, minimization" << std::endl;
 const char * method[] = { "DFS" , "BFS" , "BestFS" , "BestFSDive" ,
                          "BestEstimate" };
 for( int m = 0 ; m < 5 ; ++m )
  for( int p = 0 ; p < 2 ; ++p ) {
   auto inst = make( 7 , false );
   auto bnx = new_bnx();
   set_int( bnx , "intSolveMethod" , m );
   set_int( bnx , "intBoundingProtocol" , p );
   inst.block->register_Solver( bnx );
   solve_and_check( bnx , inst , std::string( method[ m ] ) +
                    ( p ? ", lazy" : ", eager" ) );
   destroy( inst );
   }
 }

/*--------------------------------------------------------------------------*/
/// every exploration strategy on a maximization

void test_maximization( void )
{
 std::cout << "strategies, maximization" << std::endl;
 const char * method[] = { "DFS" , "BFS" , "BestFS" , "BestFSDive" ,
                          "BestEstimate" };
 for( int m = 0 ; m < 5 ; ++m ) {
  auto inst = make( 11 , true );
  auto bnx = new_bnx();
  set_int( bnx , "intSolveMethod" , m );
  inst.block->register_Solver( bnx );
  solve_and_check( bnx , inst , std::string( "max, " ) + method[ m ] );
  destroy( inst );
  }
 }

/*--------------------------------------------------------------------------*/
/// the branch-and-bound of PBeFS [see intTreeExplorer]: every
/// strategy under both bounding protocols, in both senses, by every scheme
/// that strPBeFSScheme names, a node budget stopping it, and the
/// infeasible boxes, while an unknown scheme is refused; one worker only,
/// since the Box Solver write their bounds in the Block that the workers
/// share (the tests with several workers are the batteries of the
/// BinaryKnapsackBlock)

void test_tree_explorer( void )
{
 std::cout << "the branch-and-bound of PBeFS" << std::endl;
 const char * method[] = { "DFS" , "BFS" , "BestFS" , "BestFSDive" ,
                          "BestEstimate" };
 const auto label = []( const std::string & t ) {
  return( "PBeFS " + t + ", " );
  };
 const auto pbefs = []( Solver * bnx , const std::string & t ) {
  set_int( bnx , "intTreeExplorer" , BranchAndXSolver::PBeFS );
  set_str( bnx , "strPBeFSScheme" , t );
  };
 const std::string explorers[] = { "sequential" , "shared" , "distributed" ,
                                   "deterministic" , "selfsplit" };

 bool refused = false;
 auto none = new_bnx();
 try {
  set_str( none , "strPBeFSScheme" , "nonesuch" );
  }
 catch( const std::invalid_argument & ) {
  refused = true;
  }
 delete none;
 check( refused , "PBeFS: an unknown scheme is refused" );

 for( const auto & t : explorers )
  for( bool max : { false , true } )
   for( int m = 0 ; m < 5 ; ++m )
    for( int p = 0 ; p < 2 ; ++p ) {
     auto inst = make( max ? 11 : 7 , max );
     auto bnx = new_bnx();
     pbefs( bnx , t );
     set_int( bnx , "intSolveMethod" , m );
     set_int( bnx , "intBoundingProtocol" , p );
     inst.block->register_Solver( bnx );
     solve_and_check( bnx , inst , label( t ) +
                      ( max ? "max, " : "min, " ) + method[ m ] +
                      ( p ? ", lazy" : ", eager" ) );
     destroy( inst );
     }

 for( const auto & t : explorers )
  for( bool max : { false , true } ) {
   const std::string tag = max ? ", max" : ", min";
   auto inst = make( 7 , max );
   auto bnx = new_bnx();
   pbefs( bnx , t );
   set_int( bnx , "intMaxNodes" , 2 );
   inst.block->register_Solver( bnx );
   solve_and_check_stop( bnx , inst , Solver::kStopIter ,
                         label( t ) + "two nodes" + tag );
   set_int( bnx , "intMaxNodes" , 0 );
   solve_and_check( bnx , inst , label( t ) + "no node limit" + tag );
   destroy( inst );
   }

 // an empty box, and a box with no integer point whose relaxation is not
 // empty, as in test_infeasible(): an infeasible node is no infeasible root
 for( const auto & t : explorers )
  for( int kind = 0 ; kind < 2 ; ++kind )
   for( bool max : { false , true } ) {
    Instance inst;
    inst.data = random_data( 7 , max );
    inst.max = max;
    if( kind == 0 )
     inst.data[ 2 ].l = 1 , inst.data[ 2 ].u = 0;
    else
     inst.data[ 2 ].l = 0.2 , inst.data[ 2 ].u = 0.8;
    build( inst );
    auto bnx = new_bnx();
    pbefs( bnx , t );
    inst.block->register_Solver( bnx );
    const std::string name = label( t ) +
                             ( kind ? "no integer point" : "empty box" ) +
                             ( max ? ", max" : ", min" );
    const int status = bnx->compute();
    std::cout << std::left << std::setw( 40 ) << name << std::right
              << " status " << status << std::endl;
    check( status == Solver::kInfeasible , name + ": status " +
           std::to_string( status ) + ", not kInfeasible" );
    destroy( inst );
    }
 }

/*--------------------------------------------------------------------------*/
/// the best-estimate search with the estimates of the children given by
/// the relaxation [see BranchEstimator], right or misleading (swapped, and
/// at times none), by both explorers, under both bounding protocols, in
/// both senses: the estimates are asked for, the optimum is the same, and
/// the order of the branchings depends on the estimates in at least one
/// case of each explorer

void test_best_estimate( void )
{
 std::cout << "best-estimate with the estimates of the relaxation"
           << std::endl;
 std::map< std::tuple< bool , int , bool , int > ,
           std::vector< double > > order;
 for( const char * cfg : { "InnerBSCfg-estimate.txt" ,
                           "InnerBSCfg-wrong-estimate.txt" } )
  for( int t : { int( BranchAndXSolver::OwnExplorer ) ,
                 int( BranchAndXSolver::PBeFS ) } )
   for( bool max : { false , true } )
    for( int p = 0 ; p < 2 ; ++p ) {
     auto inst = make( max ? 11 : 7 , max );
     auto bnx = new_bnx( cfg );
     set_int( bnx , "intTreeExplorer" , t );
     set_int( bnx , "intSolveMethod" , BranchAndXSolver::BestEstimate );
     set_int( bnx , "intBoundingProtocol" , p );
     inst.block->register_Solver( bnx );
     const int before = BoxEstimate::n_estimates;
     const bool wrong = std::string( cfg ).find( "wrong" ) !=
                        std::string::npos;
     const std::string name = std::string( wrong ? "wrong" : "right" ) +
                              " estimates, " +
                              ( t ? "PBeFS, " : "own, " ) +
                              ( max ? "max, " : "min, " ) +
                              ( p ? "lazy" : "eager" );
     BoxRelaxation::trace.clear();
     solve_and_check( bnx , inst , name );
     check( BoxEstimate::n_estimates > before ,
            name + ": no estimate asked for" );
     order[ { wrong , t , max , p } ] = BoxRelaxation::trace;
     destroy( inst );
     }

 for( int t : { int( BranchAndXSolver::OwnExplorer ) ,
                int( BranchAndXSolver::PBeFS ) } ) {
  bool differ = false;
  for( bool max : { false , true } )
   for( int p = 0 ; p < 2 ; ++p )
    differ |= order[ { false , t , max , p } ] !=
              order[ { true , t , max , p } ];
  check( differ , std::string( t ? "PBeFS" : "own" ) +
         ": the order does not depend on the estimates" );
  }
 }

/*--------------------------------------------------------------------------*/
/// the pseudocosts and their priors read from a file [see BranchingRule.h],
/// against values computed by hand, with their limit cases: no observation
/// at all, a direction that is neither 0 nor 1, an observation that is not
/// finite or that moves nothing, a negative output of the model, and a file
/// with one model only or truncated

void test_pseudocosts( void )
{
 std::cout << "pseudocosts and their learned priors" << std::endl;
 const double nan = std::numeric_limits< double >::quiet_NaN();
 PseudoCosts pc;
 check( pc.get( 3 , 0 , nan ) == 0 , "no observation: 0" );
 check( pc.get( 3 , 1 , 2.5 ) == 2.5 , "no observation: the prior" );

 pc.record( 3 , 0 , 0.5 , 2 );        // 4 per unit
 pc.record( 3 , 0 , 0.25 , 2 );       // 8 per unit
 pc.record( 7 , 1 , 0.5 , 1 );        // 2 per unit
 pc.record( 7 , 1 , 0.5 , -1 );       // worse < 0 counts as 0
 pc.record( 9 , 0 , 0 , 5 );          // moves nothing: ignored
 pc.record( 9 , 0 , 0.5 , INF );      // infeasible child: ignored
 check( ( pc.get( 3 , 0 , nan ) == 6 ) && pc.known( 3 , 0 ) &&
        ! pc.known( 3 , 1 ) , "mean of the observations" );
 check( pc.get( 7 , 1 , nan ) == 1 , "a negative worsening counts as 0" );
 check( ( ! pc.known( 9 , 0 ) ) && ( pc.observations() == 4 ) ,
        "observations that move nothing or are infinite are ignored" );
 check( pc.get( 9 , 0 , nan ) == 6 , "no observation: the mean of all" );
 check( pc.get( 9 , 1 , 0.5 ) == 0.5 , "no observation: the prior first" );
 bool thrown = false;
 try { pc.record( 3 , 2 , 1 , 1 ); } catch( std::invalid_argument & ) {
  thrown = true; }
 check( thrown , "a direction neither 0 nor 1 is refused" );

 // two candidates, 3 at 0.25 and 7 at 0.5: branching on 3 at a node of
 // bound 10, its children worsen it by 6 * 0.25 down and by the mean "up"
 // 1 * 0.75 up, plus the cheapest rounding of 7, min( 6 * 0.5 , 1 * 0.5 )
 std::vector< BranchCandidate > cands = { { 3 , 0.25 , 0.25 , 1 , 0 , 0 } ,
                                          { 7 , 0.5 , 0.5 , 1 , 0 , 0 } };
 const auto e = pc.estimates( cands , 0 , 10 , false );
 check( ( std::abs( e[ 0 ] - ( 10 + 1.5 + 0.5 ) ) < 1e-12 ) &&
        ( std::abs( e[ 1 ] - ( 10 + 0.75 + 0.5 ) ) < 1e-12 ) ,
        "estimates of a minimization" );
 const auto m = pc.estimates( cands , 0 , 10 , true );
 check( ( std::abs( m[ 0 ] - ( 10 - 2 ) ) < 1e-12 ) &&
        ( std::abs( m[ 1 ] - ( 10 - 1.25 ) ) < 1e-12 ) ,
        "estimates of a maximization" );
 const auto p = pc.estimates( cands , 0 , 10 , false ,
                              []( std::size_t ) {
                               return( std::array< double , 2 >{ 100 , 100 }
                                       ); } );
 check( std::abs( p[ 0 ] - e[ 0 ] ) < 1e-12 ,
        "a prior is not taken where there are observations" );
 thrown = false;
 try { ( void ) pc.estimates( cands , 2 , 10 , false ); }
 catch( std::invalid_argument & ) { thrown = true; }
 check( thrown , "no such candidate is refused" );

 // two models whose output is the bias: 0.5 down, -1 (hence 0) up
 const std::string file = "pseudocosts-test.txt";
 const auto model = []( double bias ) {
  std::string s = "LBRModel 8 1\n";
  for( int i = 0 ; i < 8 ; ++i )
   s += "0 ";
  return( s + "\n0\n0\n" + std::to_string( bias ) + "\n" );
  };
 std::ofstream( file ) << model( 0.5 ) << model( -1 );
 LearnedPseudoCosts lp;
 lp.load( file );
 const auto pr = lp.prior( BranchingFeatures::Features{} , -20 );
 check( ( std::abs( pr[ 0 ] - 0.5 * 20 ) < 1e-12 ) && ( pr[ 1 ] == 0 ) ,
        "learned priors, scaled, a negative output counting as 0" );
 std::ofstream( file ) << model( 0.5 );
 thrown = false;
 try { LearnedPseudoCosts().load( file ); }
 catch( std::invalid_argument & ) { thrown = true; }
 check( thrown , "a file with one model only is refused" );
 std::ofstream( file ) << "LBRModel 8 2\n1 2 3\n";
 thrown = false;
 try { LearnedPseudoCosts().load( file ); }
 catch( std::invalid_argument & ) { thrown = true; }
 check( thrown , "a truncated file is refused" );
 std::remove( file.c_str() );
 }

/*--------------------------------------------------------------------------*/
/// the rounding heuristic gives the incumbent, whose solution is kept

void test_heuristic( void )
{
 std::cout << "relaxation and rounding heuristic" << std::endl;
 for( bool max : { false , true } ) {
  auto inst = make( 13 , max );
  auto bnx = new_bnx( "InnerBSCfg-rounding.txt" );
  inst.block->register_Solver( bnx );
  solve_and_check( bnx , inst , max ? "rounding, max" : "rounding, min" );
  destroy( inst );
  }
 }

/*--------------------------------------------------------------------------*/
/// an empty box, and a box with no integer point whose relaxation is not
/// empty: the problem is infeasible, and the bounds say so

void test_infeasible( void )
{
 std::cout << "infeasible boxes" << std::endl;
 for( int kind = 0 ; kind < 2 ; ++kind )
  for( bool max : { false , true } ) {
   Instance inst;
   inst.data = random_data( 7 , max );
   inst.max = max;
   if( kind == 0 )
    inst.data[ 2 ].l = 1 , inst.data[ 2 ].u = 0;
   else
    inst.data[ 2 ].l = 0.2 , inst.data[ 2 ].u = 0.8;
   build( inst );
   auto bnx = new_bnx();
   inst.block->register_Solver( bnx );
   const std::string name = std::string( kind ? "no integer point" :
                                         "empty box" ) + ( max ? ", max" :
                                                           ", min" );
   int status = Solver::kError;
   try {
    status = bnx->compute();
    }
   catch( std::exception & e ) {
    check( false , name + ": compute() threw " + e.what() );
    }
   const double worst = max ? - INF : INF;
   std::cout << std::left << std::setw( 40 ) << name << std::right
             << " status " << status << " lb " << bnx->get_lb() << " ub "
             << bnx->get_ub() << std::endl;
   check( status == Solver::kInfeasible , name + ": status " +
          std::to_string( status ) + ", not kInfeasible" );
   check( bnx->get_var_value() == worst , name + ": finite value" );
   // infeasibility is proven, so both bounds are the "optimum"
   check( ( bnx->get_lb() == worst ) && ( bnx->get_ub() == worst ) ,
          name + ": the bounds are not both the infinity of the sense" );
   destroy( inst );
   }
 }

/*--------------------------------------------------------------------------*/
/// a root whose relaxation is already integer: solved at the root, which
/// a budget of a single node is enough for

void test_integral_root( void )
{
 std::cout << "integer root" << std::endl;
 for( int nodes : { 0 , 1 } ) {
  Instance inst;
  inst.data = random_data( 7 , false );
  // the unconstrained optimum of each coordinate is an integer in the box
  for( auto & d : inst.data )
   d.b = - 2 * d.q * std::ceil( d.l );
  build( inst );
  auto bnx = new_bnx();
  set_int( bnx , "intMaxNodes" , nodes );
  inst.block->register_Solver( bnx );
  solve_and_check( bnx , inst , std::string( "integer root, " ) +
                   ( nodes ? "one node" : "no node limit" ) );
  destroy( inst );
  }
 }

/*--------------------------------------------------------------------------*/
/// a node budget and a time budget that stop the enumeration

void test_limits( void )
{
 std::cout << "node and time limits" << std::endl;
 // the next compute() goes on, nothing having changed, once the budget is
 // lifted: from scratch, or from the retained tree [see intReoptimize]
 for( int reopt : { 0 , 1 } )
  for( bool max : { false , true } ) {
   const std::string tag = std::string( max ? ", max" : ", min" ) +
                           ( reopt ? ", retained tree" : "" );
   auto inst = make( 7 , max );
   auto bnx = new_bnx();
   set_int( bnx , "intMaxNodes" , 2 );
   set_int( bnx , "intReoptimize" , reopt );
   inst.block->register_Solver( bnx );
   solve_and_check_stop( bnx , inst , Solver::kStopIter , "two nodes" + tag );
   set_int( bnx , "intMaxNodes" , 0 );
   solve_and_check( bnx , inst , "no node limit" + tag );
   destroy( inst );
   }

 for( bool max : { false , true } ) {
  const std::string tag = max ? ", max" : ", min";
  auto inst = make( 7 , max );
  auto bnx = new_bnx();
  const auto t = bnx->dbl_par_str2idx( "dblMaxTime" );
  bnx->set_par( t , 0.0 );
  inst.block->register_Solver( bnx );
  solve_and_check_stop( bnx , inst , Solver::kStopTime , "no time" + tag );
  bnx->set_par( t , INF );
  solve_and_check( bnx , inst , "no time limit" + tag );
  destroy( inst );
  }
 }

/*--------------------------------------------------------------------------*/
/// the parameters: names, defaults, ComputeConfig, invalid values

void test_parameters( void )
{
 std::cout << "parameters" << std::endl;
 auto bnx = new_bnx();

 // every parameter of BranchAndXSolver goes from index to name and back
 const int n_int = bnx->get_num_int_par();
 for( int i = BranchAndXSolver::intSolveMethod ; i < n_int ; ++i )
  check( bnx->int_par_str2idx( bnx->int_par_idx2str( i ) ) == i ,
         "int parameter " + std::to_string( i ) + " name round trip" );
 const int s = BranchAndXSolver::strNameOfBlockSolverConfigurationFile;
 check( bnx->str_par_str2idx( bnx->str_par_idx2str( s ) ) == s ,
        "string parameter name round trip" );
 check( n_int == BranchAndXSolver::intLastBXSPar ,
        "number of int parameters" );

 // the defaults
 check( bnx->get_dflt_int_par( BranchAndXSolver::intSolveMethod ) ==
        BranchAndXSolver::BestFS , "default intSolveMethod" );
 check( bnx->get_dflt_int_par( BranchAndXSolver::intMaxNodes ) == INT_MAX ,
        "default intMaxNodes" );
 check( bnx->get_dflt_int_par( BranchAndXSolver::intTreeExplorer ) ==
        BranchAndXSolver::OwnExplorer , "default intTreeExplorer" );

 // a ComputeConfig sets them by name, the double ones included
 ComputeConfig cc;
 cc.int_pars.push_back( { "intSolveMethod" , BranchAndXSolver::DFS } );
 cc.int_pars.push_back( { "intMaxNodes" , 7 } );
 cc.int_pars.push_back( { "intBoundingProtocol" , BranchAndXSolver::Lazy } );
 cc.dbl_pars.push_back( { "dblMaxTime" , 12.5 } );
 cc.dbl_pars.push_back( { "dblRelAcc" , 1e-9 } );
 cc.dbl_pars.push_back( { "dblAbsAcc" , 1e-7 } );
 cc.str_pars.push_back( { "strNameOfBlockSolverConfigurationFile" ,
                          "InnerBSCfg-rounding.txt" } );
 try {
  bnx->set_ComputeConfig( & cc );
  }
 catch( std::exception & e ) {
  check( false , std::string( "set_ComputeConfig threw " ) + e.what() );
  }
 check( get_int( bnx , "intSolveMethod" ) == BranchAndXSolver::DFS ,
        "intSolveMethod by name" );
 check( get_int( bnx , "intMaxNodes" ) == 7 , "intMaxNodes by name" );
 check( get_int( bnx , "intBoundingProtocol" ) == BranchAndXSolver::Lazy ,
        "intBoundingProtocol by name" );
 check( get_dbl( bnx , "dblMaxTime" ) == 12.5 , "dblMaxTime by name" );
 check( get_dbl( bnx , "dblRelAcc" ) == 1e-9 , "dblRelAcc by name" );
 check( get_dbl( bnx , "dblAbsAcc" ) == 1e-7 , "dblAbsAcc by name" );
 check( get_str( bnx , "strNameOfBlockSolverConfigurationFile" ) ==
        "InnerBSCfg-rounding.txt" , "the BlockSolverConfig file by name" );

 // intMaxNodes <= 0 means no limit
 set_int( bnx , "intMaxNodes" , -3 );
 check( get_int( bnx , "intMaxNodes" ) == INT_MAX ,
        "intMaxNodes <= 0 is no limit" );

 // invalid values are refused when they are set, not at compute() time
 const auto refused = [ & ]( const std::string & name , int value ) {
  try {
   set_int( bnx , name , value );
   }
  catch( std::invalid_argument & ) {
   return( true );
   }
  return( false );
  };
 check( refused( "intThreadForDifferentSolvers" , 0 ) ,
        "intThreadForDifferentSolvers = 0 accepted" );
 check( refused( "intSolveMethod" , 9 ) , "intSolveMethod = 9 accepted" );
 check( refused( "intBoundingProtocol" , -1 ) ,
        "intBoundingProtocol = -1 accepted" );

 // a missing BlockSolverConfig file is refused
 bool thrown = false;
 try {
  bnx->set_par( s , "no-such-file.txt" );
  }
 catch( std::invalid_argument & ) {
  thrown = true;
  }
 check( thrown , "a missing BlockSolverConfig file accepted" );

 // the factory reset brings every parameter back to its default, the name
 // of the BlockSolverConfig file (empty) included
 try {
  bnx->set_ComputeConfig( nullptr );
  }
 catch( std::exception & e ) {
  check( false , std::string( "the factory reset threw " ) + e.what() );
  }
 check( get_int( bnx , "intSolveMethod" ) == BranchAndXSolver::BestFS ,
        "intSolveMethod after the reset" );
 check( get_dbl( bnx , "dblMaxTime" ) == INF , "dblMaxTime after the reset"
        );
 check( bnx->get_str_par( s ).empty() , "the file name after the reset" );

 // a new BlockSolverConfig on an attached Solver replaces the inner Solver;
 // the name is given as a temporary, i.e., to the rvalue version
 auto inst = make( 7 , false );
 inst.block->register_Solver( bnx );
 bnx->set_par( s , "InnerBSCfg-rounding.txt" );
 bnx->set_par( s , "InnerBSCfg.txt" );
 check( get_str( bnx , "strNameOfBlockSolverConfigurationFile" ) ==
        "InnerBSCfg.txt" , "the file name given as a temporary" );
 bnx->set_par( bnx->dbl_par_str2idx( "dblRelAcc" ) , 1e-12 );
 solve_and_check( bnx , inst , "BlockSolverConfig replaced" );
 destroy( inst );
 }

/*--------------------------------------------------------------------------*/
/// the Solver is detached and attached to another Block, and back

void test_detach_reattach( void )
{
 std::cout << "detach and reattach" << std::endl;
 auto a = make( 7 , false );
 auto b = make( 11 , true );
 auto bnx = new_bnx();
 try {
  a.block->register_Solver( bnx );
  solve_and_check( bnx , a , "on the first Block" );
  a.block->unregister_Solver( bnx );
  b.block->register_Solver( bnx );
  solve_and_check( bnx , b , "on the second Block" );
  b.block->unregister_Solver( bnx );
  a.block->register_Solver( bnx );
  solve_and_check( bnx , a , "back on the first Block" );
  }
 catch( std::exception & e ) {
  check( false , std::string( "detach / reattach threw " ) + e.what() );
  }
 a.block->unregister_Solver( bnx );
 b.block->unregister_Solver( bnx );
 delete bnx;
 destroy( a );
 destroy( b );
 }

/*--------------------------------------------------------------------------*/
/// a Modification of the Objective between two compute(), from scratch and
/// reoptimizing the retained tree; a compute() with nothing changed

void test_objective_change( void )
{
 std::cout << "change of the Objective" << std::endl;
 for( int reopt : { 0 , 1 } ) {
  const std::string tag = reopt ? ", reoptimizing" : ", from scratch";
  auto inst = make( 7 , false );
  auto bnx = new_bnx();
  set_int( bnx , "intReoptimize" , reopt );
  inst.block->register_Solver( bnx );
  solve_and_check( bnx , inst , "before the change" + tag );

  // nothing changed: the same answer
  solve_and_check( bnx , inst , "nothing changed" + tag );

  // the linear coefficients of two coordinates change (not their sign
  // alone, which in a box around 0 may leave the optimal value as it is)
  for( int i : { 0 , 3 } ) {
   inst.data[ i ].b += 5;
   inst.f->modify_term( i , inst.data[ i ].b , inst.data[ i ].q );
   }
  solve_and_check( bnx , inst , "after the change" + tag );

  // and a second change, a quadratic coefficient this time
  inst.data[ 1 ].q *= 3;
  inst.f->modify_term( 1 , inst.data[ 1 ].b , inst.data[ 1 ].q );
  solve_and_check( bnx , inst , "after a second change" + tag );
  destroy( inst );
  }
 }

/*--------------------------------------------------------------------------*/
/// a Variable that becomes integer, a restriction which the retained tree
/// survives, and one that becomes continuous, which it does not, its
/// leaves no longer covering the region; the branching treats every
/// coordinate as integer, so the optimum is the same, and what tells the
/// two cases apart is the number of nodes evaluated by the re-solve, fewer
/// than from scratch on the same problem when the tree is reused

void test_integrality_change( void )
{
 std::cout << "change of the integrality" << std::endl;
 // the nodes evaluated by a re-solve after the given change, from scratch
 // or reoptimizing the retained tree
 const auto resolve = [ & ]( int reopt , bool integer ,
                             const std::string & tag ) {
  auto inst = make( 7 , false );
  if( ! integer )         // to be made continuous, it is integer first
   ( *inst.x )[ 2 ].set_type( ColVariable::kInteger , eNoMod );
  auto bnx = new_bnx();
  set_int( bnx , "intReoptimize" , reopt );
  inst.block->register_Solver( bnx );
  solve_and_check( bnx , inst , "before the change" + tag );
  ( *inst.x )[ 2 ].set_type( integer ? ColVariable::kInteger
                                     : ColVariable::kContinuous );
  const int before = BoxBase::n_compute;
  solve_and_check( bnx , inst , "after the change" + tag );
  const int nodes = BoxBase::n_compute - before;
  destroy( inst );
  return( nodes );
  };

 const int scratch_i = resolve( 0 , true , ", integer, from scratch" );
 const int reopt_i = resolve( 1 , true , ", integer, reoptimizing" );
 const int scratch_c = resolve( 0 , false , ", continuous, from scratch" );
 const int reopt_c = resolve( 1 , false , ", continuous, reoptimizing" );
 std::cout << "  nodes after the change: integer " << scratch_i << " / "
           << reopt_i << ", continuous " << scratch_c << " / " << reopt_c
           << " (from scratch / reoptimizing)" << std::endl;
 check( reopt_i < scratch_i ,
        "a Variable becoming integer has discarded the retained tree" );
 check( reopt_c == scratch_c ,
        "a Variable becoming continuous has kept the retained tree" );
 }

}  // namespace

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::vector< void ( * )( void ) > cases = {
  test_closed_form , test_strategies , test_maximization ,
  test_tree_explorer , test_best_estimate , test_pseudocosts ,
  test_heuristic ,
  test_infeasible , test_integral_root , test_limits , test_parameters ,
  test_detach_reattach , test_objective_change , test_integrality_change };

 for( auto c : cases ) {
  try {
   c();
   }
  catch( std::exception & e ) {
   check( false , std::string( "exception " ) + e.what() );
   }
  }

 std::cout << ( n_fail ? "Some test FAILED!!" : "All tests passed!!" )
           << std::endl;
 return( n_fail ? 1 : 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/

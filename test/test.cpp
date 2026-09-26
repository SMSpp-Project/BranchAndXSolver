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
 *   the solution of the relaxation to the nearest integer inside the box.
 *
 * The inner Solver are given to BranchAndXSolver by the BlockSolverConfig
 * in InnerBSCfg.txt (the relaxation alone) and InnerBSCfg-rounding.txt
 * (the relaxation and the heuristic), read from the working directory.
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
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "AbstractBlock.h"
#include "BoxSolver.h"
#include "BranchAndXSolver.h"
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

 std::vector< ColVariable > * f_x = nullptr;     ///< the ColVariable
 std::vector< BoxConstraint > * f_box = nullptr;  ///< their bounds
 std::vector< double > f_val;  ///< the solution of the last compute()

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
  if( v - std::floor( v ) < 0.5 )
   return( std::vector< Change * >{ down , up } );
  return( std::vector< Change * >{ up , down } );
  }

 /// a change of the Objective leaves the feasible region as it is
 int classify( const sp_Mod & mod ) override {
  if( auto fm = std::dynamic_pointer_cast< const FunctionMod >( mod ) )
   if( dynamic_cast< Objective * >( fm->function()->get_Observer() ) )
    return( eModObjective );
  return( eModEverything );
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

  auto obj = static_cast< FRealObjective * >( f_Block->get_objective() );
  auto qf = static_cast< DQuadFunction * >( obj->get_function() );
  f_value = qf->get_constant_term();
  for( const auto & el : qf->get_v_var() ) {
   const double r = f_round[ std::get< 0 >( el ) - & f_x->front() ];
   f_value += std::get< 2 >( el ) * r * r + std::get< 1 >( el ) * r;
   }
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
 // the enumeration has not been completed: no bound on the side of the
 // nodes not yet explored
 check( inst.max ? ub == INF : lb == - INF ,
        name + ": a dual bound is claimed after a stop" );
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
 const char * method[] = { "DFS" , "BFS" , "BestFS" , "BestFSDive" };
 for( int m = 0 ; m < 4 ; ++m )
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
 const char * method[] = { "DFS" , "BFS" , "BestFS" , "BestFSDive" };
 for( int m = 0 ; m < 4 ; ++m ) {
  auto inst = make( 11 , true );
  auto bnx = new_bnx();
  set_int( bnx , "intSolveMethod" , m );
  inst.block->register_Solver( bnx );
  solve_and_check( bnx , inst , std::string( "max, " ) + method[ m ] );
  destroy( inst );
  }
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

}  // namespace

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::vector< void ( * )( void ) > cases = {
  test_closed_form , test_strategies , test_maximization , test_heuristic ,
  test_infeasible , test_integral_root , test_limits , test_parameters ,
  test_detach_reattach , test_objective_change };

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

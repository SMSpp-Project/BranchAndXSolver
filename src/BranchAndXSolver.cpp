/*--------------------------------------------------------------------------*/
/*--------------------- File BranchAndXSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class BranchAndXSolver, a generic
 * Branch-and-Bound Solver built on top of the ChangeSolver / RelaxationSolver
 * concepts [see ChangeSolver.h]: the attached RelaxationSolver(s) provide
 * relaxation bounds, true solutions and the branching Changes, the attached
 * heuristic ChangeSolver(s) further primal bounds, and the enumeration tree is
 * navigated by applying (undo) Changes to the Solver. See the file-level
 * comment of BranchAndXSolver.h for an overview.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Filippo Magi \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Filippo Magi, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <chrono>
#include <iostream>
#include <cmath>
#include <deque>
#include <memory>
#include <stack>
#include <thread>
#include <functional>

#include "BranchAndXSolver.h"

#include "Change.h"

#include "Objective.h"

#include "PBeFS/bnb/BranchAndBound.h"

#include "PBeFS/bnb/DeterministicBranchAndBound.h"

#include "PBeFS/bnb/ParallelBranchAndBound.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BranchAndXSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0( BranchAndXSolver );

/*--------------------------------------------------------------------------*/
/*-------------------------- INTERNAL FUNCTIONS ----------------------------*/
/*--------------------------------------------------------------------------*/

/// tells if a relaxation bound cannot improve the incumbent within tolerance
/** True if \p bound cannot improve \p best by more than
 * max( absAcc , relAcc * max( | best | , 1 ) ), i.e., the node can be
 * pruned within the required optimality tolerances. */

static bool cannot_improve( double bound , double best , bool minimizing ,
                            double relAcc , double absAcc )
{
 if( std::isinf( best ) )      // no incumbent yet: everything can improve
  return( false );
 // an absAcc at its default +Inf means "not active" [see Solver::dblAbsAcc]
 const double eps = std::max( absAcc == Inf< double >() ? 0.0 : absAcc ,
                              relAcc * std::max( std::abs( best ) , 1.0 ) );
 return( minimizing ? bound >= best - eps : bound <= best + eps );
 }

/// whether \p value is strictly better than the incumbent \p best
static bool improves( double value , double best , bool minimizing )
{
 return( minimizing ? value < best : value > best );
 }

/// the Solution corresponding to the current solution of \p slvr
/** Writes the current solution of \p slvr into its Block (under lock) and
 * returns the corresponding newly minted Solution object, whose ownership
 * is transferred to the caller. This is how a Solution is obtained from a
 * heuristic :ChangeSolver, which unlike a :RelaxationSolver [see
 * RelaxationSolver::get_true_solution()] has no way of its own to produce
 * one. */

static Solution * solution_of( Solver * slvr ,
                               Configuration * solc = nullptr )
{
 auto blck = slvr->get_Block();
 // TODO: check that lock()-ing / unlock()-ing here is appropriate
 blck->lock( slvr );
 slvr->get_var_solution( solc );
 auto sol = blck->get_Solution( solc , false );  // loaded, not empty
 blck->unlock( slvr );
 return( sol );
 }

/*--------------------------------------------------------------------------*/
/// write on the log one line per explored node, in CSV form
/** When a log stream is set [see Solver::set_log()], each node is described
 * by one line as soon as its fate is decided (right after its evaluation,
 * whenever that happens [see BranchAndXSolver::intBoundingProtocol]), so that
 * the whole exploration can be reconstructed offline: which nodes were
 * generated, in which order they were evaluated, what their relaxation bound
 * was, how the incumbent moved, and why each node was fenced. \p evalTime is
 * the time spent evaluating this very node, \p elapsed the time since the
 * beginning of the solve. */

static void logNode( std::ostream * log , int iter , ExploringNode * node ,
                     double bestBound , double elapsed , double evalTime ,
                     bool boundPruned )
{
 if( ! log )
  return;
 auto father = node->get_parent();
 ( *log ) << iter << ',' << node->get_name() << ','
	  << ( father ? father->get_name() : -1 ) << ','
	  << node->get_level() << ',' << bestBound << ','
	  << node->get_bound() << ',' << elapsed << ',' << evalTime
	  << ',' << boundPruned << ',' << node->is_infeasible() << std::endl;
 }

/*--------------------------------------------------------------------------*/
/// the header of the per-node CSV log [see logNode()]

static void logHeader( std::ostream * log )
{
 if( log )
  ( *log ) << "iter,node,father,level,incumbent,dualBound,"
	      "elapsedTime,evaluationTime,boundPruned,infeasible" << std::endl;
 }

/*--------------------------------------------------------------------------*/
/// initialize the common per-solve variables
/** @param solvers list to be filled with all the :ChangeSolver to drive
 *  @param f_RelaxationSolvers the relaxation solvers to insert in solvers
 *  @param f_HeuristicSolvers the heuristic solvers to insert in solvers
 *  @param minimizing set to true if the problem is a minimization one */

static void initializeVariables(
                  std::list< ChangeSolver * > * solvers ,
                  std::vector< std::pair< RelaxationSolver * ,
                               Solver * > > * f_RelaxationSolvers ,
                  std::vector< std::pair< ChangeSolver * ,
                               Solver * > > * f_HeuristicSolvers ,
                  bool & minimizing )
{
 if( ! f_HeuristicSolvers->empty() )
  minimizing = f_HeuristicSolvers->front().second->get_Block()
                                  ->get_objective_sense() == Objective::eMin;
 else if( ! f_RelaxationSolvers->empty() )
  minimizing = f_RelaxationSolvers->front().second->get_Block()
                                  ->get_objective_sense() == Objective::eMin;
 else
  throw( std::logic_error( "BranchAndXSolver::initializeVariables: "
         "both the HeuristicSolvers and the RelaxationSolvers are empty" ) );

 // the solver set moved along the tree is made of the traits (the moves
 // only ever apply() Changes)
 for( auto & p : * f_RelaxationSolvers )
  solvers->push_back( p.first );
 for( auto & p : * f_HeuristicSolvers )
  solvers->push_back( p.first );
 }

/*--------------------------------------------------------------------------*/
/// move the given :ChangeSolver down to the son node currentNode
/** The first Solver produces (and stores in the node) the undo Change, the
 * others just apply the node Change. */

static void moveSolverToSon( Node * currentNode ,
                             std::list< ChangeSolver * > * solvers )
{
 if( ! currentNode->get_f_change() )
  return;
 bool modified = false;
 for( const auto s : *solvers ) {
  if( ! modified )
   currentNode->set_toFather( s->apply( currentNode->get_f_change() ,
                                        true ) );
  else
   s->apply( currentNode->get_f_change() , false );
  modified = true;
  }
 }

/*--------------------------------------------------------------------------*/
/// move the given :ChangeSolver back to the father of currentNode

static void moveSolverToFather( Node * currentNode ,
                                std::list< ChangeSolver * > * solvers )
{
 if( currentNode->get_toFather() )
  for( const auto s : *solvers )
   s->apply( currentNode->get_toFather() , false );
 }

/*--------------------------------------------------------------------------*/
// runs every Solver's compute() concurrently [defined below]; used by the
// parallel evaluation path of computeRelaxations() / computeHeuristic()
template< typename SolverPtr >
static bool computeAllParallel( const std::vector< SolverPtr > & slvrs ,
                                std::vector< int > & zs , int maxThreads );

/*--------------------------------------------------------------------------*/
/// fence the current node whose relaxation stopped at the cutoff
/** The relaxation \p s has returned kCutOff, i.e., it has proven its bound
 * beyond the cutoff it was given, the incumbent [see computeRelaxations()]:
 * the node is fenced by bound, with the bound the relaxation reports, if
 * better than the one it has. */

static void fence_by_cutoff( Node * currentNode , Solver * s ,
                             const bool minimizing , bool & toPrune )
{
 const double bound = minimizing ? s->get_lb() : s->get_ub();
 if( minimizing ? bound > currentNode->get_bound()
                : bound < currentNode->get_bound() )
  currentNode->set_bound( bound );
 toPrune = true;
 }

/*--------------------------------------------------------------------------*/
/// compute the relaxations at the current node
/** Computes every RelaxationSolver at the current node, updating the node
 * relaxation bound, the branching solver, and, when a true solution improves
 * it, the incumbent bestBound / bestSol; sets \p toPrune when the node can be
 * discarded. The incumbent is given to every relaxation as its cutoff
 * (dblUpCutOff when minimizing, dblLwCutOff when maximizing), so that one
 * proven beyond it stops there: returning kCutOff, and the node is fenced by
 * bound [see fence_by_cutoff()], or kInfeasible, and the node is fenced as
 * infeasible. With
 * \p nThreads > 1 (and no shared incumbent lock) the several relaxations are
 * computed concurrently, bit-identically to the serial reduction [see
 * computeAllParallel()].
 *  @return the sol_type [see Solver.h] of the computation */

static int computeRelaxations(
                  std::vector< std::pair< RelaxationSolver * ,
                               Solver * > > * f_RelaxationSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol ,
                  std::atomic< double > * incumbentCell ,
                  bool & toPrune , RelaxationSolver * & branchSolver ,
                  double relAcc = 0 , double absAcc = 0 ,
                  int nThreads = 1 , std::mutex * incumbentMutex = nullptr )
{
 // the incumbent is the cutoff of the relaxations [see dblUpCutOff in
 // Solver.h]: a relaxation proven beyond it says that the node cannot
 // improve on it, which for the search is as good as infeasible (and then,
 // at the root, it proves the incumbent optimal rather than the problem
 // infeasible)
 const double incumbent = incumbentCell->load();
 const bool has_incumbent = minimizing ? incumbent < Inf< double >()
                                       : incumbent > - Inf< double >();
 for( auto & pr : *f_RelaxationSolvers )
  pr.second->set_par( minimizing ? Solver::dblUpCutOff : Solver::dblLwCutOff ,
                      incumbent );

 // parallel evaluation of the (several) relaxations of this node: only when more
 // than one thread is asked, there are at least two solvers, and no shared
 // incumbent lock is in force (a parallel-tree worker evaluates its per-node
 // solvers serially); the reduction is bit-identical to the serial path below
 if( ( nThreads > 1 ) && ( f_RelaxationSolvers->size() > 1 ) &&
     ( ! incumbentMutex ) ) {
  const std::size_t n = f_RelaxationSolvers->size();
  std::vector< int > zs( n , INT_MIN );   // INT_MIN: not computed (early stop)
  if( computeAllParallel( *f_RelaxationSolvers , zs , nThreads ) ) {
   if( ( ! currentNode->get_toFather() ) && ( ! has_incumbent ) )
    return( Solver::kInfeasible );       // the root is infeasible
   toPrune = true;
   currentNode->set_infeasible( true );
   // an infeasible node has the worst conceivable bound: this is what tells
   // it apart, in the log and in the fenced frontier, from one fenced by
   // bound, which carries the finite bound that fenced it
   currentNode->set_bound( minimizing ? Inf< double >()
                                           : - Inf< double >() );
   return( Solver::kOK );
   }
  for( std::size_t i = 0 ; i < n ; ++i ) {
   auto & [ rs , s ] = (*f_RelaxationSolvers)[ i ];
   if( zs[ i ] == Solver::kCutOff ) {  // beyond the incumbent [see above]
    fence_by_cutoff( currentNode , s , minimizing , toPrune );
    return( Solver::kOK );
    }
   if( zs[ i ] != ThinComputeInterface::kOK )
    return( zs[ i ] );
   if( rs->has_true_var_solution() ) {
    double primal_bound = minimizing ? rs->get_true_ub() : rs->get_true_lb();
    if( improves( primal_bound , incumbentCell->load() , minimizing ) ) {
     bestBound = primal_bound;
     incumbentCell->store( bestBound );
     delete bestSol;
     bestSol = rs->get_true_solution();
     }
    }
   // the bound is recorded before the fate of the node is decided, so that
   // a node fenced by bound carries the bound that fenced it rather than the
   // one it was initialized with
   // the node may already have the (inherited) bound of its parent: the
   // bound is improved only if this one is better, while the branching is
   // done by the relaxation with the best bound, or by the first one
   auto relaxBound = minimizing ? s->get_lb() : s->get_ub();
   const bool better = minimizing ? relaxBound > currentNode->get_bound()
                                  : relaxBound < currentNode->get_bound();
   if( better )
    currentNode->set_bound( relaxBound );
   if( better || ( ! branchSolver ) )
    branchSolver = rs;
   if( cannot_improve( relaxBound , incumbentCell->load() , minimizing ,
                       relAcc , absAcc ) ) {
    toPrune = true;
    return( Solver::kOK );
    }
   }
  return( Solver::kOK );
  }

 // serial evaluation
 for( auto & [ rs , s ] : *f_RelaxationSolvers ) {
  auto z = s->compute();
  if( z == Solver::kInfeasible ) {
   if( ( ! currentNode->get_toFather() ) && ( ! has_incumbent ) )
    return( Solver::kInfeasible );       // the root is infeasible
   toPrune = true;
   currentNode->set_infeasible( true );
   currentNode->set_bound( minimizing ? Inf< double >()
                                           : - Inf< double >() );
   break;
   }
  if( z == Solver::kCutOff ) {  // beyond the incumbent [see above]
   fence_by_cutoff( currentNode , s , minimizing , toPrune );
   return( Solver::kOK );
   }
  if( z != ThinComputeInterface::kOK )
   return( z );

  // see if the primal bound improves thanks to a true solution
  if( rs->has_true_var_solution() ) {
   double primal_bound = minimizing ? rs->get_true_ub() : rs->get_true_lb();
   if( improves( primal_bound , incumbentCell->load() , minimizing ) ) {
    // in the parallel exploration the incumbent is shared between the
    // workers: re-check the improvement under the mutex
    std::unique_lock< std::mutex > guard;
    if( incumbentMutex )
     guard = std::unique_lock< std::mutex >( *incumbentMutex );
    if( improves( primal_bound , incumbentCell->load() , minimizing ) ) {
     bestBound = primal_bound;
     incumbentCell->store( bestBound );
     delete bestSol;
     bestSol = rs->get_true_solution();
     }
    }
   }

  // update the relaxation bound of the node and select the branching solver;
  // the bound is recorded before the fate of the node is decided [see above]
  auto relaxBound = minimizing ? s->get_lb() : s->get_ub();
  const bool better = minimizing ? relaxBound > currentNode->get_bound()
                                 : relaxBound < currentNode->get_bound();
  if( better )
   currentNode->set_bound( relaxBound );
  if( better || ( ! branchSolver ) )  // see above
   branchSolver = rs;
  if( cannot_improve( relaxBound , incumbentCell->load() , minimizing ,
                      relAcc , absAcc ) ) {
   toPrune = true;
   return( Solver::kOK );
   }
  }
 return( Solver::kOK );
 }

/*--------------------------------------------------------------------------*/
/// compute the heuristics at the current node
/** Computes every heuristic ChangeSolver at the current node, updating the
 * incumbent bestBound / bestSol when improved; sets \p toPrune when the
 * node can be discarded.
 *  @return the sol_type [see Solver.h] of the computation */

static int computeHeuristic(
                  std::vector< std::pair< ChangeSolver * ,
                               Solver * > > * f_HeuristicSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol ,
                  std::atomic< double > * incumbentCell , bool & toPrune ,
                  int nThreads = 1 , std::mutex * incumbentMutex = nullptr )
{
 // parallel evaluation of the (several) heuristics of this node, bit-identical
 // to the serial path below [see computeRelaxations()]
 if( ( nThreads > 1 ) && ( f_HeuristicSolvers->size() > 1 ) &&
     ( ! incumbentMutex ) ) {
  const std::size_t n = f_HeuristicSolvers->size();
  std::vector< int > zs( n , INT_MIN );   // INT_MIN: not computed (early stop)
  if( computeAllParallel( *f_HeuristicSolvers , zs , nThreads ) ) {
   if( ! currentNode->get_toFather() )    // the root is infeasible
    return( Solver::kInfeasible );
   toPrune = true;
   return( Solver::kOK );
   }
  for( std::size_t i = 0 ; i < n ; ++i ) {
   auto & [ hs , s ] = (*f_HeuristicSolvers)[ i ];
   if( zs[ i ] != ThinComputeInterface::kOK )
    return( zs[ i ] );
   double primal_bound = minimizing ? s->get_ub() : s->get_lb();
   if( s->has_var_solution() && s->is_var_feasible() &&
       improves( primal_bound , incumbentCell->load() , minimizing ) ) {
    bestBound = primal_bound;
    incumbentCell->store( bestBound );
    delete bestSol;
    bestSol = solution_of( s );
    }
   }
  return( Solver::kOK );
  }

 // serial evaluation
 for( auto & [ hs , s ] : *f_HeuristicSolvers ) {
  auto z = s->compute();
  if( z == Solver::kInfeasible ) {
   if( ! currentNode->get_toFather() )    // the root is infeasible
    return( Solver::kInfeasible );
   toPrune = true;
   break;
   }
  if( z != ThinComputeInterface::kOK )
   return( z );

  // see if the primal bound improves
  double primal_bound = minimizing ? s->get_ub() : s->get_lb();
  if( s->has_var_solution() && s->is_var_feasible() &&
      improves( primal_bound , incumbentCell->load() , minimizing ) ) {
   std::unique_lock< std::mutex > guard;
   if( incumbentMutex )
    guard = std::unique_lock< std::mutex >( *incumbentMutex );
   if( improves( primal_bound , incumbentCell->load() , minimizing ) ) {
    bestBound = primal_bound;
    incumbentCell->store( bestBound );
    delete bestSol;
    bestSol = solution_of( s );
    }
   }
  }
 return( Solver::kOK );
 }

/*--------------------------------------------------------------------------*/
/// run the compute() of a set of Solver in parallel over \p maxThreads
/** Runs s->compute() for every Solver \p s in \p slvrs concurrently on (up
 * to) \p maxThreads threads (the calling thread is one of them), storing the
 * return codes in \p zs. Only compute() is run here, never any result
 * extraction (get_lb() / get_true_solution() / ...): those touch the shared
 * Block and are left to the serial reduction in the caller. This is safe
 * because the inner Solver of the BranchAndXSolver are independent clones
 * with private state [see createWorkerSolvers() / the per-node Solver sets],
 * so their compute() do not race on the Block.
 *
 * As soon as one Solver returns kInfeasible the node is provably infeasible
 * (a relaxation with an empty feasible region proves the original one empty),
 * so the remaining Solver are not started: the entries of \p zs they would
 * have filled stay untouched. The method returns true in this case, telling
 * the caller to prune the node without reducing \p zs (which is incomplete);
 * when it returns false every Solver was computed and \p zs is complete. */

template< typename SolverPtr >
static bool computeAllParallel( const std::vector< SolverPtr > & slvrs ,
                                std::vector< int > & zs , int maxThreads )
{
 const int n = int( slvrs.size() );
 const int K = std::max( 1 , std::min( maxThreads , n ) );
 std::atomic< int > next( 0 );
 std::atomic< bool > infeasible( false );
 auto work = [ & ]() {
  for( int i ; ( ! infeasible.load( std::memory_order_relaxed ) ) &&
               ( ( i = next.fetch_add( 1 ) ) < n ) ; ) {
   zs[ i ] = slvrs[ i ].second->compute();
   if( zs[ i ] == Solver::kInfeasible )
    infeasible.store( true , std::memory_order_relaxed );
   }
  };
 std::vector< std::thread > pool;
 pool.reserve( K - 1 );
 for( int t = 0 ; t < K - 1 ; ++t )
  pool.emplace_back( work );
 work();                       // the calling thread is a worker too
 for( auto & th : pool )
  th.join();
 return( infeasible.load() );
 }

/*--------------------------------------------------------------------------*/
/// evaluate the root node and produce its branching list
/** The root is fenced with the same tolerances \p relAcc and \p absAcc as
 * any other node [see cannot_improve()]: with none, a root whose bound is
 * the value of its own true solution up to the accuracy of the relaxation
 * would be branched, with nothing left to branch upon. */

static int initializeRoot(
                  std::vector< std::pair< RelaxationSolver * ,
                               Solver * > > * f_RelaxationSolvers ,
                  std::vector< std::pair< ChangeSolver * ,
                               Solver * > > * f_HeuristicSolvers ,
                  Node * rootNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol ,
                  std::atomic< double > * incumbentCell ,
                  RelaxationSolver * & branchSolver ,
                  std::mutex & globalMutex ,
                  double relAcc , double absAcc )
{
 rootNode->initializeBound( minimizing );
 rootNode->set_evaluated( true );
 bool toPrune = false;
 int res = computeRelaxations( f_RelaxationSolvers , rootNode , minimizing ,
                               bestBound , bestSol , incumbentCell ,
                               toPrune , branchSolver , relAcc , absAcc ,
                               1 , &globalMutex );
 if( toPrune || ( res != ThinComputeInterface::kOK ) )
  return( res );
 res = computeHeuristic( f_HeuristicSolvers , rootNode , minimizing ,
                         bestBound , bestSol , incumbentCell , toPrune , 1 ,
                         &globalMutex );
 if( toPrune || ( res != ThinComputeInterface::kOK ) )
  return( res );
 rootNode->obtainBranchList( branchSolver );
 return( Solver::kOK );
 }

/*--------------------------------------------------------------------------*/
/// evaluate a freshly generated node: relaxations, heuristics, separation
/** The per-node work every exploration performs the same way, factored out so
 * the serial explore(), the parallel ramp-up and its workers share it instead
 * of each repeating it. It moves the :ChangeSolver onto \p node (storing its
 * undo Change), solves the relaxations and heuristics, and leaves the solvers
 * ON the node (the caller then keeps it and branches it, or prunes it and
 * moves them back). \p nThreads > 1 evaluates the several solvers of the node
 * concurrently (serial-tree explorations only); \p incumbentMutex != nullptr
 * serializes the incumbent updates (parallel-tree exploration). Sets \p toPrune
 * / \p branchSolver, and marks the node infeasible [see Node::set_infeasible()]
 * when its relaxation has no solution.
 *  @return the sol_type [see Solver.h] of the evaluation */

static int evaluateNode(
                  Node * node , std::list< ChangeSolver * > & solvers ,
                  std::vector< std::pair< RelaxationSolver * ,
                               Solver * > > * relaxation ,
                  std::vector< std::pair< ChangeSolver * ,
                               Solver * > > * heuristic ,
                  bool minimizing , double & bestBound , Solution * & bestSol ,
                  std::atomic< double > * incumbentCell ,
                  bool & toPrune , RelaxationSolver * & branchSolver ,
                  int nThreads ,
                  double relAcc , double absAcc , std::mutex & globalMutex ,
                  std::mutex * incumbentMutex = nullptr ,
                  bool moveToNode = true )
{
 if( moveToNode )                  // with the lazy protocol the solvers have
  moveSolverToSon( node , &solvers );  // already been moved to the node
 // the relaxation bound the node has inherited from its parent is valid within
 // this solve, the data being the same: the evaluation can only improve it
 // (a re-seeded frontier resets it, the data having changed, see
 // reseedFrontier())
 const double inherited = node->get_bound();
 node->initializeBound( minimizing );
 node->set_bound( inherited );
 node->set_evaluated( true );

 // computeRelaxations() / computeHeuristic() dispatch internally on nThreads:
 // a worker of the parallel tree exploration (incumbentMutex set) evaluates its
 // per-node solvers serially, otherwise nThreads > 1 runs them in parallel
 int res = computeRelaxations( relaxation , node , minimizing , bestBound ,
                               bestSol , incumbentCell , toPrune ,
                               branchSolver , relAcc ,
                               absAcc , nThreads , incumbentMutex );
 if( ( res != Solver::kOK ) || toPrune )
  return( res );

 return( computeHeuristic( heuristic , node , minimizing , bestBound , bestSol ,
                           incumbentCell , toPrune , nThreads ,
                           incumbentMutex ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OPEN-LIST DISCIPLINES -------------------------*/
/*--------------------------------------------------------------------------*/
// the three open-list disciplines, as thin adapters that give the standard
// container adapters a common interface [see OpenList]: a LIFO stack (depth-
// first), a FIFO queue (breadth-first), a relaxation-bound priority queue
// (best- first). The storage is entirely std::stack / std::queue /
// std::priority_- queue, there is no hand-rolled data structure here

namespace {

/// LIFO open list: a std::stack, giving depth-first exploration

class StackOpenList final : public OpenList {

 public:

 [[nodiscard]] bool empty( void ) const override { return( c.empty() ); }

 void push( ExploringNode * node ) override { c.push( node ); }

 ExploringNode * pop( void ) override {
  auto node = c.top();
  c.pop();
  return( node );
  }

 [[nodiscard]] bool isLIFO( void ) const override { return( true ); }

 private:

 std::stack< ExploringNode * > c;

 };  // end( class( StackOpenList ) )

/*--------------------------------------------------------------------------*/
/// FIFO open list: a std::queue, giving breadth-first exploration

class QueueOpenList final : public OpenList {

 public:

 [[nodiscard]] bool empty( void ) const override { return( c.empty() ); }

 void push( ExploringNode * node ) override { c.push( node ); }

 ExploringNode * pop( void ) override {
  auto node = c.front();
  c.pop();
  return( node );
  }

 private:

 std::queue< ExploringNode * > c;

 };  // end( class( QueueOpenList ) )

/*--------------------------------------------------------------------------*/
/// relaxation-bound priority open list: a std::priority_queue, best-first

class PriorityOpenList final : public OpenList {

 using Cmp = std::function< bool( ExploringNode * , ExploringNode * ) >;

 public:

 explicit PriorityOpenList( Cmp cmp ) : c( std::move( cmp ) ) {}

 [[nodiscard]] bool empty( void ) const override { return( c.empty() ); }

 void push( ExploringNode * node ) override { c.push( node ); }

 ExploringNode * pop( void ) override {
  auto node = c.top();
  c.pop();
  return( node );
  }

 private:

 std::priority_queue< ExploringNode * , std::vector< ExploringNode * > , Cmp >
  c;

 };  // end( class( PriorityOpenList ) )

/*--------------------------------------------------------------------------*/
/// best-first open list with depth-first dives
/** The global frontier is a relaxation-bound priority queue, but once a node
 * is taken the search dives straight into its most promising child down to a
 * leaf, leaving the siblings to the global frontier. The dive reaches a
 * complete (feasible) solution quickly, so the incumbent tightens early and
 * the pruning bites sooner. The children of the just-expanded node arrive
 * through push() before the next pop(): pop() routes the most promising of
 * them onward (continuing the dive) and spills the rest to the priority queue;
 * when no child arrives the dive has bottomed out and the next global best is
 * taken. */

class DiveOpenList final : public OpenList {

 using Cmp = std::function< bool( ExploringNode * , ExploringNode * ) >;

 public:

 explicit DiveOpenList( Cmp cmp ) : best( std::move( cmp ) ) {}

 [[nodiscard]] bool empty( void ) const override {
  return( best.empty() && children.empty() );
  }

 void push( ExploringNode * node ) override { children.push_back( node ); }

 ExploringNode * pop( void ) override {
  if( ! children.empty() ) {
   // the loop is isLIFO(), so it pushed the children in reverse branching
   // order: back() is the first (most promising) branch, the dive follows it
   auto node = children.back();
   children.pop_back();
   for( auto sibling : children )
    best.push( sibling );
   children.clear();
   return( node );
   }
  auto node = best.top();
  best.pop();
  return( node );
  }

 [[nodiscard]] bool isLIFO( void ) const override { return( true ); }

 private:

 std::priority_queue< ExploringNode * , std::vector< ExploringNode * > , Cmp >
  best;

 std::vector< ExploringNode * > children;  // pushed since the last pop()

 };  // end( class( DiveOpenList ) )

}  // anonymous namespace

/*--------------------------------------------------------------------------*/
/*------------------- THE TREE AS A PROBLEM OF PBeFS -----------------------*/
/*--------------------------------------------------------------------------*/
// the tree of the BranchAndXSolver as a BranchAndBound problem of
// PBeFS [see PBeFSSolve()]: the classes live here, and not in
// the member function, since the evaluator has a member template

namespace {

namespace PB = PBeFS_di_unipi_it;

/// a node of the tree for the branch-and-bound of PBeFS
/** The Change that generates the node from its parent, the node of the
 * parent and the depth; its evaluation (the relaxation bound, whether it is
 * pruned, i.e., infeasible or fenced, and its branching Changes) is filled
 * once by the worker that evaluates it, before the node is handed to any
 * other worker. The Changes are plain data, applied by the Solver of any
 * worker [see ParallelDFSSolve()]. */

struct BXLink {
 std::shared_ptr< Change > change;              ///< from the parent
 std::shared_ptr< BXLink > parent;              ///< nullptr at the root
 int depth = 0;                                 ///< the depth
 double bound = 0;                              ///< in the sense of the Block
 /// the estimate of the best value under it, in the sense of the Block,
 /// NaN if none [see BranchEstimator]
 double estimate = std::numeric_limits< double >::quiet_NaN();
 bool evaluated = false;                        ///< evaluated already
 bool pruned = false;                           ///< infeasible or fenced
 std::vector< std::shared_ptr< Change > > branches;  ///< to the children
 std::vector< double > estimates;               ///< of the children, if any
 /// the value of the feasible solution its evaluation has found, in the
 /// sense of the Block, NaN if none, and the solution (only by synchronous
 /// rounds [see BXContext::deterministic])
 double found = std::numeric_limits< double >::quiet_NaN();
 std::unique_ptr< Solution > solution;          ///< that solution
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// what the workers share: the data of the BranchAndXSolver

struct BXContext {

 using Relaxations = std::vector< std::pair< RelaxationSolver * ,
                                             Solver * > >;
 using Heuristics = std::vector< std::pair< ChangeSolver * , Solver * > >;

 BXContext( double & bb , Solution * & bs , std::atomic< double > * ic ,
            double rt , double at , bool min , bool e , bool l , bool d ,
            std::mutex & gm )
  : bestBound( bb ) , bestSolution( bs ) , incumbentCell( ic ) ,
    relTol( rt ) , absTol( at ) , minimizing( min ) , eager( e ) ,
    lifo( l ) , deterministic( d ) , globalMutex( gm ) {}

 double & bestBound;                     ///< the incumbent value
 Solution * & bestSolution;              ///< the incumbent solution
 std::atomic< double > * incumbentCell;  ///< the incumbent, lock-free
 double relTol;                          ///< the relative tolerance
 double absTol;                          ///< the absolute tolerance
 bool minimizing;                        ///< the sense of the Block
 bool eager;                             ///< the eager bounding protocol
 bool lifo;                              ///< the first branch last
 /// by synchronous rounds: each worker keeps the incumbent of the round
 /// [see BXWorker::cutoff()] and the solutions stay in the nodes
 bool deterministic;
 std::mutex & globalMutex;               ///< of the BranchAndXSolver
 std::mutex incumbentMutex;              ///< serializes the incumbent
 std::vector< Relaxations * > relaxation;  ///< of each worker
 std::vector< Heuristics * > heuristic;    ///< of each worker
 std::atomic< int > error{ Solver::kOK };  ///< the first error, if any

 /// the sign that turns a value of the Block into one to minimize
 [[nodiscard]] double sign( void ) const { return( minimizing ? 1 : -1 ); }

 };  // end( struct( BXContext ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the evaluator of a worker: its Solver, and the path of the node they
/// sit on, with the undo Changes of its own

class BXWorker {

 public:

 BXWorker( BXContext & c , unsigned i )
  : f_c( c ) , f_i( i ) , f_best( c.sign() * Inf< double >() ) ,
    f_cell( f_best ) {
  for( auto & p : * c.relaxation[ i ] )
   f_solvers.push_back( p.first );
  for( auto & p : * c.heuristic[ i ] )
   f_solvers.push_back( p.first );
  }

 BXWorker( BXWorker && w ) noexcept
  : f_c( w.f_c ) , f_i( w.f_i ) , f_solvers( std::move( w.f_solvers ) ) ,
    v_path( std::move( w.v_path ) ) , f_best( w.f_best ) ,
    f_sol( w.f_sol ) , f_cell( w.f_cell.load() ) {
  w.v_path.clear();
  w.f_sol = nullptr;
  }

 /// back to the root, so that the Solver are where they started
 ~BXWorker() {
  move_to( nullptr );
  delete f_sol;
  }

 /// the incumbent z of the search, in the sense of PBeFS, which by
 /// synchronous rounds is the one the worker evaluates against
 void cutoff( double z ) {
  if( ! f_c.deterministic )
   return;
  f_best = f_c.sign() * z;
  f_cell.store( f_best );
  }

 /// the evaluation of the node n, done now unless already done; by
 /// synchronous rounds, the feasible value is the one found by the node
 PB::Evaluation evaluate( std::shared_ptr< BXLink > & n ) {
  if( ! n->evaluated )
   evaluate_here( n );
  PB::Evaluation e;
  e.bound = f_c.sign() * n->bound;
  if( f_c.deterministic ) {
   e.infeasible = n->pruned && std::isnan( n->found );
   if( ! std::isnan( n->found ) )
    e.feasible = f_c.sign() * n->found;
   return( e );
   }
  e.infeasible = n->pruned;
  const double z = f_c.incumbentCell->load();
  if( std::isfinite( z ) )
   e.feasible = f_c.sign() * z;
  return( e );
  }

 /// the children of the evaluated node n, evaluated now with the eager
 /// protocol (and then given only if they survive), the first branch
 /// given last when the open set is a stack
 template< class F >
 void branch( const std::shared_ptr< BXLink > & n , F && f ) {
  std::vector< std::shared_ptr< BXLink > > kids;
  for( std::size_t i = 0 ; i < n->branches.size() ; ++i ) {
   auto k = std::make_shared< BXLink >();
   k->change = n->branches[ i ];
   k->parent = n;
   k->depth = n->depth + 1;
   k->bound = n->bound;                 // inherited
   if( i < n->estimates.size() )
    k->estimate = n->estimates[ i ];
   if( f_c.eager ) {
    evaluate_here( k );
    if( k->pruned && std::isnan( k->found ) )
     continue;
    }
   kids.push_back( std::move( k ) );
   }
  n->branches.clear();
  n->estimates.clear();
  if( f_c.lifo )
   std::reverse( kids.begin() , kids.end() );
  for( auto & k : kids ) {
   const double b = f_c.sign() * k->bound;
   f( std::move( k ) , b );
   }
  }

 private:

 // evaluates the node n [see evaluateNode()], and asks for its branching
 // Changes while the Solver are on it
 void evaluate_here( const std::shared_ptr< BXLink > & n ) {
  move_to( n );
  DFSNode tmp( nullptr );
  tmp.set_bound( n->bound );
  // a node with no undo Change is the root to the evaluation, which tells
  // an infeasible root from any other infeasible node: the node borrows the
  // undo Change of the worker, and gives it back before it dies
  if( n->change )
   tmp.set_toFather( v_path.back().second );
  bool toPrune = false;
  RelaxationSolver * bs = nullptr;
  // by synchronous rounds, the incumbent of the worker [see cutoff()]
  const bool d = f_c.deterministic;
  std::atomic< double > * cell = d ? & f_cell : f_c.incumbentCell;
  const int res = evaluateNode( & tmp , f_solvers , f_c.relaxation[ f_i ] ,
                                f_c.heuristic[ f_i ] , f_c.minimizing ,
                                d ? f_best : f_c.bestBound ,
                                d ? f_sol : f_c.bestSolution , cell ,
                                toPrune , bs , 1 , f_c.relTol , f_c.absTol ,
                                f_c.globalMutex ,
                                d ? nullptr : & f_c.incumbentMutex , false );
  tmp.set_toFather( nullptr );
  n->evaluated = true;
  n->bound = tmp.get_bound();
  if( f_sol ) {                         // found by this node
   n->solution.reset( f_sol );
   n->found = f_best;
   f_sol = nullptr;
   }
  if( res != Solver::kOK ) {
   // an infeasible root is no error: the search ends with no incumbent,
   // and compute() tells that the problem is infeasible
   if( ( res != Solver::kInfeasible ) || n->change ) {
    int ok = Solver::kOK;
    f_c.error.compare_exchange_strong( ok , res );
    }
   n->pruned = true;
   return;
   }
  n->pruned = toPrune || tmp.is_infeasible() ||
              cannot_improve( n->bound , cell->load() ,
                              f_c.minimizing , f_c.relTol , f_c.absTol );
  if( ( ! n->pruned ) && bs ) {
   for( Change * c : bs->branch() )
    n->branches.emplace_back( c );
   if( auto be = dynamic_cast< BranchEstimator * >( bs ) ) {
    n->estimates = be->branch_estimates();
    if( n->estimates.size() != n->branches.size() )
     n->estimates.clear();
    }
   }
  }

 // moves the Solver from their node to n (nullptr is the root): up to the
 // deepest common ancestor by the undo Changes, then down by the Changes,
 // the first Solver giving the undo Change [see moveBetweenNodes()]
 void move_to( const std::shared_ptr< BXLink > & n ) {
  std::vector< std::shared_ptr< BXLink > > target;
  for( auto x = n ; x && x->change ; x = x->parent )
   target.push_back( x );
  std::reverse( target.begin() , target.end() );
  std::size_t k = 0;
  while( ( k < v_path.size() ) && ( k < target.size() ) &&
         ( v_path[ k ].first == target[ k ] ) )
   ++k;
  while( v_path.size() > k ) {
   Change * undo = v_path.back().second;
   for( auto s : f_solvers )
    s->apply( undo , false );
   delete undo;
   v_path.pop_back();
   }
  for( ; k < target.size() ; ++k ) {
   Change * undo = nullptr;
   bool first = true;
   for( auto s : f_solvers ) {
    if( first )
     undo = s->apply( target[ k ]->change.get() , true );
    else
     s->apply( target[ k ]->change.get() , false );
    first = false;
    }
   v_path.emplace_back( target[ k ] , undo );
   }
  }

 BXContext & f_c;                               ///< the shared data
 unsigned f_i;                                  ///< the index of the worker
 std::list< ChangeSolver * > f_solvers;         ///< its Solver
 /// the nodes from the root to that of the Solver, with the undo Changes
 std::vector< std::pair< std::shared_ptr< BXLink > , Change * > > v_path;
 /// by synchronous rounds, the incumbent of the worker [see cutoff()]: its
 /// value, the solution found by the node being evaluated, if any, and the
 /// value as the cell that evaluateNode() reads
 double f_best;
 Solution * f_sol = nullptr;
 std::atomic< double > f_cell;

 };  // end( class( BXWorker ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// the tree of the BranchAndXSolver as a BranchAndBound problem; the values
/// are those of the Block times BXContext::sign(), so that it minimizes

struct BXTree {

 using node_type = std::shared_ptr< BXLink >;
 static constexpr PB::Sense sense = PB::Sense::minimize;

 BXContext & c;                                 ///< the shared data

 [[nodiscard]] node_type root( void ) const {
  auto r = std::make_shared< BXLink >();
  r->bound = c.minimizing ? - Inf< double >() : Inf< double >();
  return( r );
  }

 [[nodiscard]] BXWorker make_worker( unsigned i ) const {
  return( BXWorker( c , i ) );
  }

 /// the estimate of the best value under n, which orders the nodes of a
 /// best-estimate search: the one given by the RelaxationSolver that
 /// branched, never better than the bound of n, or the bound if none
 [[nodiscard]] double estimate( const node_type & n ) const {
  double e = n->bound;
  if( ! std::isnan( n->estimate ) )
   e = c.minimizing ? std::max( n->estimate , n->bound )
                    : std::min( n->estimate , n->bound );
  return( c.sign() * e );
  }

 };  // end( struct( BXTree ) )

}  // anonymous namespace

/*--------------------------------------------------------------------------*/
/*----------------- METHODS OF BranchAndXSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

double BranchAndXSolver::no_incumbent( void ) const
{
 return( f_Block->get_objective_sense() == Objective::eMax
	 ? - Inf< double >() : Inf< double >() );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue BranchAndXSolver::get_lb( void )
{
 if( f_Block->get_objective_sense() == Objective::eMax )
  return( bestBound );                 // the incumbent is a lower bound
 // the whole tree explored, which a proof of infeasibility also requires,
 // proves the incumbent; otherwise the nodes left open bound it
 if( ( f_state == kOK ) || ( f_state == kInfeasible ) )
  return( bestBound );
 return( std::min( f_open_bound , bestBound ) );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue BranchAndXSolver::get_ub( void )
{
 if( f_Block->get_objective_sense() == Objective::eMin )
  return( bestBound );                 // the incumbent is an upper bound
 // the whole tree explored, which a proof of infeasibility also requires,
 // proves the incumbent; otherwise the nodes left open bound it
 if( ( f_state == kOK ) || ( f_state == kInfeasible ) )
  return( bestBound );
 return( std::max( f_open_bound , bestBound ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int BranchAndXSolver::compute( bool changedvars )
{
 lock();
 process_outstanding_Modification();
 int old_state = f_state;
 if( f_state == kStillRunning )
  return( kError );
 f_state = kStillRunning;
 // no bound until an exploration says otherwise [see explore()]
 f_open_bound = ( f_Block->get_objective_sense() == Objective::eMax )
                ? Inf< double >() : - Inf< double >();
 std::mutex globalMutex;

 // per-solve budgets and tolerances from the inherited standard parameters
 nodeBudget = get_int_par( intMaxNodes );
 timeBudget = get_dbl_par( dblMaxTime );
 relTol = get_dbl_par( dblRelAcc );
 absTol = get_dbl_par( dblAbsAcc );

 // the tree is retained by any serial exploration, but not by the parallel
 // depth-first one, whose worker forest is not kept [see ParallelDFSSolve()]
 const bool parallelDFS = ( solveType == DFS ) &&
                          ( get_int_par( intMaxThread ) > 1 );
 // the branch-and-bound of PBeFS, which retains no tree either
 const bool astar = ( treeExplorer != OwnExplorer );

 // incumbent-dependent local fixing folded into the branching Changes is
 // unsafe when the tree is retained across re-solves with different
 // incumbents: forbid it in that case, allow it otherwise [see
 // GlobalInformation::local_fixing_allowed()]
 f_lfaCell->store( ! ( ( reoptimize > 0 ) && ( ! parallelDFS ) &&
                       ( ! astar ) ) );

 // nothing changed since the last solve, which has ended: its answer
 // stands; one stopped by a budget goes on instead [see Solver::kStopTime]
 if( ( changes == 0 ) && ( ( old_state == kOK ) ||
                           ( old_state == kInfeasible ) ||
                           ( old_state == kUnbounded ) ) )
  f_state = old_state;
 else {
  // the outstanding changes say what the Modification have changed [see
  // process_outstanding_Modification()], where one that the tree does not
  // survive [see keeps_tree()] has made them anything: when they are not, a
  // serial re-solve with retention enabled reoptimizes out of the fenced
  // frontier of the previous tree rather than solving from scratch [see
  // intReoptimize / reseedFrontier()]; otherwise the retained tree (if any)
  // is discarded; the test is not keeps_tree( changes ), since the bits of
  // different Modification or-ed together (an added Variable and an added
  // Constraint, say) may read as one that the tree does not survive
  if( ( changes == Modification::eModAll ) || ( ! reoptimize ) ||
      parallelDFS || astar )
   discardRetainedTree();

  bestBound = no_incumbent();
  // no incumbent yet [see no_incumbent() and GlobalInformation::str_Incumbent]
  f_incumbentCell->store( bestBound );
  // one exploration, parameterized by the open set that the strategy asks
  // for [see treeSolve()]; the only separate one is the parallel depth-first,
  // which runs its own workers [see ParallelDFSSolve()]
  if( astar )
   f_state = PBeFSSolve( globalMutex , get_int_par( intMaxThread ) );
  else
   f_state = parallelDFS ? ParallelDFSSolve( globalMutex ,
                                             get_int_par( intMaxThread ) )
                         : treeSolve( globalMutex );

  // an exploration that has been completed without ever finding a feasible
  // solution proves that there is none
  if( ( f_state == kOK ) && ( bestBound == no_incumbent() ) )
   f_state = kInfeasible;
  }

 changes = 0;
 unlock();
 return( f_state );

 }  // end( BranchAndXSolver::compute )


/*--------------------------------------------------------------------------*/

void BranchAndXSolver::discardTree( ExploringNode * root , OpenList & open ,
                                    std::list< ChangeSolver * > * solvers ,
                                    const ExploringNode * current )
{
 // the global relaxation bound: the best among those of the nodes left open
 // and of the one being evaluated, whose region covers that of its children;
 // with none of them there is no bound
 const bool minimizing =
  ( f_Block->get_objective_sense() != Objective::eMax );
 std::vector< double > bounds;
 if( current )
  bounds.push_back( current->get_bound() );

 // the nodes still in the open set are also children in the tree, so the
 // recursive deletion of the tree covers them; emptying the open set first
 // (without deleting) avoids a double delete
 while( ! open.empty() )
  bounds.push_back( open.pop()->get_bound() );

 if( bounds.empty() )
  f_open_bound = minimizing ? - Inf< double >() : Inf< double >();
 else
  f_open_bound = minimizing ? *std::min_element( bounds.begin() ,
                                                 bounds.end() )
                            : *std::max_element( bounds.begin() ,
                                                 bounds.end() );
 std::function< void( ExploringNode * ) > deleteTree =
  [ & ]( ExploringNode * node ) {
   if( ! node )
    return;
   for( auto child : node->get_children() )
    deleteTree( child );
   node->get_children().clear();
   delete node;
   };
 deleteTree( root );
 delete solvers;
 f_treeRoot = nullptr;
 subOptimalNodes.clear();
 infeasibleNodes.clear();
 integerNodes.clear();
 }

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::reseedFrontier( OpenList & open ,
                                     std::list< ChangeSolver * > * solvers ,
                                     bool minimizing ,
                                     ExploringNode * & rootNode ,
                                     ExploringNode * & currentNode )
{
 // the retained tree is taken back: compute() guarantees it is only reused
 // under changes that leave it valid [see keeps_tree()], so its interior is
 // NOT re-derived, only its fenced frontier is re-evaluated
 rootNode = f_treeRoot;
 f_treeRoot = nullptr;              // ownership back to this solve
 currentNode = rootNode;

 std::list< ExploringNode * > frontier;
 frontier.swap( subOptimalNodes );
 if( Modification::may_grow_region( changes ) ) {
  // only a feasible region that may grow can un-fence an infeasible node:
  // otherwise (say, when only the Objective changes) that part of the
  // frontier stays fenced, with no re-evaluation at all
  frontier.splice( frontier.end() , infeasibleNodes );
  infeasibleNodes.clear();
  }

 // best (previous) bound first: the optimum almost surely lives in the first
 // few nodes, so the incumbent warms up immediately and the rest of the
 // frontier mostly just re-fences
 frontier.sort( [ minimizing ]( ExploringNode * a , ExploringNode * b ) {
  return( minimizing ? a->get_bound() < b->get_bound()
                     : a->get_bound() > b->get_bound() );
  } );

 int res = Solver::kOK;
 std::vector< ExploringNode * > reopened;  // pushed after the loop
 for( auto F : frontier ) {
  ExploringNode::moveBetweenNodes( currentNode , F , solvers );
  currentNode = F;
  F->initializeBound( minimizing );
  bool toPrune = false;
  RelaxationSolver * nodeBranchSolver = nullptr;
  res = computeRelaxations( &f_RelaxationSolvers , F , minimizing , bestBound ,
                            bestSolution , f_incumbentCell , toPrune ,
                            nodeBranchSolver ,
                            relTol , absTol , maxThreadForSolvers );
  if( ( res == Solver::kOK ) && ( ! toPrune ) )
   res = computeHeuristic( &f_HeuristicSolvers , F , minimizing , bestBound ,
                           bestSolution , f_incumbentCell , toPrune ,
                           maxThreadForSolvers );
  if( res != Solver::kOK )
   break;
  if( ( ! toPrune ) &&
      ( ! cannot_improve( F->get_bound() , f_incumbentCell->load() ,
                          minimizing , relTol , absTol ) ) ) {
   // re-opened: the branching Changes of the previous solve, if any, are
   // still valid (any branching is), so they are reused rather than leaked
   if( F->getBranches().empty() )
    F->obtainBranchList( nodeBranchSolver );
   reopened.push_back( F );
   }
  else if( F->is_infeasible() )      // fenced again, by reason
   infeasibleNodes.push_back( F );
  else
   subOptimalNodes.push_back( F );
  }

 if( res != Solver::kOK ) {          // leave the solvers at the root
  ExploringNode::moveBetweenNodes( currentNode , rootNode , solvers );
  currentNode = rootNode;
  return( res );
  }

 // hand the re-opened nodes to the open set so that the most promising is
 // explored first whatever the discipline: a LIFO stack receives them in
 // reverse, the others in bound order [see explore()]
 if( open.isLIFO() )
  for( auto it = reopened.rbegin() ; it != reopened.rend() ; ++it )
   open.push( *it );
 else
  for( auto n : reopened )
   open.push( n );

 return( Solver::kOK );

 }  // end( BranchAndXSolver::reseedFrontier )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::explore( OpenList & open , std::mutex & globalMutex ,
                               std::list< ChangeSolver * > * solvers ,
                               bool minimizing , int & counter ,
                               ExploringNode * rootNode ,
                               ExploringNode * currentNode , bool retain ,
                               std::chrono::high_resolution_clock::time_point
                                                                       start )
{
 RelaxationSolver * branchSolver = nullptr;
 int res = Solver::kOK;
 ExploringNode * oldNode = nullptr;
 int iterations = 0;      // extractions from the open set, for the log

 while( ( ! open.empty() ) && ( nodeBudget > 1 ) &&
        ( timeBudget > std::chrono::duration< double >(
           std::chrono::high_resolution_clock::now() - start ).count() ) ) {
  nodeBudget--;
  ++iterations;
  oldNode = currentNode;
  currentNode = open.pop();
  if( currentNode->get_parent() )   // the root: the solvers are already there
   ExploringNode::moveBetweenNodes( oldNode , currentNode , solvers );

  // with the lazy protocol the node is evaluated now, upon extraction,
  // rather than when it was created [see intBoundingProtocol]; a node is
  // evaluated at most once, whichever the protocol
  bool prunedHere = false;
  if( ! currentNode->is_evaluated() ) {
   auto evalStart = std::chrono::high_resolution_clock::now();
   res = evaluateNode( currentNode , *solvers , &f_RelaxationSolvers ,
                       &f_HeuristicSolvers , minimizing , bestBound ,
                       bestSolution , f_incumbentCell , prunedHere ,
                       branchSolver ,
                       maxThreadForSolvers , relTol , absTol , globalMutex ,
                       nullptr , false );
   if( res != Solver::kOK ) {       // an infeasible root ends up here
    // leave the solvers at the root, the Changes of the node undone
    ExploringNode::moveBetweenNodes( currentNode , rootNode , solvers );
    discardTree( rootNode , open , solvers , currentNode );
    return( res );
    }
   if( ! prunedHere )
    currentNode->obtainBranchList( branchSolver );
   logNode( f_log , iterations , currentNode , bestBound ,
	    std::chrono::duration< double >(
	     std::chrono::high_resolution_clock::now() - start ).count() ,
	    std::chrono::duration< double >(
	     std::chrono::high_resolution_clock::now() - evalStart ).count() ,
	    prunedHere || cannot_improve( currentNode->get_bound() ,
					  f_incumbentCell->load() , minimizing , relTol ,
					  absTol ) );
   }

  // branching and evaluation of the new children
  if( ( ! prunedHere ) &&
      ( ! cannot_improve( currentNode->get_bound() , f_incumbentCell->load() ,
                          minimizing , relTol , absTol ) ) ) {
   auto & branches = currentNode->getBranches();
   std::vector< ExploringNode * > kept;  // survivors, pushed after the loop
   for( std::size_t i = 0 ; i < branches.size() ; ++i ) {
    ExploringNode * new_node = new ExploringNode( branches[ i ] ,
                                                  currentNode ,
                                                  currentNode->get_level()
                                                  + 1 , ++counter );
    // the relaxation bound of the parent is valid for the child, whose region
    // is a part of that of the parent: the child starts from it, and its own
    // evaluation can only improve it; the estimate is the one given by the
    // RelaxationSolver that branched, if any [see BranchEstimator]
    new_node->set_bound( currentNode->get_bound() );
    new_node->set_estimate( currentNode->child_estimate( i ) );
    if( boundingProtocol == Lazy ) {
     // the child is not evaluated now: it goes into the open set carrying the
     // relaxation bound of its parent, and will be evaluated if and when it is
     // extracted [see intBoundingProtocol]
     currentNode->get_children().push_back( new_node );
     // the child now owns the branching Change as its f_change: null the
     // entry so that ~Node does not double-delete it on teardown
     *( std::find( branches.begin() , branches.end() ,
                   new_node->get_f_change() ) ) = nullptr;
     kept.push_back( new_node );
     continue;
     }
    bool toPrune = false;
    auto evalStart = std::chrono::high_resolution_clock::now();
    res = evaluateNode( new_node , *solvers , &f_RelaxationSolvers ,
                        &f_HeuristicSolvers , minimizing , bestBound ,
                        bestSolution , f_incumbentCell , toPrune ,
                        branchSolver ,
                        maxThreadForSolvers , relTol , absTol , globalMutex );
    if( res != Solver::kOK ) {
     // leave the solvers at the root, the Changes of the child undone
     ExploringNode::moveBetweenNodes( new_node , rootNode , solvers );
     discardTree( rootNode , open , solvers , currentNode );
     return( res );
     }
    logNode( f_log , iterations , new_node , bestBound ,
	     std::chrono::duration< double >(
	      std::chrono::high_resolution_clock::now() - start ).count() ,
	     std::chrono::duration< double >(
	      std::chrono::high_resolution_clock::now() - evalStart ).count() ,
	     toPrune || cannot_improve( new_node->get_bound() ,
					f_incumbentCell->load() , minimizing , relTol ,
					absTol ) );
    // keep the node only if its relaxation bound can improve the incumbent
    if( ( ! toPrune ) &&
        ( ! cannot_improve( new_node->get_bound() , f_incumbentCell->load() ,
                            minimizing , relTol , absTol ) ) ) {
     new_node->obtainBranchList( branchSolver );
     currentNode->get_children().push_back( new_node );
     // the kept child now owns the branching Change as its f_change: null
     // the entry so that ~Node does not double-delete it on teardown
     *( std::find( branches.begin() , branches.end() ,
                   new_node->get_f_change() ) ) = nullptr;
     moveSolverToFather( new_node , solvers );
     kept.push_back( new_node );  // pushed to the open set after the loop
     }
    else {                       // fenced or pruned child
     moveSolverToFather( new_node , solvers );
     *( std::find( branches.begin() , branches.end() ,
                   new_node->get_f_change() ) ) = nullptr;
     if( retain ) {              // it belongs to the fenced frontier:
      currentNode->get_children().push_back( new_node );
      if( new_node->is_infeasible() )  // infeasibility certificates survive
       infeasibleNodes.push_back( new_node );  //  a region that cannot grow
      else
       subOptimalNodes.push_back( new_node );
      }
     else
      delete new_node;
     }
    }
   // hand the surviving children to the open set so that the most promising
   // (the first from branch()) is explored first whatever the discipline:
   // a LIFO stack receives them in reverse, the others in branching order
   if( open.isLIFO() )
    for( auto it = kept.rbegin() ; it != kept.rend() ; ++it )
     open.push( *it );
   else
    for( auto n : kept )
     open.push( n );
   branches.erase( std::remove( branches.begin() , branches.end() ,
                                nullptr ) , branches.end() );
   currentNode->clear_estimates();   // they went with their children
   if( currentNode->get_children().empty() ) {
    if( retain )                 // fathomed leaf: fenced frontier
     subOptimalNodes.push_back( currentNode );
    if( currentNode->get_toFather() ) {   // the root has nowhere to climb
     if( retain ) {
      for( const auto s : *solvers )
       s->apply( currentNode->get_toFather() , false );
      currentNode = currentNode->get_parent();
      }
     else
      currentNode = ExploringNode::prune( currentNode , solvers );
     }
    }
   }
  else {
   // fenced at pop: it belongs to the fenced frontier, the root included -
   // a retained tree whose root is fathomed would otherwise have an empty
   // frontier, and the next re-solve would have nothing to re-explore
   if( retain )
    subOptimalNodes.push_back( currentNode );
   if( currentNode->get_toFather() ) {   // the root has nowhere to climb
    if( retain ) {
     for( const auto s : *solvers )
      s->apply( currentNode->get_toFather() , false );
     currentNode = currentNode->get_parent();
     }
    else
     currentNode = ExploringNode::prune( currentNode , solvers );
    }
   }
  }

 // the nodes left open when a budget stops the loop are taken out of the
 // open set: the exploration is complete only if none of them can improve
 // the incumbent, whatever the budgets that remain
 bool complete = true;
 std::vector< ExploringNode * > leftOpen;
 // the global relaxation bound: the best among those of the open nodes
 f_open_bound = minimizing ? Inf< double >() : - Inf< double >();
 while( ! open.empty() ) {
  auto node = open.pop();
  if( ! cannot_improve( node->get_bound() , f_incumbentCell->load() ,
                        minimizing , relTol , absTol ) )
   complete = false;
  f_open_bound = minimizing ? std::min( f_open_bound ,
                                        node->get_bound() )
                            : std::max( f_open_bound ,
                                        node->get_bound() );
  leftOpen.push_back( node );
  }

 // move the :ChangeSolver back to the root
 while( currentNode->get_toFather() ) {
  for( const auto s : *solvers )
   s->apply( currentNode->get_toFather() , false );
  currentNode = currentNode->get_parent();
  }

 if( retain ) {
  // retain the tree for future reoptimizations: the nodes left in the open
  // set (early stops) are open work, hence part of the frontier to re-seed
  subOptimalNodes.insert( subOptimalNodes.end() , leftOpen.begin() ,
                          leftOpen.end() );
  f_treeRoot = rootNode;
  }
 else {
  // the nodes left in the open set are also children in the tree: the
  // recursive deletion of the tree covers them
  std::function< void( ExploringNode * ) > deleteTree =
   [ & ]( ExploringNode * node ) {
    if( ! node )
     return;
    for( auto child : node->get_children() )
     deleteTree( child );
    node->get_children().clear();
    delete node;
    };
  deleteTree( rootNode );
  }
 delete solvers;

 // open nodes that can improve are left only when the loop has run out of
 // one of the budgets: the node one is gone when a single node is left
 if( ! complete )
  return( nodeBudget <= 1 ? Solver::kStopIter : Solver::kStopTime );
 return( kOK );

 }  // end( BranchAndXSolver::explore )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::treeSolve( std::mutex & globalMutex )
{
 auto * solvers = new std::list< ChangeSolver * >();
 bool minimizing;
 initializeVariables( solvers , &f_RelaxationSolvers , &f_HeuristicSolvers ,
                      minimizing );

 auto start = std::chrono::high_resolution_clock::now();
 int counter = 0;

 // the open set *is* the exploration strategy: a LIFO stack explores
 // depth-first, a FIFO queue breadth-first, a relaxation-bound priority queue
 // best-first, and the dive variant wraps the latter so that each best node
 // is followed depth-first down to a leaf [see DiveOpenList]; best-estimate
 // is a priority queue on the estimates [see Node::estimate_key()]
 auto cmp = [ minimizing ]( ExploringNode * a , ExploringNode * b ) {
  return( minimizing ? a->get_bound() > b->get_bound()
                     : a->get_bound() < b->get_bound() );
  };
 auto ecmp = [ minimizing ]( ExploringNode * a , ExploringNode * b ) {
  const double ka = a->estimate_key( minimizing );
  const double kb = b->estimate_key( minimizing );
  return( minimizing ? ka > kb : ka < kb );
  };
 std::unique_ptr< OpenList > openPtr;
 switch( solveType ) {
  case( DFS ):
   openPtr.reset( new StackOpenList() );
   break;
  case( BFS ):
   openPtr.reset( new QueueOpenList() );
   break;
  case( BestFS ):
   openPtr.reset( new PriorityOpenList( cmp ) );
   break;
  case( BestFSDive ):
   openPtr.reset( new DiveOpenList( cmp ) );
   break;
  case( BestEstimate ):
   openPtr.reset( new PriorityOpenList( ecmp ) );
   break;
  default:
   delete solvers;
   throw( std::invalid_argument( "BranchAndXSolver::treeSolve: invalid "
                                 "intSolveMethod" ) );
  }
 OpenList & open = *openPtr;

 logHeader( f_log );

 ExploringNode * rootNode = new ExploringNode( nullptr , nullptr , 0 );
 RelaxationSolver * branchSolver = nullptr;
 // the tree is retained for reoptimization whatever the strategy [see
 // intReoptimize]: a re-solve re-seeds from the frontier of the retained
 // tree, in the order that the open set of this very strategy dictates
 const bool retain = ( reoptimize > 0 );
 int res;
 ExploringNode * currentNode;
 if( f_treeRoot ) {
  delete rootNode;                  // the fresh root is not needed
  res = reseedFrontier( open , solvers , minimizing , rootNode , currentNode );
  }
 else {
  if( boundingProtocol == Lazy ) {
   // like any other node, the root is evaluated when it is extracted
   rootNode->initializeBound( minimizing );
   res = Solver::kOK;
   }
  else
   res = initializeRoot( &f_RelaxationSolvers , &f_HeuristicSolvers ,
                         rootNode , minimizing , bestBound , bestSolution ,
                         f_incumbentCell , branchSolver , globalMutex ,
                         relTol , absTol );
  if( res == Solver::kOK )
   open.push( rootNode );
  currentNode = rootNode;
  }

 if( res != Solver::kOK ) {
  discardTree( rootNode , open , solvers );
  return( res );
  }

 return( explore( open , globalMutex , solvers , minimizing , counter ,
                  rootNode , currentNode , retain , start ) );

 }  // end( BranchAndXSolver::treeSolve )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

int BranchAndXSolver::workerDFS( Node * currentNode ,
                                 std::list< ChangeSolver * > & solvers ,
                                 std::vector< std::pair< RelaxationSolver * ,
                                              Solver * > > & relaxation ,
                                 std::vector< std::pair< ChangeSolver * ,
                                              Solver * > > & heuristic ,
                                 bool minimizing ,
                                 std::mutex & incumbentMutex ,
                                 std::atomic< int > & nodeBdg ,
                                 std::chrono::steady_clock::time_point
                                                                  deadline ,
                                 int & nameCounter )
{
 if( --nodeBdg <= 0 )
  return( Solver::kStopIter );
 if( std::chrono::steady_clock::now() > deadline )
  return( Solver::kStopTime );

 RelaxationSolver * branchSolver = nullptr;
 bool toPrune = false;
 // out-of-mutex reads of the shared incumbent may be slightly stale, only
 // making the pruning marginally less aggressive; every incumbent update
 // happens under incumbentMutex with a double check. The solvers run serially
 // (the per-node-solvers parallelism does not lock the incumbent), so the
 // globalMutex argument is unused here
 int res = evaluateNode( currentNode , solvers , &relaxation , &heuristic ,
                         minimizing , bestBound , bestSolution ,
                         f_incumbentCell , toPrune ,
                         branchSolver , 1 ,
                         relTol , absTol , incumbentMutex , &incumbentMutex );
 if( res != ThinComputeInterface::kOK )
  return( res );

 if( ( ! toPrune ) &&
     ( ! cannot_improve( currentNode->get_bound() , f_incumbentCell->load() ,
                         minimizing , relTol , absTol ) ) ) {
  auto branches = branchSolver->branch();
  for( auto br : branches ) {
   DFSNode * new_node = new DFSNode( br , ++nameCounter );
   new_node->set_bound( currentNode->get_bound() );  // inherited
   int RV = workerDFS( new_node , solvers , relaxation , heuristic ,
                       minimizing , incumbentMutex , nodeBdg , deadline ,
                       nameCounter );
   delete new_node;
   if( RV != Solver::kOK )
    return( RV );
   }
  }

 moveSolverToFather( currentNode , &solvers );
 return( Solver::kOK );

 }  // end( BranchAndXSolver::workerDFS )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::ParallelDFSSolve( std::mutex & globalMutex , int K )
{
 createWorkerSolvers( K );

 // the serial Solver set drives the ramp-up
 std::list< ChangeSolver * > solvers;
 bool minimizing;
 initializeVariables( &solvers , &f_RelaxationSolvers , &f_HeuristicSolvers ,
                      minimizing );

 // an infinite time budget must not be duration_cast (it overflows into a
 // deadline in the past): it simply means no deadline at all
 const auto deadline = timeBudget == Inf< double >() ?
  std::chrono::steady_clock::time_point::max() :
  std::chrono::steady_clock::now() +
   std::chrono::duration_cast< std::chrono::steady_clock::duration >(
                              std::chrono::duration< double >( timeBudget ) );
 std::atomic< int > nodeBdg( nodeBudget );

 // ramp-up- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // ordered (FIFO = discovery order) BestFS-style expansion with the serial
 // Solver set: every node put in the pool has its relaxation bound, branching
 // Changes (see obtainBranchList()) and undo (toFather) already computed,
 // so the workers can later claim it by just re-applying the Changes found
 // on its path, which are plain data usable by any Solver
 int nameCounter = 0;
 ExploringNode * rootNode = new ExploringNode( nullptr , nullptr , 0 );
 RelaxationSolver * branchSolver = nullptr;
 int res = initializeRoot( &f_RelaxationSolvers , &f_HeuristicSolvers ,
                           rootNode , minimizing , bestBound , bestSolution ,
                           f_incumbentCell , branchSolver , globalMutex ,
                           relTol , absTol );
 if( res != Solver::kOK ) {
  delete rootNode;
  return( res );
  }

 std::deque< ExploringNode * > pool;
 pool.push_back( rootNode );
 ExploringNode * currentNode = rootNode;
 const std::size_t target = std::size_t( 4 * K );

 while( ( ! pool.empty() ) && ( pool.size() < target ) &&
        ( --nodeBdg > 0 ) &&
        ( std::chrono::steady_clock::now() < deadline ) ) {
  ExploringNode * oldNode = currentNode;
  currentNode = pool.front();
  pool.pop_front();
  if( currentNode->get_parent() )
   ExploringNode::moveBetweenNodes( oldNode , currentNode , &solvers );

  if( ! cannot_improve( currentNode->get_bound() , f_incumbentCell->load() ,
                        minimizing , relTol , absTol ) ) {
   auto & branches = currentNode->getBranches();
   for( Change * br : branches ) {
    auto new_node = new ExploringNode( br , currentNode ,
                                       currentNode->get_level() + 1 ,
                                       ++nameCounter );
    new_node->set_bound( currentNode->get_bound() );  // inherited
    bool toPrune = false;
    RelaxationSolver * nodeBranchSolver = nullptr;
    res = evaluateNode( new_node , solvers , &f_RelaxationSolvers ,
                        &f_HeuristicSolvers , minimizing , bestBound ,
                        bestSolution , f_incumbentCell , toPrune ,
                        nodeBranchSolver ,
                        1 , relTol , absTol , globalMutex );
    if( res != Solver::kOK ) {
     moveSolverToFather( new_node , &solvers );
     delete new_node;
     break;
     }
    // either way the ownership of br leaves the branches of currentNode
    // (the kept child owns it as its f_change, the discarded one died
    // with it): null the entry so that ~Node does not double-delete it
    *( std::find( branches.begin() , branches.end() ,
                  new_node->get_f_change() ) ) = nullptr;
    if( ( ! toPrune ) &&
        ( ! cannot_improve( new_node->get_bound() , f_incumbentCell->load() ,
                            minimizing , relTol , absTol ) ) ) {
     new_node->obtainBranchList( nodeBranchSolver );
     pool.push_back( new_node );
     currentNode->get_children().push_back( new_node );
     moveSolverToFather( new_node , &solvers );
     }
    else {                       // fenced or pruned: discard the child
     moveSolverToFather( new_node , &solvers );
     delete new_node;
     }
    }
   if( res != Solver::kOK )
    break;
   }
  }

 // back to the root, ready for the workers
 ExploringNode::moveBetweenNodes( currentNode , rootNode , &solvers );

 // workers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // each worker repeatedly claims the OLDEST open subtree of the pool
 // (preserving the sequential search order), positions its own Solver set
 // on it by re-applying the Changes of its path, explores it depth-first
 // and un-winds back to the root via the (already computed) toFather
 std::atomic< int > result( res );
 if( ( res == Solver::kOK ) && ( ! pool.empty() ) ) {
  std::mutex poolMutex , incumbentMutex;
  std::atomic< std::size_t > poolSize( pool.size() );
  std::atomic< int > busy( 0 );

  auto workerLoop = [ & ]( int w ) {
   auto & ws = v_workerSolvers[ w ];
   std::list< ChangeSolver * > wsolvers;
   for( auto & pr : ws.relaxation )
    wsolvers.push_back( pr.first );
   for( auto & pr : ws.heuristic )
    wsolvers.push_back( pr.first );
   int wNameCounter = ( w + 1 ) * 10000000;  // disjoint per-worker names

   while( result.load() == Solver::kOK ) {
    ExploringNode * node;
    {
     std::lock_guard< std::mutex > guard( poolMutex );
     if( pool.empty() ) {
      // the pool may be re-fed by the work-sharing of a busy worker: only
      // an empty pool with NO busy worker means there is no work left
      if( busy.load() == 0 )
       break;
      node = nullptr;
      }
     else {
      node = pool.front();
      pool.pop_front();
      --poolSize;
      ++busy;
      }
     }
    if( ! node ) {              // pool empty but someone is still working
     std::this_thread::sleep_for( std::chrono::microseconds( 100 ) );
     continue;
     }

    // claim the subtree: re-apply the path of Changes from the root
    std::list< Change * > path;
    for( auto n = node ; n ; n = n->get_parent() )
     if( n->get_f_change() )
      path.push_front( n->get_f_change() );
    for( auto chg : path )
     for( auto s : wsolvers )
      s->apply( chg , false );

    // explore the subtree depth-first, unless meanwhile fenced
    if( ! cannot_improve( node->get_bound() , f_incumbentCell->load() ,
                          minimizing , relTol , absTol ) ) {
     bool first = true;
     for( auto & br : node->getBranches() ) {
      if( ! br )                 // child discarded during the ramp-up
       continue;

      // work-sharing: if the pool is REALLY starving (less than half the
      // workers could claim something), donate this branch to it (eagerly
      // evaluated like in the ramp-up, so it is claimable) rather than
      // exploring it; the donated child hangs off the claimed node, which
      // always outlives it, so the path stays replayable. The first branch
      // is always explored locally, and the threshold is conservative: on
      // small subtrees an aggressive donation policy costs more (one eager
      // evaluation per donation) than it parallelizes
      if( ( ! first ) && ( 2 * poolSize.load() < std::size_t( K ) ) ) {
       auto new_node = new ExploringNode( br , node ,
                                          node->get_level() + 1 ,
                                          ++wNameCounter );
       br = nullptr;
       moveSolverToSon( new_node , &wsolvers );
       new_node->set_bound( node->get_bound() );  // inherited
       bool toPrune = false;
       RelaxationSolver * nodeBranchSolver = nullptr;
       int RV = computeRelaxations( &ws.relaxation , new_node , minimizing ,
                                    bestBound , bestSolution ,
                                    f_incumbentCell , toPrune ,
                                    nodeBranchSolver , relTol , absTol ,
                                    1 , &incumbentMutex );
       if( RV == Solver::kOK && ! toPrune )
        RV = computeHeuristic( &ws.heuristic , new_node , minimizing ,
                               bestBound , bestSolution ,
                               f_incumbentCell , toPrune ,
                               1 , &incumbentMutex );
       // the branches are asked while the Solver are still at the child,
       // whose relaxation they come out of, as everywhere else [see
       // obtainBranchList()]: only then the Solver goes back to the father
       const bool donate = ( RV == Solver::kOK ) && ( ! toPrune ) &&
                           ( ! cannot_improve( new_node->get_bound() ,
                                               f_incumbentCell->load() ,
                                               minimizing , relTol ,
                                               absTol ) );
       if( donate )
        new_node->obtainBranchList( nodeBranchSolver );
       moveSolverToFather( new_node , &wsolvers );
       if( RV != Solver::kOK ) {
        delete new_node;
        int expected = Solver::kOK;
        result.compare_exchange_strong( expected , RV );
        break;
        }
       if( donate ) {
        node->get_children().push_back( new_node );
        std::lock_guard< std::mutex > guard( poolMutex );
        pool.push_back( new_node );
        ++poolSize;
        }
       else                      // pruned or fenced: nothing to donate
        delete new_node;
       continue;
       }

      // the DFSNode takes the ownership of the branching Change: null the
      // entry so that ~Node does not double-delete it on the final cleanup
      DFSNode * new_node = new DFSNode( br , ++wNameCounter );
      new_node->set_bound( node->get_bound() );  // inherited
      br = nullptr;
      first = false;
      int RV = workerDFS( new_node , wsolvers , ws.relaxation ,
                          ws.heuristic , minimizing , incumbentMutex ,
                          nodeBdg , deadline , wNameCounter );
      delete new_node;
      if( RV != Solver::kOK ) {
       int expected = Solver::kOK;
       result.compare_exchange_strong( expected , RV );
       break;
       }
      }
     }

    // un-wind back to the root via the toFather of the path
    for( auto n = node ; n ; n = n->get_parent() )
     if( n->get_toFather() )
      for( auto s : wsolvers )
       s->apply( n->get_toFather() , false );

    --busy;
    }
   };

  std::vector< std::thread > workers;
  workers.reserve( K );
  for( int w = 0 ; w < K ; ++w )
   workers.emplace_back( workerLoop , w );
  for( auto & t : workers )
   t.join();
  }



 // cleanup of the ramp-up tree- - - - - - - - - - - - - - - - - - - - - - - -
 std::function< void( ExploringNode * ) > deleteTree =
  [ & ]( ExploringNode * node ) {
   for( auto child : node->get_children() )
    deleteTree( child );
   node->get_children().clear();
   delete node;
   };
 deleteTree( rootNode );

 if( ( result.load() == Solver::kOK ) && ( nodeBdg <= 0 ) )
  return( Solver::kStopIter );
 return( result.load() );

 }  // end( BranchAndXSolver::ParallelDFSSolve )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::PBeFSSolve( std::mutex & globalMutex , int K )
{
 K = std::max( 1 , K );
 createWorkerSolvers( Index( K ) );
 const bool minimizing =
  ( f_Block->get_objective_sense() != Objective::eMax );
 BXContext c( bestBound , bestSolution , f_incumbentCell , relTol , absTol ,
              minimizing , boundingProtocol == Eager ,
              ( solveType == DFS ) || ( solveType == BestFSDive ) ,
              treeExplorer == PBeFSDeterministic , globalMutex );
 for( auto & ws : v_workerSolvers ) {
  c.relaxation.push_back( & ws.relaxation );
  c.heuristic.push_back( & ws.heuristic );
  }

 // the strategy as a selection of PBeFS: the dive of BestFSDive goes
 // down to a leaf whatever the bound of the child
 PB::ParallelBBOptions o;
 switch( solveType ) {
  case( DFS ):    o.selection = PB::Selection::depth_first; break;
  case( BFS ):    o.selection = PB::Selection::breadth_first; break;
  case( BestFS ): o.selection = PB::Selection::best_first; break;
  case( BestEstimate ): o.selection = PB::Selection::best_estimate; break;
  default:
   o.selection = PB::Selection::plunging;
   o.plunge_gap = Inf< double >();
  }
 o.rel_gap = relTol;
 o.abs_gap = absTol == Inf< double >() ? 0 : absTol;
 if( nodeBudget < INT_MAX )
  o.node_limit = std::size_t( std::max( 1 , nodeBudget ) );
 if( timeBudget < Inf< double >() )
  o.time_limit = timeBudget;
 o.threads = unsigned( K );

 const BXTree t{ c };
 PB::BBSolution< BXTree::node_type > s;
 if( c.deterministic )
  s = PB::DeterministicBranchAndBound< BXTree >( t , o ).solve();
 else
  if( K == 1 )
   s = PB::branch_and_bound( t , o );
  else
   if( treeExplorer == PBeFSDistributed )
    s = PB::DistributedBranchAndBound< BXTree >( t , o ).solve();
   else
    s = PB::SharedBranchAndBound< BXTree >( t , o ).solve();

 if( c.error.load() != kOK )
  return( c.error.load() );
 // by synchronous rounds, the solution is that of the incumbent node
 if( c.deterministic && s.node && ( *s.node )->solution ) {
  bestBound = ( *s.node )->found;
  f_incumbentCell->store( bestBound );
  delete bestSolution;
  bestSolution = ( *s.node )->solution.release();
  }
 switch( s.status ) {
  case( PB::BBStatus::node_limit ):
   f_open_bound = c.sign() * s.bound;
   return( kStopIter );
  case( PB::BBStatus::time_limit ):
   f_open_bound = c.sign() * s.bound;
   return( kStopTime );
  default:
   return( kOK );
  }
 }  // end( BranchAndXSolver::PBeFSSolve )

/*--------------------------------------------------------------------------*/

void BranchAndXSolver::discardRetainedTree( void )
{
 if( ! f_treeRoot )
  return;

 std::function< void( ExploringNode * ) > deleteTree =
  [ & ]( ExploringNode * node ) {
   for( auto child : node->get_children() )
    deleteTree( child );
   node->get_children().clear();
   delete node;
   };
 deleteTree( f_treeRoot );

 f_treeRoot = nullptr;
 subOptimalNodes.clear();
 infeasibleNodes.clear();
 integerNodes.clear();

 }  // end( BranchAndXSolver::discardRetainedTree )

/*--------------------------------------------------------------------------*/

void BranchAndXSolver::process_outstanding_Modification( void )
{
 Lst_sp_Mod v_mod_tmp;              // temporary list of modifications

 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 for( auto mod : v_mod )
  v_mod_tmp.push_back( mod );       // copy v_mod in v_mod_tmp

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 if( v_mod_tmp.empty() )
  return;

 // what the Modification change, composed by or-ing their bits: the
 // Modification say it themselves [see Modification::changes()], whatever
 // the RelaxationSolver; one at a time, those in a GroupModification
 // included, since whether the tree survives depends on the kind and the
 // direction of each change [see keeps_tree()], and this is lost once its
 // bits are or-ed with those of others: one that the tree does not survive,
 // or a physical one that says nothing of the model, counts as anything
 std::function< void( const Modification & ) > add =
  [ & ]( const Modification & mod ) {
   if( auto gmod = dynamic_cast< const GroupModification * >( & mod ) ) {
    for( const auto & sub : gmod->sub_Modifications() )
     add( *sub );
    return;
    }
   const auto c = mod.changes();
   if( ( Modification::is_physical( c ) &&
         ( ! Modification::changes_variables( c ) ) &&
         ( ! Modification::changes_constraints( c ) ) &&
         ( ! Modification::changes_objective( c ) ) ) ||
       ( ! keeps_tree( c ) ) )
    changes = Modification::eModAll;
   else
    changes |= c;
   };

 for( const auto & mod : v_mod_tmp )
  add( *mod );

 }  // end( BranchAndXSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/

bool BranchAndXSolver::keeps_tree( Modification::ModConcern c )
{
 // the branching Changes act on the Variable and on the Constraint (fixing
 // the former, changing the bounds of either): a change of the data of the
 // Variable may be undone by them when moving across the tree, a removed
 // Variable or Constraint may be one they act upon, and a Variable that
 // becomes continuous may have been branched upon as an integer one
 // (x <= 0 and x >= 1 do not cover ( 0 , 1 )), so that the leaves no longer
 // cover the region; anything else is a relaxation or a restriction that
 // the re-evaluation of the frontier takes care of [see reseedFrontier()],
 // an added Variable being free in every node and an added Constraint
 // holding in every node
 if( c & Modification::eModVarData )
  return( false );
 if( Modification::changes_integrality( c ) &&
     Modification::may_grow_region( c ) )
  return( false );
 if( ( c & Modification::eModVarSet ) &&
     Modification::may_shrink_region( c ) )
  return( false );
 if( ( c & Modification::eModCnsSet ) &&
     Modification::may_grow_region( c ) )
  return( false );
 return( true );

 }  // end( BranchAndXSolver::keeps_tree )

/*--------------------------------------------------------------------------*/
/*------------------ End File BranchAndXSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/

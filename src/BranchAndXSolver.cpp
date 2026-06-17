/*--------------------------------------------------------------------------*/
/*--------------------- File BranchAndXSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class BranchAndXSolver, a generic
 * Branch-and-Bound Solver built on top of the ChangeSolver /
 * RelaxationSolver concepts [see ChangeSolver.h]: the attached
 * RelaxationSolver(s) provide dual bounds, true solutions and the branching
 * Changes, the attached heuristic ChangeSolver(s) further primal bounds,
 * and the enumeration tree is navigated by applying (undo) Changes to the
 * Solver. See the file-level comment of BranchAndXSolver.h for an overview.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Federica Di Pasquale \n
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
 * Copyright &copy by Antonio Frangioni, Federica Di Pasquale, Filippo
 * Magi, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <chrono>
#include <cmath>
#include <deque>
#include <stack>
#include <thread>
#include <functional>

#include "BranchAndXSolver.h"

#include "GroupChange.h"

#include "Objective.h"

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

/// tells if a dual bound cannot improve the incumbent beyond the tolerances
/** True if \p dual cannot improve \p best by more than
 * max( absAcc , relAcc * max( | best | , 1 ) ), i.e., the node can be
 * pruned within the required optimality tolerances. */

static bool cannot_improve( double dual , double best , bool minimizing ,
                            double relAcc , double absAcc )
{
 if( std::isinf( best ) )      // no incumbent yet: everything can improve
  return( false );
 // an absAcc at its default +Inf means "not active" [see Solver::dblAbsAcc]
 const double eps = std::max( absAcc == Inf< double >() ? 0.0 : absAcc ,
                              relAcc * std::max( std::abs( best ) , 1.0 ) );
 return( minimizing ? dual >= best - eps : dual <= best + eps );
 }

/// initialize the common per-solve variables
/** @param solvers list to be filled with all the :ChangeSolver to drive
 *  @param f_RelaxationSolvers the relaxation solvers to insert in solvers
 *  @param f_HeuristicSolvers the heuristic solvers to insert in solvers
 *  @param minimizing set to true if the problem is a minimization one */

static void initializeVariables(
                  std::list< ChangeSolver * > * solvers ,
                  std::vector< RelaxationSolver * > * f_RelaxationSolvers ,
                  std::vector< ChangeSolver * > * f_HeuristicSolvers ,
                  bool & minimizing )
{
 if( ! f_HeuristicSolvers->empty() )
  minimizing = f_HeuristicSolvers->front()->get_Block()
                                  ->get_objective_sense() == Objective::eMin;
 else if( ! f_RelaxationSolvers->empty() )
  minimizing = f_RelaxationSolvers->front()->get_Block()
                                  ->get_objective_sense() == Objective::eMin;
 else
  throw( std::logic_error( "BranchAndXSolver::initializeVariables: "
         "both the HeuristicSolvers and the RelaxationSolvers are empty" ) );

 solvers->insert( solvers->end() , f_RelaxationSolvers->begin() ,
                  f_RelaxationSolvers->end() );
 solvers->insert( solvers->end() , f_HeuristicSolvers->begin() ,
                  f_HeuristicSolvers->end() );
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
/// compute the relaxations at the current node
/** Computes every RelaxationSolver at the current node, updating the node
 * dual bound, the branching solver, and, when a true solution improves it,
 * the incumbent bestBound / bestSol; sets \p toPrune when the node can
 * be discarded.
 *  @return the sol_type [see Solver.h] of the computation */

static int computeRelaxations(
                  std::vector< RelaxationSolver * > * f_RelaxationSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol ,
                  bool & toPrune , RelaxationSolver * & branchSolver ,
                  double relAcc = 0 , double absAcc = 0 ,
                  std::mutex * incumbentMutex = nullptr )
{
 for( auto s : *f_RelaxationSolvers ) {
  auto z = s->compute();
  if( z == Solver::kInfeasible ) {
   if( ! currentNode->get_toFather() )    // the root is infeasible
    return( Solver::kInfeasible );
   toPrune = true;
   currentNode->set_infeasible( true );
   break;
   }
  if( z != ThinComputeInterface::kOK )
   return( z );

  // see if the primal bound improves thanks to a true solution
  if( s->has_true_var_solution() ) {
   double primal_bound = minimizing ? s->get_true_ub() : s->get_true_lb();
   if( minimizing ? primal_bound < bestBound : primal_bound > bestBound ) {
    // in the parallel exploration the incumbent is shared between the
    // workers: re-check the improvement under the mutex
    std::unique_lock< std::mutex > guard;
    if( incumbentMutex )
     guard = std::unique_lock< std::mutex >( *incumbentMutex );
    if( minimizing ? primal_bound < bestBound : primal_bound > bestBound ) {
     bestBound = primal_bound;
     delete bestSol;
     bestSol = s->get_true_solution();
     }
    }
   }

  // update the dual bound of the node and select the branching solver
  auto dualBound = minimizing ? s->get_lb() : s->get_ub();
  if( cannot_improve( dualBound , bestBound , minimizing , relAcc ,
                      absAcc ) ) {
   toPrune = true;
   return( Solver::kOK );
   }
  if( minimizing ? dualBound > currentNode->get_dual_bound()
                 : dualBound < currentNode->get_dual_bound() ) {
   currentNode->set_dual_bound( dualBound );
   branchSolver = s;
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
                  std::vector< ChangeSolver * > * f_HeuristicSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol , bool & toPrune ,
                  std::mutex * incumbentMutex = nullptr )
{
 for( auto s : *f_HeuristicSolvers ) {
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
      ( minimizing ? primal_bound < bestBound
                   : primal_bound > bestBound ) ) {
   std::unique_lock< std::mutex > guard;
   if( incumbentMutex )
    guard = std::unique_lock< std::mutex >( *incumbentMutex );
   if( minimizing ? primal_bound < bestBound : primal_bound > bestBound ) {
    bestBound = primal_bound;
    delete bestSol;
    bestSol = s->get_Solution();
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
   zs[ i ] = slvrs[ i ]->compute();
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
/// parallel version of computeRelaxations()
/** Computes every RelaxationSolver at the current node in parallel (only the
 * compute() is parallel, see computeAllParallel()), then reduces the results
 * with EXACTLY the serial logic of computeRelaxations(), same order of the
 * incumbent updates, same dual bound, same pruning, so that the outcome is
 * bit-identical to the serial path, only faster when several (expensive)
 * relaxations are attached. */

static int computeRelaxationsParallel(
                  std::vector< RelaxationSolver * > * f_RelaxationSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol , bool & toPrune ,
                  RelaxationSolver * & branchSolver ,
                  double relAcc , double absAcc ,
                  std::mutex & globalMutex , int maxThreads )
{
 const std::size_t n = f_RelaxationSolvers->size();
 if( n <= 1 )                  // nothing to parallelize
  return( computeRelaxations( f_RelaxationSolvers , currentNode , minimizing ,
                              bestBound , bestSol , toPrune , branchSolver ,
                              relAcc , absAcc , nullptr ) );

 std::vector< int > zs( n , INT_MIN );      // INT_MIN: not computed (early stop)
 if( computeAllParallel( *f_RelaxationSolvers , zs , maxThreads ) ) {
  // one relaxation proved the node infeasible: prune it, the others were
  // skipped and their entries of zs are not to be read
  if( ! currentNode->get_toFather() )       // the root is infeasible
   return( Solver::kInfeasible );
  toPrune = true;
  currentNode->set_infeasible( true );
  return( Solver::kOK );
  }

 // no infeasibility: every relaxation was computed, reduce with EXACTLY the
 // serial logic of computeRelaxations() but reading the pre-computed codes
 for( std::size_t i = 0 ; i < n ; ++i ) {
  auto s = (*f_RelaxationSolvers)[ i ];
  const auto z = zs[ i ];
  if( z != ThinComputeInterface::kOK )
   return( z );

  if( s->has_true_var_solution() ) {
   double primal_bound = minimizing ? s->get_true_ub() : s->get_true_lb();
   if( minimizing ? primal_bound < bestBound : primal_bound > bestBound ) {
    bestBound = primal_bound;
    delete bestSol;
    bestSol = s->get_true_solution();
    }
   }

  auto dualBound = minimizing ? s->get_lb() : s->get_ub();
  if( cannot_improve( dualBound , bestBound , minimizing , relAcc ,
                      absAcc ) ) {
   toPrune = true;
   return( Solver::kOK );
   }
  if( minimizing ? dualBound > currentNode->get_dual_bound()
                 : dualBound < currentNode->get_dual_bound() ) {
   currentNode->set_dual_bound( dualBound );
   branchSolver = s;
   }
  }
 return( Solver::kOK );
 }

/*--------------------------------------------------------------------------*/
/// parallel version of computeHeuristic()
/** Computes every heuristic ChangeSolver at the current node in parallel
 * (only the compute() is parallel), then reduces with the serial logic of
 * computeHeuristic(), bit-identical to the serial path. */

static int computeHeuristicParallel(
                  std::vector< ChangeSolver * > * f_HeuristicSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol , bool & toPrune ,
                  std::mutex & globalMutex , int maxThreads )
{
 const std::size_t n = f_HeuristicSolvers->size();
 if( n <= 1 )                  // nothing to parallelize
  return( computeHeuristic( f_HeuristicSolvers , currentNode , minimizing ,
                            bestBound , bestSol , toPrune ) );

 std::vector< int > zs( n , INT_MIN );      // INT_MIN: not computed (early stop)
 if( computeAllParallel( *f_HeuristicSolvers , zs , maxThreads ) ) {
  // one heuristic proved the node infeasible: prune it, the others were
  // skipped and their entries of zs are not to be read
  if( ! currentNode->get_toFather() )       // the root is infeasible
   return( Solver::kInfeasible );
  toPrune = true;
  return( Solver::kOK );
  }

 for( std::size_t i = 0 ; i < n ; ++i ) {
  auto s = (*f_HeuristicSolvers)[ i ];
  const auto z = zs[ i ];
  if( z != ThinComputeInterface::kOK )
   return( z );

  double primal_bound = minimizing ? s->get_ub() : s->get_lb();
  if( s->has_var_solution() && s->is_var_feasible() &&
      ( minimizing ? primal_bound < bestBound
                   : primal_bound > bestBound ) ) {
   bestBound = primal_bound;
   delete bestSol;
   bestSol = s->get_Solution();
   }
  }
 return( Solver::kOK );
 }

/*--------------------------------------------------------------------------*/
/// evaluate the node relaxations, serially or in parallel per \p nThreads
/** Thin dispatcher used by every (serial-tree) exploration: with
 * \p nThreads <= 1 it calls computeRelaxations(), otherwise
 * computeRelaxationsParallel(); both give the same result [see the latter]. */

static int evalRelaxations(
                  int nThreads ,
                  std::vector< RelaxationSolver * > * f_RelaxationSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol , bool & toPrune ,
                  RelaxationSolver * & branchSolver , double relAcc ,
                  double absAcc , std::mutex & globalMutex )
{
 if( nThreads <= 1 )
  return( computeRelaxations( f_RelaxationSolvers , currentNode , minimizing ,
                              bestBound , bestSol , toPrune , branchSolver ,
                              relAcc , absAcc , nullptr ) );
 return( computeRelaxationsParallel( f_RelaxationSolvers , currentNode ,
                                     minimizing , bestBound , bestSol ,
                                     toPrune , branchSolver , relAcc , absAcc ,
                                     globalMutex , nThreads ) );
 }

/*--------------------------------------------------------------------------*/
/// evaluate the node heuristics, serially or in parallel per \p nThreads

static int evalHeuristics(
                  int nThreads ,
                  std::vector< ChangeSolver * > * f_HeuristicSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol , bool & toPrune ,
                  std::mutex & globalMutex )
{
 if( nThreads <= 1 )
  return( computeHeuristic( f_HeuristicSolvers , currentNode , minimizing ,
                            bestBound , bestSol , toPrune ) );
 return( computeHeuristicParallel( f_HeuristicSolvers , currentNode ,
                                   minimizing , bestBound , bestSol , toPrune ,
                                   globalMutex , nThreads ) );
 }


/*--------------------------------------------------------------------------*/
/// evaluate the root node and produce its branching list

static int initializeRoot(
                  std::vector< RelaxationSolver * > * f_RelaxationSolvers ,
                  std::vector< ChangeSolver * > * f_HeuristicSolvers ,
                  Node * rootNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol ,
                  RelaxationSolver * & branchSolver ,
                  std::mutex & globalMutex )
{
 rootNode->initializeBound( minimizing );
 bool toPrune = false;
 int res = computeRelaxations( f_RelaxationSolvers , rootNode , minimizing ,
                               bestBound , bestSol , toPrune , branchSolver );
 if( toPrune || ( res != ThinComputeInterface::kOK ) )
  return( res );
 res = computeHeuristic( f_HeuristicSolvers , rootNode , minimizing ,
                         bestBound , bestSol , toPrune );
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
                  std::vector< RelaxationSolver * > * relaxation ,
                  std::vector< ChangeSolver * > * heuristic ,
                  bool minimizing , double & bestBound , Solution * & bestSol ,
                  bool & toPrune , RelaxationSolver * & branchSolver ,
                  int nThreads ,
                  double relAcc , double absAcc , std::mutex & globalMutex ,
                  std::mutex * incumbentMutex = nullptr )
{
 moveSolverToSon( node , &solvers );
 node->initializeBound( minimizing );

 // the per-node-solvers parallelism does not lock the incumbent, so a worker
 // of the parallel tree exploration (incumbentMutex set) runs the solvers
 // serially; otherwise nThreads decides
 int res;
 if( incumbentMutex || ( nThreads <= 1 ) )
  res = computeRelaxations( relaxation , node , minimizing , bestBound ,
                            bestSol , toPrune , branchSolver , relAcc ,
                            absAcc , incumbentMutex );
 else
  res = computeRelaxationsParallel( relaxation , node , minimizing , bestBound ,
                                    bestSol , toPrune , branchSolver , relAcc ,
                                    absAcc , globalMutex , nThreads );
 if( ( res != Solver::kOK ) || toPrune )
  return( res );

 if( incumbentMutex || ( nThreads <= 1 ) )
  res = computeHeuristic( heuristic , node , minimizing , bestBound , bestSol ,
                          toPrune , incumbentMutex );
 else
  res = computeHeuristicParallel( heuristic , node , minimizing , bestBound ,
                                  bestSol , toPrune , globalMutex , nThreads );
 return( res );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OPEN-LIST DISCIPLINES -------------------------*/
/*--------------------------------------------------------------------------*/
// the three open-list disciplines, as thin adapters that give the standard
// container adapters a common interface [see OpenList]: a LIFO stack (depth-
// first), a FIFO queue (breadth-first), a dual-bound priority queue (best-
// first). The storage is entirely std::stack / std::queue / std::priority_-
// queue, there is no hand-rolled data structure here

namespace {

class StackOpenList final : public OpenList {
 std::stack< ExploringNode * > c;
 public:
 [[nodiscard]] bool empty( void ) const override { return( c.empty() ); }
 void push( ExploringNode * node ) override { c.push( node ); }
 ExploringNode * pop( void ) override
  { auto node = c.top(); c.pop(); return( node ); }
 [[nodiscard]] bool isLIFO( void ) const override { return( true ); }
 };

class QueueOpenList final : public OpenList {
 std::queue< ExploringNode * > c;
 public:
 [[nodiscard]] bool empty( void ) const override { return( c.empty() ); }
 void push( ExploringNode * node ) override { c.push( node ); }
 ExploringNode * pop( void ) override
  { auto node = c.front(); c.pop(); return( node ); }
 };

class PriorityOpenList final : public OpenList {
 using Cmp = std::function< bool( ExploringNode * , ExploringNode * ) >;
 std::priority_queue< ExploringNode * , std::vector< ExploringNode * > , Cmp >
  c;
 public:
 explicit PriorityOpenList( Cmp cmp ) : c( std::move( cmp ) ) {}
 [[nodiscard]] bool empty( void ) const override { return( c.empty() ); }
 void push( ExploringNode * node ) override { c.push( node ); }
 ExploringNode * pop( void ) override
  { auto node = c.top(); c.pop(); return( node ); }
 };

}  // anonymous namespace

/*--------------------------------------------------------------------------*/
/*----------------- METHODS OF BranchAndXSolver ------------------------*/
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
 std::mutex globalMutex;

 // per-solve budgets and tolerances from the inherited standard parameters
 nodeBudget = get_int_par( intMaxNodes );
 timeBudget = get_dbl_par( dblMaxTime );
 relTol = get_dbl_par( dblRelAcc );
 absTol = get_dbl_par( dblAbsAcc );

 // incumbent-dependent local fixing folded into the branching Changes is
 // unsafe when the tree is retained across re-solves with different
 // incumbents: forbid it in that case, allow it otherwise [see
 // GlobalInformation::local_fixing_allowed()]
 f_globalInfo.set_local_fixing_allowed(
                       ! ( ( reoptimize > 0 ) && ( solveType == BestFS ) ) );

 if( changes == 0 )         // nothing changed since the last solve
  f_state = old_state;
 else {
  // the outstanding changes are now properly classified (see
  // process_outstanding_Modification(): 1 = objective only, 2 = feasible
  // region only, 3 = both, 4 = anything); classes 1 / 2 / 3 reoptimize out
  // of the fenced frontier of the previous tree (e.g., an objective-only
  // change leaves every infeasibility certificate valid) rather than solving
  // from scratch, but the tree retained for reoptimization (if any) can only
  // be reused by a BestFS re-solve under class 1-3 changes when retention is
  // enabled [see intReoptimize / BestFirstSolve()]; otherwise it is discarded
  if( ( changes == 4 ) || ( solveType != BestFS ) || ( ! reoptimize ) )
   discardRetainedTree();

  bestBound = ( f_Block->get_objective_sense() == Objective::eMax )
              ? - Inf< double >() : Inf< double >();
  switch( solveType ) {
   case( DFS ): {
    if( const int K = get_int_par( intMaxThread ) ; K > 1 )
     f_state = ParallelDFSSolve( globalMutex , K );
    else
     f_state = DFSSolve( globalMutex );
    break;
    }
   case( BFS ):
    f_state = BFSSolve( globalMutex );
    break;
   case( BestFS ):
    f_state = BestFirstSolve( globalMutex );
    break;
   default:
    throw( std::invalid_argument( "BranchAndXSolver::compute: invalid "
                                  "intSolveMethod" ) );
   }
  }

 changes = 0;
 unlock();
 return( f_state );

 }  // end( BranchAndXSolver::compute )


/*--------------------------------------------------------------------------*/

void BranchAndXSolver::discardTree( ExploringNode * root , OpenList & open ,
                                    std::list< ChangeSolver * > * solvers )
{
 // the nodes still in the open set are also children in the tree, so the
 // recursive deletion of the tree covers them; emptying the open set first
 // (without deleting) avoids a double delete
 while( ! open.empty() )
  open.pop();
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

 while( ( ! open.empty() ) && ( nodeBudget > 1 ) &&
        ( timeBudget > std::chrono::duration< double >(
           std::chrono::high_resolution_clock::now() - start ).count() ) ) {
  nodeBudget--;
  oldNode = currentNode;
  currentNode = open.pop();
  if( currentNode->get_parent() )   // the root: the solvers are already there
   ExploringNode::moveBetweenNodes( oldNode , currentNode , solvers );

  // branching and evaluation of the new children
  if( ! cannot_improve( currentNode->get_dual_bound() , bestBound ,
                        minimizing , relTol , absTol ) ) {
   auto & branches = currentNode->getBranches();
   std::vector< ExploringNode * > kept;  // survivors, pushed after the loop
   for( auto br : branches ) {
    ExploringNode * new_node = new ExploringNode( br , currentNode ,
                                                  currentNode->get_level()
                                                  + 1 , ++counter );
    bool toPrune = false;
    res = evaluateNode( new_node , *solvers , &f_RelaxationSolvers ,
                        &f_HeuristicSolvers , minimizing , bestBound ,
                        bestSolution , toPrune , branchSolver ,
                        maxThreadForSolvers , relTol , absTol , globalMutex );
    if( res != Solver::kOK ) {
     discardTree( rootNode , open , solvers );
     return( res );
     }
    // keep the node only if its dual bound can improve the incumbent
    if( ( ! toPrune ) &&
        ( ! cannot_improve( new_node->get_dual_bound() , bestBound ,
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
       infeasibleNodes.push_back( new_node );  //  objective-only changes
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
   if( currentNode->get_children().empty() && currentNode->get_toFather() ) {
    if( retain ) {               // fathomed leaf: fenced frontier
     subOptimalNodes.push_back( currentNode );
     for( const auto s : *solvers )
      s->apply( currentNode->get_toFather() , false );
     currentNode = currentNode->get_parent();
     }
    else
     currentNode = ExploringNode::prune( currentNode , solvers );
    }
   }
  else if( currentNode->get_toFather() ) {
   if( retain ) {                // fenced at pop: fenced frontier
    subOptimalNodes.push_back( currentNode );
    for( const auto s : *solvers )
     s->apply( currentNode->get_toFather() , false );
    currentNode = currentNode->get_parent();
    }
   else
    currentNode = ExploringNode::prune( currentNode , solvers );
   }
  }

 // move the :ChangeSolver back to the root
 while( currentNode->get_toFather() ) {
  for( const auto s : *solvers )
   s->apply( currentNode->get_toFather() , false );
  currentNode = currentNode->get_parent();
  }

 if( retain ) {
  // retain the tree for future reoptimizations: the nodes still in the open
  // set (early stops) are open work, hence part of the frontier to re-seed
  while( ! open.empty() )
   subOptimalNodes.push_back( open.pop() );
  f_treeRoot = rootNode;
  }
 else {
  // the nodes still in the open set are also children in the tree: the
  // recursive deletion of the tree covers them (deleting them from the open
  // set too would be a double delete)
  while( ! open.empty() )
   open.pop();
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

 if( nodeBudget <= 0 )
  return( Solver::kStopIter );
 if( timeBudget <= std::chrono::duration< double >(
                 std::chrono::high_resolution_clock::now() - start ).count()
     )
  return( Solver::kStopTime );
 return( kOK );

 }  // end( BranchAndXSolver::explore )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::BestFirstSolve( std::mutex & globalMutex )
{
 auto * solvers = new std::list< ChangeSolver * >();
 bool minimizing;
 initializeVariables( solvers , &f_RelaxationSolvers , &f_HeuristicSolvers ,
                      minimizing );

 auto start = std::chrono::high_resolution_clock::now();
 int counter = 0;
 // the open set is a priority queue, best dual bound first
 auto cmp = [ minimizing ]( ExploringNode * a , ExploringNode * b ) {
  return( minimizing ? a->get_dual_bound() > b->get_dual_bound()
                     : a->get_dual_bound() < b->get_dual_bound() );
  };
 PriorityOpenList open( cmp );
 ExploringNode * rootNode = new ExploringNode( nullptr , nullptr , 0 );
 RelaxationSolver * branchSolver = nullptr;

 const bool retain = ( reoptimize > 0 );
 int res;
 ExploringNode * currentNode;
 if( f_treeRoot ) {
  // reoptimization: re-seed from the fenced frontier of the previous tree
  // (compute() guarantees the retained tree is only reused for class 1-3
  // changes under BestFS): the interior of the tree is NOT re-derived;
  // every frontier node is re-evaluated under the new data and either
  // re-fenced or re-opened into the open set. Class-specific saving: an
  // objective-only change cannot un-fence an infeasible node, so the
  // infeasible part of the frontier is left untouched (see below)
  delete rootNode;                  // the fresh root is not needed
  rootNode = f_treeRoot;
  f_treeRoot = nullptr;             // ownership back to this solve
  currentNode = rootNode;
  std::list< ExploringNode * > frontier;
  frontier.swap( subOptimalNodes );
  if( changes == RelaxationSolver::eModObjective ) {
   // an objective-only change cannot un-fence an infeasible node: that
   // part of the frontier stays fenced, with no re-evaluation at all
   }
  else {
   frontier.splice( frontier.end() , infeasibleNodes );
   infeasibleNodes.clear();
   }
  // best (previous) bound first: the optimum almost surely lives in the
  // first few nodes, so the incumbent warms up immediately and the rest of
  // the frontier mostly just re-fences
  frontier.sort( [ minimizing ]( ExploringNode * a , ExploringNode * b ) {
   return( minimizing ? a->get_dual_bound() < b->get_dual_bound()
                      : a->get_dual_bound() > b->get_dual_bound() );
   } );
  res = Solver::kOK;
  for( auto F : frontier ) {
   ExploringNode::moveBetweenNodes( currentNode , F , solvers );
   currentNode = F;
   F->initializeBound( minimizing );
   bool toPrune = false;
   RelaxationSolver * nodeBranchSolver = nullptr;
   res = evalRelaxations( maxThreadForSolvers , &f_RelaxationSolvers , F ,
                          minimizing , bestBound , bestSolution , toPrune ,
                          nodeBranchSolver , relTol , absTol , globalMutex );
   if( ( res == Solver::kOK ) && ( ! toPrune ) )
    res = evalHeuristics( maxThreadForSolvers , &f_HeuristicSolvers , F ,
                          minimizing , bestBound , bestSolution , toPrune ,
                          globalMutex );
   if( res != Solver::kOK )
    break;
   if( ( ! toPrune ) &&
       ( ! cannot_improve( F->get_dual_bound() , bestBound , minimizing ,
                           relTol , absTol ) ) ) {
    // re-opened: the branching Changes of the previous solve, if any, are
    // still valid (any branching is), so they are reused rather than leaked
    if( F->getBranches().empty() )
     F->obtainBranchList( nodeBranchSolver );
    open.push( F );
    }
   else if( F->is_infeasible() )     // fenced again, by reason
    infeasibleNodes.push_back( F );
   else
    subOptimalNodes.push_back( F );
   }
  if( res != Solver::kOK ) {
   ExploringNode::moveBetweenNodes( currentNode , rootNode , solvers );
   discardTree( rootNode , open , solvers );
   return( res );
   }
  }
 else {
  res = initializeRoot( &f_RelaxationSolvers , &f_HeuristicSolvers ,
                        rootNode , minimizing , bestBound , bestSolution ,
                        branchSolver , globalMutex );
  if( res != Solver::kOK ) {
   discardTree( rootNode , open , solvers );
   return( res );
   }
  open.push( rootNode );
  currentNode = rootNode;
  }

 return( explore( open , globalMutex , solvers , minimizing , counter ,
                  rootNode , currentNode , retain , start ) );

 }  // end( BranchAndXSolver::BestFirstSolve )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::BFSSolve( std::mutex & globalMutex )
{
 auto * solvers = new std::list< ChangeSolver * >();
 bool minimizing;
 initializeVariables( solvers , &f_RelaxationSolvers , &f_HeuristicSolvers ,
                      minimizing );
 auto start = std::chrono::high_resolution_clock::now();
 int counter = 0;
 // the open set is a plain FIFO queue
 QueueOpenList open;
 ExploringNode * rootNode = new ExploringNode( nullptr , nullptr , 0 );
 RelaxationSolver * branchSolver = nullptr;
 int res = initializeRoot( &f_RelaxationSolvers , &f_HeuristicSolvers ,
                           rootNode , minimizing , bestBound , bestSolution ,
                           branchSolver , globalMutex );
 if( res != Solver::kOK ) {
  discardTree( rootNode , open , solvers );
  return( res );
  }
 open.push( rootNode );

 return( explore( open , globalMutex , solvers , minimizing , counter ,
                  rootNode , rootNode , false , start ) );

 }  // end( BranchAndXSolver::BFSSolve )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::DFSSolve( std::mutex & globalMutex )
{
 auto * solvers = new std::list< ChangeSolver * >();
 bool minimizing;
 initializeVariables( solvers , &f_RelaxationSolvers , &f_HeuristicSolvers ,
                      minimizing );
 auto start = std::chrono::high_resolution_clock::now();
 int counter = 0;
 // the open set is a LIFO stack: same explore() loop, depth-first order
 StackOpenList open;
 ExploringNode * rootNode = new ExploringNode( nullptr , nullptr , 0 );
 RelaxationSolver * branchSolver = nullptr;
 int res = initializeRoot( &f_RelaxationSolvers , &f_HeuristicSolvers ,
                           rootNode , minimizing , bestBound , bestSolution ,
                           branchSolver , globalMutex );
 if( res != Solver::kOK ) {
  discardTree( rootNode , open , solvers );
  return( res );
  }
 open.push( rootNode );

 return( explore( open , globalMutex , solvers , minimizing , counter ,
                  rootNode , rootNode , false , start ) );

 }  // end( BranchAndXSolver::DFSSolve )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

int BranchAndXSolver::workerDFS( Node * currentNode ,
                                 std::list< ChangeSolver * > & solvers ,
                                 std::vector< RelaxationSolver * > &
                                                                 relaxation ,
                                 std::vector< ChangeSolver * > & heuristic ,
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
 // out-of-mutex reads of the shared bestBound may be slightly stale, only
 // making the pruning marginally less aggressive; every incumbent update
 // happens under incumbentMutex with a double check. The solvers run serially
 // (the per-node-solvers parallelism does not lock the incumbent), so the
 // globalMutex argument is unused here
 int res = evaluateNode( currentNode , solvers , &relaxation , &heuristic ,
                         minimizing , bestBound , bestSolution , toPrune ,
                         branchSolver , 1 ,
                         relTol , absTol , incumbentMutex , &incumbentMutex );
 if( res != ThinComputeInterface::kOK )
  return( res );

 if( ( ! toPrune ) &&
     ( ! cannot_improve( currentNode->get_dual_bound() , bestBound ,
                         minimizing , relTol , absTol ) ) ) {
  auto branches = branchSolver->branch();
  for( auto br : branches ) {
   DFSNode * new_node = new DFSNode( br , ++nameCounter );
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
 // Solver set: every node put in the pool has its dual bound, branching
 // Changes (see obtainBranchList()) and undo (toFather) already computed,
 // so the workers can later claim it by just re-applying the Changes found
 // on its path, which are plain data usable by any Solver
 int nameCounter = 0;
 ExploringNode * rootNode = new ExploringNode( nullptr , nullptr , 0 );
 RelaxationSolver * branchSolver = nullptr;
 int res = initializeRoot( &f_RelaxationSolvers , &f_HeuristicSolvers ,
                           rootNode , minimizing , bestBound , bestSolution ,
                           branchSolver , globalMutex );
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

  if( ! cannot_improve( currentNode->get_dual_bound() , bestBound ,
                        minimizing , relTol , absTol ) ) {
   auto & branches = currentNode->getBranches();
   for( Change * br : branches ) {
    auto new_node = new ExploringNode( br , currentNode ,
                                       currentNode->get_level() + 1 ,
                                       ++nameCounter );
    bool toPrune = false;
    RelaxationSolver * nodeBranchSolver = nullptr;
    res = evaluateNode( new_node , solvers , &f_RelaxationSolvers ,
                        &f_HeuristicSolvers , minimizing , bestBound ,
                        bestSolution , toPrune , nodeBranchSolver ,
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
        ( ! cannot_improve( new_node->get_dual_bound() , bestBound ,
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
   std::list< ChangeSolver * > wsolvers( ws.relaxation.begin() ,
                                         ws.relaxation.end() );
   wsolvers.insert( wsolvers.end() , ws.heuristic.begin() ,
                    ws.heuristic.end() );
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
    if( ! cannot_improve( node->get_dual_bound() , bestBound , minimizing ,
                          relTol , absTol ) ) {
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
       new_node->initializeBound( minimizing );
       bool toPrune = false;
       RelaxationSolver * nodeBranchSolver = nullptr;
       int RV = computeRelaxations( &ws.relaxation , new_node , minimizing ,
                                    bestBound , bestSolution , toPrune ,
                                    nodeBranchSolver , relTol , absTol ,
                                    &incumbentMutex );
       if( RV == Solver::kOK && ! toPrune )
        RV = computeHeuristic( &ws.heuristic , new_node , minimizing ,
                               bestBound , bestSolution , toPrune ,
                               &incumbentMutex );
       moveSolverToFather( new_node , &wsolvers );
       if( RV != Solver::kOK ) {
        delete new_node;
        int expected = Solver::kOK;
        result.compare_exchange_strong( expected , RV );
        break;
        }
       if( ( ! toPrune ) &&
           ( ! cannot_improve( new_node->get_dual_bound() , bestBound ,
                               minimizing , relTol , absTol ) ) ) {
        new_node->obtainBranchList( nodeBranchSolver );
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

 // classify the changes by asking a RelaxationSolver: what a Modification
 // does to the fencing certificates of the tree is problem-specific
 // knowledge [see RelaxationSolver::classify()], so this Solver never
 // looks into the Modification itself; the classes compose bitwise, and
 // anything the RelaxationSolver cannot vouch for invalidates everything
 if( f_RelaxationSolvers.empty() ) {
  changes = 4;
  return;
  }
 auto rs = f_RelaxationSolvers.front();
 int acc = changes == 4 ? RelaxationSolver::eModEverything : changes;
 for( const auto & mod : v_mod_tmp ) {
  const int c = rs->classify( mod );
  if( c == RelaxationSolver::eModEverything ) {
   acc = c;
   break;
   }
  acc |= c;
  }
 changes = char( acc );

 }  // end( BranchAndXSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/
/*------------------ End File BranchAndXSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/

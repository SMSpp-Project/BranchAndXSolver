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
#include <numeric>
#include <set>
#include <stack>
#include <thread>
#include <functional>
#include <condition_variable>

#include "BranchAndXSolver.h"

#include "Change.h"

#include "Objective.h"

#include "PBeFS/bnb/BranchAndBoundScheme.h"

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

/// whether the evaluation of a node has ended inconclusively
/** True for the return codes with which a relaxation says that it has
 * neither solved nor failed (stopped by a budget of its own, or unable to
 * reach the required accuracy): the relaxation bound that the node has is
 * still valid, but the node can be neither branched upon nor fenced, and it
 * is set aside [see BranchAndXSolver::v_unresolved]. */

static bool inconclusive( int res )
{
 return( ( res == Solver::kStopTime ) || ( res == Solver::kStopIter ) ||
         ( res == Solver::kLowPrecision ) );
 }

/// the time by which a solve has to end, max() if none
using Deadline = std::chrono::steady_clock::time_point;

/// run the compute() of \p s within the time left before \p deadline
/** The dblMaxTime of \p s is the smallest between its own and the time
 * left, and it gets its own back after; with no time left \p s is not run,
 * \p ran is set to false and the result is kStopTime. */

static int compute_by( Solver * s , const Deadline deadline , bool & ran )
{
 ran = true;
 if( deadline == Deadline::max() )
  return( s->compute() );
 const double left = std::chrono::duration< double >(
			    deadline - std::chrono::steady_clock::now() ).count();
 if( left <= 0 ) {
  ran = false;
  return( Solver::kStopTime );
  }
 const double own = s->get_dbl_par( Solver::dblMaxTime );
 if( own <= left )
  return( s->compute() );
 s->set_par( Solver::dblMaxTime , left );
 int z;
 try {
  z = s->compute();
  }
 catch( ... ) {
  s->set_par( Solver::dblMaxTime , own );
  throw;
  }
 s->set_par( Solver::dblMaxTime , own );
 return( z );
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
                                std::vector< int > & zs ,
                                std::vector< char > & ran , int maxThreads ,
                                Deadline deadline );

/*--------------------------------------------------------------------------*/
/// improve the bound of the current node with that of the relaxation \p s
/** Also for a relaxation that has not solved the node, whose bound is valid
 * whatever the status of its compute() [see Solver::get_lb()]: a node set
 * aside carries it [see BranchAndXSolver::v_unresolved]. */

static void improve_bound( Node * currentNode , Solver * s ,
                           const bool minimizing )
{
 const double bound = minimizing ? s->get_lb() : s->get_ub();
 if( minimizing ? bound > currentNode->get_bound()
                : bound < currentNode->get_bound() )
  currentNode->set_bound( bound );
 }

/*--------------------------------------------------------------------------*/
/// fence the current node whose relaxation stopped at the cutoff
/** The relaxation \p s has returned kCutOff, i.e., it has proven its bound
 * beyond the cutoff it was given, the incumbent [see computeRelaxations()]:
 * the node is fenced by bound, with the bound the relaxation reports, if
 * better than the one it has. With no incumbent the cutoff given is
 * infinite, and so a kCutOff comes from a cutoff the relaxation has from
 * elsewhere, which proves nothing about the incumbent: the node keeps the
 * bound and is not fenced, kLowPrecision saying that it is neither solved
 * nor fenced [see inconclusive()].
 *  @return kOK if the node is fenced, kLowPrecision otherwise */

static int fence_by_cutoff( Node * currentNode , Solver * s ,
                            const bool minimizing , const bool has_incumbent ,
                            bool & toPrune )
{
 improve_bound( currentNode , s , minimizing );
 if( ! has_incumbent )
  return( Solver::kLowPrecision );
 toPrune = true;
 return( Solver::kOK );
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
                  const Deadline deadline ,
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

 // parallel evaluation of the (several) relaxations of this node: only when
 // more than one thread is asked, there are at least two solvers, and no
 // shared incumbent lock is in force (a parallel-tree worker evaluates its
 // per-node solvers serially); the reduction is bit-identical to the serial
 // path below
 if( ( nThreads > 1 ) && ( f_RelaxationSolvers->size() > 1 ) &&
     ( ! incumbentMutex ) ) {
  const std::size_t n = f_RelaxationSolvers->size();
  std::vector< int > zs( n , INT_MIN );   // INT_MIN: not computed (early stop)
  std::vector< char > ran( n , 0 );
  if( computeAllParallel( *f_RelaxationSolvers , zs , ran , nThreads ,
                          deadline ) ) {
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
   if( zs[ i ] == Solver::kCutOff )  // beyond the incumbent [see above]
    return( fence_by_cutoff( currentNode , s , minimizing , has_incumbent ,
                             toPrune ) );
   if( zs[ i ] != ThinComputeInterface::kOK ) {
    if( ran[ i ] && inconclusive( zs[ i ] ) )
     improve_bound( currentNode , s , minimizing );
    return( zs[ i ] );
    }
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
  bool ran;
  auto z = compute_by( s , deadline , ran );
  if( z == Solver::kInfeasible ) {
   if( ( ! currentNode->get_toFather() ) && ( ! has_incumbent ) )
    return( Solver::kInfeasible );       // the root is infeasible
   toPrune = true;
   currentNode->set_infeasible( true );
   currentNode->set_bound( minimizing ? Inf< double >()
                                           : - Inf< double >() );
   break;
   }
  if( z == Solver::kCutOff )  // beyond the incumbent [see above]
   return( fence_by_cutoff( currentNode , s , minimizing , has_incumbent ,
                            toPrune ) );
  if( z != ThinComputeInterface::kOK ) {
   if( ran && inconclusive( z ) )
    improve_bound( currentNode , s , minimizing );
   return( z );
   }

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
                  std::atomic< double > * incumbentCell ,
                  const Deadline deadline , bool & toPrune ,
                  int nThreads = 1 , std::mutex * incumbentMutex = nullptr )
{
 // parallel evaluation of the (several) heuristics of this node, bit-identical
 // to the serial path below [see computeRelaxations()]
 if( ( nThreads > 1 ) && ( f_HeuristicSolvers->size() > 1 ) &&
     ( ! incumbentMutex ) ) {
  const std::size_t n = f_HeuristicSolvers->size();
  std::vector< int > zs( n , INT_MIN );   // INT_MIN: not computed (early stop)
  std::vector< char > ran( n , 0 );
  if( computeAllParallel( *f_HeuristicSolvers , zs , ran , nThreads ,
                          deadline ) ) {
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
  bool ran;
  auto z = compute_by( s , deadline , ran );
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
 * to) \p maxThreads threads (the calling thread is one of them), within the
 * time left before \p deadline [see compute_by()], storing the return codes
 * in \p zs and in \p ran whether each was run. Only compute() is run here,
 * never any result extraction (get_lb() / get_true_solution() / ...): those
 * touch the shared Block and are left to the serial reduction in the
 * caller. This is safe
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
                                std::vector< int > & zs ,
                                std::vector< char > & ran , int maxThreads ,
                                const Deadline deadline )
{
 const int n = int( slvrs.size() );
 const int K = std::max( 1 , std::min( maxThreads , n ) );
 std::atomic< int > next( 0 );
 std::atomic< bool > infeasible( false );
 auto work = [ & ]() {
  for( int i ; ( ! infeasible.load( std::memory_order_relaxed ) ) &&
               ( ( i = next.fetch_add( 1 ) ) < n ) ; ) {
   bool r;
   zs[ i ] = compute_by( slvrs[ i ].second , deadline , r );
   ran[ i ] = r;
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
 * would be branched, with nothing left to branch upon. The root starts
 * from the bound \p rootBound [see BranchAndXSolver::f_probeBound]. */

static int initializeRoot(
                  std::vector< std::pair< RelaxationSolver * ,
                               Solver * > > * f_RelaxationSolvers ,
                  std::vector< std::pair< ChangeSolver * ,
                               Solver * > > * f_HeuristicSolvers ,
                  Node * rootNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol ,
                  std::atomic< double > * incumbentCell ,
                  const Deadline deadline ,
                  RelaxationSolver * & branchSolver ,
                  std::mutex & globalMutex ,
                  double relAcc , double absAcc , double rootBound )
{
 rootNode->initializeBound( minimizing );
 rootNode->set_bound( rootBound );
 rootNode->set_evaluated( true );
 bool toPrune = false;
 int res = computeRelaxations( f_RelaxationSolvers , rootNode , minimizing ,
                               bestBound , bestSol , incumbentCell ,
                               deadline , toPrune , branchSolver , relAcc ,
                               absAcc , 1 , &globalMutex );
 if( toPrune || ( res != ThinComputeInterface::kOK ) )
  return( res );
 res = computeHeuristic( f_HeuristicSolvers , rootNode , minimizing ,
                         bestBound , bestSol , incumbentCell , deadline ,
                         toPrune , 1 , &globalMutex );
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
 * concurrently (serial-tree explorations only); \p incumbentMutex !=
 * nullptr serializes the incumbent updates (parallel-tree exploration). Sets
 * \p toPrune / \p branchSolver, and marks the node infeasible [see
 * Node::set_infeasible()] when its relaxation has no solution.
 *  @return the sol_type [see Solver.h] of the evaluation */

static int evaluateNode(
                  Node * node , std::list< ChangeSolver * > & solvers ,
                  std::vector< std::pair< RelaxationSolver * ,
                               Solver * > > * relaxation ,
                  std::vector< std::pair< ChangeSolver * ,
                               Solver * > > * heuristic ,
                  bool minimizing , double & bestBound , Solution * & bestSol ,
                  std::atomic< double > * incumbentCell ,
                  const Deadline deadline ,
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

 // computeRelaxations() / computeHeuristic() dispatch internally on
 // nThreads: a worker of the parallel tree exploration (incumbentMutex set)
 // evaluates its per-node solvers serially, otherwise nThreads > 1 runs them
 // in parallel
 int res = computeRelaxations( relaxation , node , minimizing , bestBound ,
                               bestSol , incumbentCell , deadline , toPrune ,
                               branchSolver , relAcc ,
                               absAcc , nThreads , incumbentMutex );
 if( ( res != Solver::kOK ) || toPrune )
  return( res );

 return( computeHeuristic( heuristic , node , minimizing , bestBound ,
                           bestSol , incumbentCell , deadline , toPrune ,
                           nThreads , incumbentMutex ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OPEN-LIST DISCIPLINES -------------------------*/
/*--------------------------------------------------------------------------*/
// the three open-list disciplines, as thin adapters that give the standard
// container adapters a common interface [see OpenList]: a LIFO stack (depth-
// first), a FIFO queue (breadth-first), a relaxation-bound priority queue
// (best- first). The storage is entirely std::stack / std::queue /
// std::priority_- queue, besides what each keeps to give the best bound of
// its nodes [see OpenList::best_bound()]

namespace {

/// the order of a best-first open list: the worse relaxation bound first,
/// or the worse estimate first in a best-estimate search
/** A plain function object rather than a std::function, so that the heap
 * calls it directly; worse() compares two values, and also serves the
 * lists that keep the best bound of their nodes. */

struct BoundCmp {
 bool minimizing;
 bool estimate = false;  ///< by Node::estimate_key() rather than the bound
 [[nodiscard]] bool worse( double ka , double kb ) const {
  return( minimizing ? ka > kb : ka < kb );
  }
 bool operator()( ExploringNode * a , ExploringNode * b ) const {
  const double ka = estimate ? a->estimate_key( minimizing ) : a->get_bound();
  const double kb = estimate ? b->estimate_key( minimizing ) : b->get_bound();
  return( worse( ka , kb ) );
  }
 };

/*--------------------------------------------------------------------------*/
/// LIFO open list: a std::stack, giving depth-first exploration
/** Each entry carries the best bound of the nodes from the bottom of the
 * stack up to it, so that the best bound is that of the top. */

class StackOpenList final : public OpenList {

 public:

 explicit StackOpenList( BoundCmp cmp ) : f_cmp( cmp ) {}

 [[nodiscard]] bool empty( void ) const override { return( c.empty() ); }

 void push( ExploringNode * node ) override {
  const double b = node->get_bound();
  c.emplace( node , ( c.empty() || f_cmp.worse( c.top().second , b ) )
                    ? b : c.top().second );
  }

 ExploringNode * pop( void ) override {
  auto node = c.top().first;
  c.pop();
  return( node );
  }

 [[nodiscard]] bool isLIFO( void ) const override { return( true ); }

 [[nodiscard]] bool best_bound( double & bound ) const override {
  if( c.empty() )
   return( false );
  bound = c.top().second;
  return( true );
  }

 private:

 BoundCmp f_cmp;  ///< the order of the bounds

 std::stack< std::pair< ExploringNode * , double > > c;

 };  // end( class( StackOpenList ) )

/*--------------------------------------------------------------------------*/
/// FIFO open list: a std::queue, giving breadth-first exploration
/** Besides the queue of the nodes, the deque of the bounds that may still
 * become the best one, i.e., of those no worse than any pushed after them:
 * the best bound is the front, pushing a node drops from the back the
 * bounds worse than its own, and popping a node drops the front if it is
 * its bound, which takes constant amortized time. */

class QueueOpenList final : public OpenList {

 public:

 explicit QueueOpenList( BoundCmp cmp ) : f_cmp( cmp ) {}

 [[nodiscard]] bool empty( void ) const override { return( c.empty() ); }

 void push( ExploringNode * node ) override {
  c.push( node );
  const double b = node->get_bound();
  while( ( ! best.empty() ) && f_cmp.worse( best.back() , b ) )
   best.pop_back();
  best.push_back( b );
  }

 ExploringNode * pop( void ) override {
  auto node = c.front();
  c.pop();
  if( best.front() == node->get_bound() )
   best.pop_front();
  return( node );
  }

 [[nodiscard]] bool best_bound( double & bound ) const override {
  if( best.empty() )
   return( false );
  bound = best.front();
  return( true );
  }

 private:

 BoundCmp f_cmp;  ///< the order of the bounds

 std::queue< ExploringNode * > c;

 std::deque< double > best;  ///< the bounds that may become the best

 };  // end( class( QueueOpenList ) )

/*--------------------------------------------------------------------------*/
/// relaxation-bound priority open list: a std::priority_queue, best-first
/** Ordered by the estimates, the list keeps also the multiset of the bounds
 * of its nodes, whose best is the best bound. */

class PriorityOpenList final : public OpenList {

 using Cmp = BoundCmp;

 public:

 explicit PriorityOpenList( Cmp cmp )
  : by_bound( ! cmp.estimate ) , minimizing( cmp.minimizing ) ,
    c( std::move( cmp ) ) {}

 [[nodiscard]] bool empty( void ) const override { return( c.empty() ); }

 void push( ExploringNode * node ) override {
  c.push( node );
  if( ! by_bound )
   bounds.insert( node->get_bound() );
  }

 ExploringNode * pop( void ) override {
  auto node = c.top();
  c.pop();
  if( ! by_bound )
   if( auto it = bounds.find( node->get_bound() ) ; it != bounds.end() )
    bounds.erase( it );
  return( node );
  }

 [[nodiscard]] bool next_bound( double & bound ) const override {
  if( ( ! by_bound ) || c.empty() )
   return( false );
  bound = c.top()->get_bound();
  return( true );
  }

 [[nodiscard]] bool best_bound( double & bound ) const override {
  if( by_bound )
   return( next_bound( bound ) );
  if( bounds.empty() )
   return( false );
  bound = minimizing ? *bounds.begin() : *bounds.rbegin();
  return( true );
  }

 private:

 bool by_bound;    ///< ordered by the bounds, not by the estimates
 bool minimizing;  ///< the sense of the bounds

 std::multiset< double > bounds;  ///< of the nodes, if not by_bound

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
 * taken. With a finite quotient [see BranchAndXSolver::dblPlungeQuot] the
 * dive also ends when the bound of the child it would follow is farther from
 * the best bound than that share of the gap to the incumbent. */

class DiveOpenList final : public OpenList {

 using Cmp = BoundCmp;

 public:

 explicit DiveOpenList( Cmp cmp , double quot = Inf< double >() ,
                        const std::atomic< double > * inc = nullptr )
  : f_cmp( cmp ) , f_quot( quot ) , f_inc( inc ) , best( std::move( cmp ) ) {}

 [[nodiscard]] bool empty( void ) const override {
  return( best.empty() && children.empty() );
  }

 void push( ExploringNode * node ) override { children.push_back( node ); }

 ExploringNode * pop( void ) override {
  if( ( ! children.empty() ) && plunge( children.back()->get_bound() ) ) {
   // the loop is isLIFO(), so it pushed the children in reverse branching
   // order: back() is the first (most promising) branch, the dive follows it
   auto node = children.back();
   children.pop_back();
   for( auto sibling : children )
    best.push( sibling );
   children.clear();
   return( node );
   }
  for( auto child : children )   // the dive ends here
   best.push( child );
  children.clear();
  auto node = best.top();
  best.pop();
  return( node );
  }

 [[nodiscard]] bool isLIFO( void ) const override { return( true ); }

 /// the best of the bound of the top of the queue, ordered by the bounds,
 /// and of those of the children pushed since the last pop()
 [[nodiscard]] bool best_bound( double & bound ) const override {
  if( empty() )
   return( false );
  bool first = best.empty();
  if( ! first )
   bound = best.top()->get_bound();
  for( auto child : children )
   if( first || f_cmp.worse( bound , child->get_bound() ) ) {
    bound = child->get_bound();
    first = false;
    }
  return( true );
  }

 private:

 /// whether the dive goes on to a child of the given bound: always with an
 /// infinite quotient or no incumbent, and otherwise if the bound is within
 /// the quotient of the gap from the best bound to the incumbent
 [[nodiscard]] bool plunge( double bound ) const {
  if( ( f_quot == Inf< double >() ) || ( ! f_inc ) )
   return( true );
  const double z = f_inc->load();
  double b0;
  if( ( ! std::isfinite( z ) ) || ( ! best_bound( b0 ) ) ||
      ( ! std::isfinite( b0 ) ) )
   return( true );
  return( f_cmp.minimizing ? bound <= b0 + f_quot * ( z - b0 )
                           : bound >= b0 - f_quot * ( b0 - z ) );
  }

 Cmp f_cmp;  ///< the order of the bounds

 double f_quot;  ///< the quotient of the dives [see dblPlungeQuot]

 const std::atomic< double > * f_inc;  ///< the incumbent, if any

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
 /// sense of the Block, NaN if none, and the solution
 double found = std::numeric_limits< double >::quiet_NaN();
 std::unique_ptr< Solution > solution;          ///< that solution
 };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// a node of the tree as PBeFS handles it: a BXLink whose copy is a node of
/// its own, with the Change, the parent and the solution of the original,
/// so that evaluating and branching a copy (as the dive of the SelfSplit
/// does) leaves the original as it was; moves keep the BXLink

struct BXNode {

 std::shared_ptr< BXLink > p;                   ///< the node

 BXNode( void ) = default;

 BXNode( std::shared_ptr< BXLink > l ) noexcept : p( std::move( l ) ) {}

 BXNode( const BXNode & n ) : p( copy( n.p ) ) {}

 BXNode( BXNode && ) noexcept = default;

 BXNode & operator=( const BXNode & n ) {
  if( this != & n )
   p = copy( n.p );
  return( *this );
  }

 BXNode & operator=( BXNode && ) noexcept = default;

 BXLink * operator->( void ) const noexcept { return( p.get() ); }

 private:

 static std::shared_ptr< BXLink > copy(
                                     const std::shared_ptr< BXLink > & l ) {
  if( ! l )
   return( nullptr );
  auto c = std::make_shared< BXLink >();
  c->change = l->change;
  c->parent = l->parent;
  c->depth = l->depth;
  c->bound = l->bound;
  c->estimate = l->estimate;
  c->evaluated = l->evaluated;
  c->pruned = l->pruned;
  c->branches = l->branches;
  c->estimates = l->estimates;
  c->found = l->found;
  if( l->solution )
   c->solution.reset( l->solution->clone() );
  return( c );
  }

 };  // end( struct( BXNode ) )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// what the workers share: the data of the BranchAndXSolver

struct BXContext {

 using Relaxations = std::vector< std::pair< RelaxationSolver * ,
                                             Solver * > >;
 using Heuristics = std::vector< std::pair< ChangeSolver * , Solver * > >;

 BXContext( double rt , double at , bool min , bool e , bool l ,
            std::mutex & gm , Deadline d , double f )
  : relTol( rt ) , absTol( at ) , minimizing( min ) , eager( e ) ,
    lifo( l ) , globalMutex( gm ) , deadline( d ) , rootBound( f ) {}

 double relTol;                          ///< the relative tolerance
 double absTol;                          ///< the absolute tolerance
 bool minimizing;                        ///< the sense of the Block
 bool eager;                             ///< the eager bounding protocol
 bool lifo;                              ///< the first branch last
 std::mutex & globalMutex;               ///< of the BranchAndXSolver
 Deadline deadline;                      ///< the end of the time budget
 double rootBound;                       ///< the bound the root starts from
 /// called with the bound of each node set aside unresolved, under
 /// unresolvedMutex [see BranchAndXSolver::set_aside()]
 std::function< void( double ) > set_aside;
 std::mutex unresolvedMutex;             ///< of set_aside
 std::vector< Relaxations * > relaxation;  ///< of each worker
 std::vector< Heuristics * > heuristic;    ///< of each worker
 std::atomic< int > error{ Solver::kOK };  ///< the first error, if any
 /// called with the bound of the root once it is evaluated, which is the
 /// bound of the search until the search gives a better one [see
 /// BranchAndXSolver::publish_bound()]
 std::function< void( double ) > root_evaluated;

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

 /// the incumbent z of the search, in the sense of PBeFS, which the
 /// worker evaluates its next node against
 void cutoff( double z ) {
  f_best = f_c.sign() * z;
  f_cell.store( f_best );
  }

 /// the evaluation of the node n, done now unless already done, whose
 /// feasible value is the one found by the node
 PB::Evaluation evaluate( BXNode & n ) {
  if( ! n->evaluated )
   evaluate_here( n.p );
  PB::Evaluation e;
  e.bound = f_c.sign() * n->bound;
  e.infeasible = n->pruned && std::isnan( n->found );
  if( ! std::isnan( n->found ) )
   e.feasible = f_c.sign() * n->found;
  return( e );
  }

 /// the children of the evaluated node n, evaluated now with the eager
 /// protocol (and then given only if they survive), the first branch
 /// given last when the open set is a stack
 template< class F >
 void branch( const BXNode & n , F && f ) {
  std::vector< std::shared_ptr< BXLink > > kids;
  kids.reserve( n->branches.size() );
  for( std::size_t i = 0 ; i < n->branches.size() ; ++i ) {
   auto k = std::make_shared< BXLink >();
   k->change = n->branches[ i ];
   k->parent = n.p;
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
   f( BXNode( std::move( k ) ) , b );
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
  // against the incumbent of the worker [see cutoff()]
  const int res = evaluateNode( & tmp , f_solvers , f_c.relaxation[ f_i ] ,
                                f_c.heuristic[ f_i ] , f_c.minimizing ,
                                f_best , f_sol , & f_cell , f_c.deadline ,
                                toPrune , bs , 1 , f_c.relTol , f_c.absTol ,
                                f_c.globalMutex , nullptr , false );
  tmp.set_toFather( nullptr );
  n->evaluated = true;
  n->bound = tmp.get_bound();
  if( ( ! n->change ) && ( res == Solver::kOK ) && f_c.root_evaluated )
   f_c.root_evaluated( n->bound );
  if( f_sol ) {                         // found by this node
   n->solution.reset( f_sol );
   n->found = f_best;
   f_sol = nullptr;
   }
  if( inconclusive( res ) ) {           // set aside, as if fenced
   {
    std::lock_guard< std::mutex > guard( f_c.unresolvedMutex );
    f_c.set_aside( n->bound );
    }
   n->pruned = true;
   return;
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
              cannot_improve( n->bound , f_cell.load() ,
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
 // the first Solver giving the undo Change [see moveBetweenNodes()]; the
 // common ancestor is found climbing from n only until a node sits in the
 // path at the position of its depth, the path being a chain from the root,
 // so that a move costs its own steps rather than the depth of n
 void move_to( const std::shared_ptr< BXLink > & n ) {
  v_down.clear();
  const std::shared_ptr< BXLink > * x = & n;
  for( ; * x && ( * x )->change ; x = & ( * x )->parent ) {
   const auto d = std::size_t( ( * x )->depth );
   if( ( d >= 1 ) && ( d <= v_path.size() ) &&
       ( v_path[ d - 1 ].first.get() == x->get() ) )
    break;
   v_down.push_back( x );
   }
  const std::size_t k = ( * x && ( * x )->change )
                        ? std::size_t( ( * x )->depth ) : 0;
  while( v_path.size() > k ) {
   Change * undo = v_path.back().second;
   for( auto s : f_solvers )
    s->apply( undo , false );
   delete undo;
   v_path.pop_back();
   }
  for( auto it = v_down.rbegin() ; it != v_down.rend() ; ++it ) {
   const std::shared_ptr< BXLink > & y = ** it;
   Change * undo = nullptr;
   bool first = true;
   for( auto s : f_solvers ) {
    if( first )
     undo = s->apply( y->change.get() , true );
    else
     s->apply( y->change.get() , false );
    first = false;
    }
   v_path.emplace_back( y , undo );
   }
  }

 BXContext & f_c;                               ///< the shared data
 unsigned f_i;                                  ///< the index of the worker
 std::list< ChangeSolver * > f_solvers;         ///< its Solver
 /// the nodes from the root to that of the Solver, with the undo Changes
 std::vector< std::pair< std::shared_ptr< BXLink > , Change * > > v_path;
 /// the nodes that a move descends, deepest first, as the pointers that
 /// hold them (n, or the parent of a node of the move) [see move_to()]
 std::vector< const std::shared_ptr< BXLink > * > v_down;
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

 using node_type = BXNode;
 static constexpr PB::Sense sense = PB::Sense::minimize;

 BXContext & c;                                 ///< the shared data

 [[nodiscard]] node_type root( void ) const {
  auto r = std::make_shared< BXLink >();
  r->bound = c.rootBound;
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

double BranchAndXSolver::start_incumbent( void ) const
{
 return( f_Block->get_objective_sense() == Objective::eMax ? lwCutOff
                                                           : upCutOff );
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
 f_evaluated_nodes.store( 0 , std::memory_order_relaxed );
 f_bounded_children = 0;
 // no bound until an exploration says otherwise [see explore()]
 f_open_bound = ( f_Block->get_objective_sense() == Objective::eMax )
                ? Inf< double >() : - Inf< double >();
 v_unresolved.clear();
 std::mutex globalMutex;

 // per-solve budgets and tolerances from the inherited standard parameters
 nodeBudget = get_int_par( intMaxNodes );
 timeBudget = get_dbl_par( dblMaxTime );
 // a time budget beyond what the clock can count (say, an infinite one)
 // must not be duration_cast, as it would overflow into a deadline in the
 // past: it simply means no deadline at all
 {
  const auto now = std::chrono::steady_clock::now();
  const double most = std::chrono::duration< double >(
   std::chrono::steady_clock::time_point::max() - now ).count();
  f_deadline = ( timeBudget < most / 2 ) ?
   now + std::chrono::duration_cast< std::chrono::steady_clock::duration >(
                             std::chrono::duration< double >( timeBudget ) ) :
   std::chrono::steady_clock::time_point::max();
  }
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
                       ( ! astar ) && ( ! probing ) ) );

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
      parallelDFS || astar || probing )
   discardRetainedTree();

  bestBound = no_incumbent();
  // no incumbent yet [see no_incumbent() and GlobalInformation::str_Incumbent]
  // but the cutoff, if any, which the nodes are pruned against as if it were
  // one [see start_incumbent()]
  f_incumbentCell->store( start_incumbent() );
  // nor any bound of the search, nor nodes set aside [see publish_bound()]
  f_bestBoundCell->store( - no_incumbent() );
  f_unresolvedBound.store( no_incumbent() , std::memory_order_relaxed );
  // no bound known but those of the nodes [see f_probeBound]
  f_probeBound = ( f_Block->get_objective_sense() == Objective::eMax )
                 ? Inf< double >() : - Inf< double >();
  // one exploration [see runExploration()], or the sequence of those of
  // the cutoff probing [see probeCutOff()]; an exception (say, a Block that
  // the parallel exploration cannot copy) leaves the Solver unlocked, and
  // in error
  try {
   f_state = probing ? probeCutOff( globalMutex )
                     : runExploration( globalMutex );
   }
  catch( ... ) {
   f_state = kError;
   unlock();
   throw;
   }

  // an exploration that has been completed without ever finding a feasible
  // solution proves that there is none, or, with a cutoff, none better than
  // it, which the bound on the side of the nodes then says [see get_lb()]
  if( ( f_state == kOK ) && ( bestBound == no_incumbent() ) ) {
   if( start_incumbent() == no_incumbent() )
    f_state = kInfeasible;
   else {
    f_state = kCutOff;
    f_open_bound = start_incumbent();
    }
   }

  // the bound of the whole solve [see get_lb()], unless it has failed
  if( f_state != kError ) {
   const bool minimizing =
    ( f_Block->get_objective_sense() != Objective::eMax );
   publish_bound( minimizing ? get_lb() : get_ub() , minimizing );
   }
  }

 // The summary also covers an eager root closed without entering explore(),
 // and the parallel search, whose worker nodes do not use the serial CSV log.
 if( f_log ) {
  ( *f_log ) << "[BranchAndXSolver] nodes_evaluated="
             << get_evaluated_nodes() << std::endl;
  if( f_bounded_children )
   ( *f_log ) << "[BranchAndXSolver] children_pruned_by_bound="
              << f_bounded_children << std::endl;
  }

 changes = 0;
 unlock();
 return( f_state );

 }  // end( BranchAndXSolver::compute )


/*--------------------------------------------------------------------------*/

int BranchAndXSolver::runExploration( std::mutex & globalMutex )
{
 // one exploration, parameterized by the open set that the strategy asks
 // for [see treeSolve()]; the only separate ones are the parallel
 // depth-first, which runs its own workers [see ParallelDFSSolve()], and
 // the branch-and-bound of PBeFS [see PBeFSSolve()]
 if( treeExplorer != OwnExplorer )
  return( PBeFSSolve( globalMutex , get_int_par( intMaxThread ) ) );
 if( ( solveType == DFS ) && ( get_int_par( intMaxThread ) > 1 ) )
  return( ParallelDFSSolve( globalMutex , get_int_par( intMaxThread ) ) );
 return( treeSolve( globalMutex ) );
 }

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::probeCutOff( std::mutex & globalMutex )
{
 // the values are those of the Block times sign, so that they are minimized
 const bool minimizing =
  ( f_Block->get_objective_sense() != Objective::eMax );
 const double sign = minimizing ? 1 : -1;

 // the bound D of the root, out of its relaxations [see computeRelaxations()]
 DFSNode root( nullptr );
 root.initializeBound( minimizing );
 root.set_evaluated( true );
 bool toPrune = false;
 RelaxationSolver * branchSolver = nullptr;
 f_evaluated_nodes.fetch_add( 1 , std::memory_order_relaxed );
 const int res = computeRelaxations( &f_RelaxationSolvers , &root ,
                                     minimizing , bestBound , bestSolution ,
                                     f_incumbentCell , f_deadline , toPrune ,
                                     branchSolver , relTol , absTol , 1 ,
                                     &globalMutex );
 const double D = sign * root.get_bound();

 // nothing to probe if the root is infeasible, fenced or inconclusive, or
 // has no finite bound: the exploration says which is the case
 if( ( res != kOK ) || toPrune || ( ! std::isfinite( D ) ) )
  return( runExploration( globalMutex ) );

 // the tolerance by which a node is pruned against the value v [see
 // cannot_improve()]
 const double absEps = ( absTol == Inf< double >() ) ? 0 : absTol;
 const auto eps = [ & ]( double v ) {
  return( std::max( absEps , relTol * std::max( std::abs( v ) , 1.0 ) ) );
  };

 double delta = probingGap;
 if( probing == RelativeGap )
  delta *= std::max( std::abs( D ) , 1.0 );
 double target = D + delta;
 double proven = D;

 for( int refuted = 0 ; ; ++refuted ) {
  // the incumbent: the best of the cutoff and of the solution found, if any
  const double zbar = sign * f_incumbentCell->load();
  const bool last = ( target >= zbar ) ||
                    ( ( probingRuns > 0 ) && ( refuted >= probingRuns ) );
  if( ! last )
   f_incumbentCell->store( sign * target );
  f_probeBound = sign * proven;
  const double found = bestBound;
  const int z = runExploration( globalMutex );
  if( last )
   return( z );
  if( bestBound != found )  // a solution better than the target: optimal
   return( z );
  // the target is no incumbent: the cell is given back the true one
  f_incumbentCell->store( sign * zbar );
  if( z != kOK ) {
   // stopped by a budget: the nodes fenced against the target have bounds
   // no better than it, the open ones say the rest
   f_open_bound = sign * std::min( sign * f_open_bound , target );
   return( z );
   }
  // no solution better than the target, which is now a valid bound but
  // for the tolerance the nodes were fenced with
  proven = target - eps( target );
  delta *= probingFactor;
  target += delta;
  }
 }  // end( BranchAndXSolver::probeCutOff )

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
 if( ! bounds.empty() )        // so do the nodes set aside unresolved
  add_unresolved_bounds( minimizing );
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
  F->set_evaluated( true );
  bool toPrune = false;
  RelaxationSolver * nodeBranchSolver = nullptr;
  f_evaluated_nodes.fetch_add( 1 , std::memory_order_relaxed );
  res = computeRelaxations( &f_RelaxationSolvers , F , minimizing , bestBound ,
                            bestSolution , f_incumbentCell , f_deadline ,
                            toPrune , nodeBranchSolver ,
                            relTol , absTol , maxThreadForSolvers );
  if( ( res == Solver::kOK ) && ( ! toPrune ) )
   res = computeHeuristic( &f_HeuristicSolvers , F , minimizing , bestBound ,
                           bestSolution , f_incumbentCell , f_deadline ,
                           toPrune , maxThreadForSolvers );
  if( inconclusive( res ) ) {    // set aside, still in the frontier
   set_aside( F->get_bound() , minimizing );
   subOptimalNodes.push_back( F );
   res = Solver::kOK;
   continue;
   }
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
 std::vector< ExploringNode * > kept;  // surviving children of a node
 // whether the loop stopped at a node whose bound cannot improve the
 // incumbent, which the open set gives best bound first: none of the open
 // nodes can, and they are all fenced as they would be at their pop
 bool fenced = false;

 while( ( ! open.empty() ) && ( nodeBudget > 1 ) &&
        ( timeBudget > std::chrono::duration< double >(
           std::chrono::high_resolution_clock::now() - start ).count() ) ) {
  // the bound of the search: that of the open nodes, the next one included,
  // which the node being evaluated has until the next extraction
  if( double best ; open.best_bound( best ) )
   publish_bound( best , minimizing );
  if( double next ; open.next_bound( next ) &&
      cannot_improve( next , f_incumbentCell->load() , minimizing , relTol ,
                      absTol ) ) {
   fenced = true;
   break;
   }
  nodeBudget--;
  ++iterations;
  oldNode = currentNode;
  currentNode = open.pop();

  // a node that can no longer improve the incumbent is pruned where it is,
  // the solvers staying where they are, unless that would delete the node
  // they are on; a retained tree keeps it in the fenced frontier instead,
  // which needs the solvers on it
  if( ( ! retain ) &&
      cannot_improve( currentNode->get_bound() , f_incumbentCell->load() ,
                      minimizing , relTol , absTol ) &&
      ExploringNode::can_prune_away( currentNode , oldNode ) ) {
   if( f_log )
    logNode( f_log , iterations , currentNode , bestBound ,
	     std::chrono::duration< double >(
	      std::chrono::high_resolution_clock::now() - start ).count() ,
	     0 , true );
   ExploringNode::prune_away( currentNode );
   currentNode = oldNode;
   continue;
   }

  if( currentNode->get_parent() )   // the root: the solvers are already there
   ExploringNode::moveBetweenNodes( oldNode , currentNode , solvers );

  // with the lazy protocol the node is evaluated now, upon extraction,
  // rather than when it was created [see intBoundingProtocol]; a node is
  // evaluated at most once, whichever the protocol, and not at all if the
  // bound it carries (that of its parent) can no longer improve the
  // incumbent
  bool prunedHere = false;
  if( ( ! currentNode->is_evaluated() ) &&
      cannot_improve( currentNode->get_bound() , f_incumbentCell->load() ,
                      minimizing , relTol , absTol ) ) {
   prunedHere = true;
   if( f_log )
    logNode( f_log , iterations , currentNode , bestBound ,
	    std::chrono::duration< double >(
	     std::chrono::high_resolution_clock::now() - start ).count() ,
	    0 , true );
   }
  else if( ! currentNode->is_evaluated() ) {
   decltype( start ) evalStart;  // the clock is read for the log only
   if( f_log )
    evalStart = std::chrono::high_resolution_clock::now();
   f_evaluated_nodes.fetch_add( 1 , std::memory_order_relaxed );
   res = evaluateNode( currentNode , *solvers , &f_RelaxationSolvers ,
                       &f_HeuristicSolvers , minimizing , bestBound ,
                       bestSolution , f_incumbentCell , f_deadline ,
                       prunedHere ,
                       branchSolver ,
                       maxThreadForSolvers , relTol , absTol , globalMutex ,
                       nullptr , false );
   if( inconclusive( res ) ) {     // set aside, as if fenced
    set_aside( currentNode->get_bound() , minimizing );
    prunedHere = true;
    res = Solver::kOK;
    }
   if( res != Solver::kOK ) {       // an infeasible root ends up here
    // leave the solvers at the root, the Changes of the node undone
    ExploringNode::moveBetweenNodes( currentNode , rootNode , solvers );
    discardTree( rootNode , open , solvers , currentNode );
    return( res );
    }
   if( ! prunedHere )
    currentNode->obtainBranchList( branchSolver );
   if( f_log )
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
   kept.clear();         // survivors, pushed after the loop
   // with the bounds of the children given by the RelaxationSolver that
   // branched [see BranchBounder], the children are taken best bound first
   // and, unless the tree is retained, those known to be empty or that
   // cannot improve the incumbent are pruned without being created
   const bool bounded = currentNode->has_child_bounds();
   std::vector< std::size_t > order( branches.size() );
   std::iota( order.begin() , order.end() , std::size_t( 0 ) );
   if( bounded )
    std::stable_sort( order.begin() , order.end() ,
                      [ & ]( std::size_t a , std::size_t b ) {
     const double ba = currentNode->child_bound( a , minimizing );
     const double bb = currentNode->child_bound( b , minimizing );
     return( minimizing ? ba < bb : ba > bb );
     } );
   for( const std::size_t i : order ) {
    const double childBound = currentNode->child_bound( i , minimizing );
    if( bounded && ( ! retain ) &&
        ( ( childBound == ( minimizing ? Inf< double >()
                                       : - Inf< double >() ) ) ||
          cannot_improve( childBound , f_incumbentCell->load() , minimizing ,
                          relTol , absTol ) ) ) {
     ++f_bounded_children;
     delete branches[ i ];
     branches[ i ] = nullptr;
     continue;
     }
    ExploringNode * new_node = new ExploringNode( branches[ i ] ,
                                                  currentNode ,
                                                  currentNode->get_level()
                                                  + 1 , ++counter );
    // the relaxation bound of the parent is valid for the child, whose region
    // is a part of that of the parent: the child starts from it, or from the
    // better one given by the RelaxationSolver that branched, if any [see
    // BranchBounder], and its own evaluation can only improve it; the
    // estimate is the one given by the RelaxationSolver that branched, if any
    // [see BranchEstimator]
    new_node->set_bound( childBound );
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
    decltype( start ) evalStart;  // the clock is read for the log only
    if( f_log )
     evalStart = std::chrono::high_resolution_clock::now();
    f_evaluated_nodes.fetch_add( 1 , std::memory_order_relaxed );
    res = evaluateNode( new_node , *solvers , &f_RelaxationSolvers ,
                        &f_HeuristicSolvers , minimizing , bestBound ,
                        bestSolution , f_incumbentCell , f_deadline , toPrune ,
                        branchSolver ,
                        maxThreadForSolvers , relTol , absTol , globalMutex );
    if( inconclusive( res ) ) {    // set aside, as if fenced
     set_aside( new_node->get_bound() , minimizing );
     toPrune = true;
     res = Solver::kOK;
     }
    if( res != Solver::kOK ) {
     // leave the solvers at the root, the Changes of the child undone
     ExploringNode::moveBetweenNodes( new_node , rootNode , solvers );
     discardTree( rootNode , open , solvers , currentNode );
     return( res );
     }
    if( f_log )
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
 // the incumbent, whatever the budgets that remain; those left when the
 // loop stopped at a fenced node are all fenced
 bool complete = true;
 std::vector< ExploringNode * > leftOpen;
 // the global relaxation bound: the best among those of the open nodes
 f_open_bound = minimizing ? Inf< double >() : - Inf< double >();
 while( ! open.empty() ) {
  auto node = open.pop();
  leftOpen.push_back( node );
  if( fenced )             // fenced, as at their pop: no open bound
   continue;
  if( ! cannot_improve( node->get_bound() , f_incumbentCell->load() ,
                        minimizing , relTol , absTol ) )
   complete = false;
  f_open_bound = minimizing ? std::min( f_open_bound ,
                                        node->get_bound() )
                            : std::max( f_open_bound ,
                                        node->get_bound() );
  }

 // move the :ChangeSolver back to the root
 while( currentNode->get_toFather() ) {
  for( const auto s : *solvers )
   s->apply( currentNode->get_toFather() , false );
  currentNode = currentNode->get_parent();
  }

 if( retain ) {
  // retain the tree for future reoptimizations: the nodes left in the open
  // set (early stops) are open work, and those fenced are of the fenced
  // frontier, hence both part of the frontier to re-seed
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
 if( ! complete ) {
  add_unresolved_bounds( minimizing );
  return( nodeBudget <= 1 ? Solver::kStopIter : Solver::kStopTime );
  }
 // some nodes set aside may have been stopped by the time budget
 if( add_unresolved_bounds( minimizing ) )
  return( out_of_time() ? Solver::kStopTime : Solver::kLowPrecision );
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
 const BoundCmp cmp{ minimizing };
 std::unique_ptr< OpenList > openPtr;
 switch( solveType ) {
  case( DFS ):
   openPtr.reset( new StackOpenList( cmp ) );
   break;
  case( BFS ):
   openPtr.reset( new QueueOpenList( cmp ) );
   break;
  case( BestFS ):
   openPtr.reset( new PriorityOpenList( cmp ) );
   break;
  case( BestFSDive ):
   openPtr.reset( new DiveOpenList( cmp , plungeQuot , f_incumbentCell ) );
   break;
  case( BestEstimate ):
   openPtr.reset( new PriorityOpenList( BoundCmp{ minimizing , true } ) );
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
 const bool retain = ( reoptimize > 0 ) && ( ! probing );
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
   rootNode->set_bound( f_probeBound );
   res = Solver::kOK;
   }
  else {
   f_evaluated_nodes.fetch_add( 1 , std::memory_order_relaxed );
   res = initializeRoot( &f_RelaxationSolvers , &f_HeuristicSolvers ,
                         rootNode , minimizing , bestBound , bestSolution ,
                         f_incumbentCell , f_deadline , branchSolver ,
                         globalMutex , relTol , absTol , f_probeBound );
   }
  if( res == Solver::kOK )
   open.push( rootNode );
  currentNode = rootNode;
  }

 if( res != Solver::kOK ) {      // the bound of an inconclusive root stands
  discardTree( rootNode , open , solvers ,
               inconclusive( res ) ? rootNode : nullptr );
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
 f_evaluated_nodes.fetch_add( 1 , std::memory_order_relaxed );
 int res = evaluateNode( currentNode , solvers , &relaxation , &heuristic ,
                         minimizing , bestBound , bestSolution ,
                         f_incumbentCell , f_deadline , toPrune ,
                         branchSolver , 1 ,
                         relTol , absTol , incumbentMutex , &incumbentMutex );
 if( inconclusive( res ) ) {       // set aside, as if fenced
  {
   std::lock_guard< std::mutex > guard( incumbentMutex );
   set_aside( currentNode->get_bound() , minimizing );
   }
  moveSolverToFather( currentNode , &solvers );
  return( Solver::kOK );
  }
 // whatever the outcome, the Solver go back to the father, so that a stopped
 // worker leaves them at the root of the subtree it has claimed
 if( res != ThinComputeInterface::kOK ) {
  moveSolverToFather( currentNode , &solvers );
  return( res );
  }

 if( ( ! toPrune ) &&
     ( ! cannot_improve( currentNode->get_bound() , f_incumbentCell->load() ,
                         minimizing , relTol , absTol ) ) ) {
  auto branches = branchSolver->branch();
  for( std::size_t i = 0 ; i < branches.size() ; ++i ) {
   // the node lives as long as its exploration, i.e., in this frame
   DFSNode new_node( branches[ i ] , ++nameCounter );
   new_node.set_bound( currentNode->get_bound() );  // inherited
   int RV = workerDFS( & new_node , solvers , relaxation , heuristic ,
                       minimizing , incumbentMutex , nodeBdg , deadline ,
                       nameCounter );
   if( RV != Solver::kOK ) {      // the branches not explored are deleted
    while( ++i < branches.size() )
     delete branches[ i ];
    moveSolverToFather( currentNode , &solvers );
    return( RV );
    }
   }
  }

 moveSolverToFather( currentNode , &solvers );
 return( Solver::kOK );

 }  // end( BranchAndXSolver::workerDFS )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::ParallelDFSSolve( std::mutex & globalMutex , int K )
{
 // the copies of the Block and their Solver last as long as this solve
 struct DropWorkers {
  BranchAndXSolver * bxs;
  ~DropWorkers() { bxs->dropWorkerSolvers(); }
  } dropWorkers{ this };
 createWorkerSolvers( K );

 // the serial Solver set drives the ramp-up
 std::list< ChangeSolver * > solvers;
 bool minimizing;
 initializeVariables( &solvers , &f_RelaxationSolvers , &f_HeuristicSolvers ,
                      minimizing );

 const auto deadline = f_deadline;
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
 f_evaluated_nodes.fetch_add( 1 , std::memory_order_relaxed );
 int res = initializeRoot( &f_RelaxationSolvers , &f_HeuristicSolvers ,
                           rootNode , minimizing , bestBound , bestSolution ,
                           f_incumbentCell , f_deadline , branchSolver ,
                           globalMutex , relTol , absTol , f_probeBound );
 if( res != Solver::kOK ) {
  if( inconclusive( res ) )        // the bound of the root, if any
   f_open_bound = rootNode->get_bound();
  delete rootNode;
  return( res );
  }

 std::deque< ExploringNode * > pool;
 pool.push_back( rootNode );
 ExploringNode * currentNode = rootNode;
 const std::size_t target = std::size_t( 4 * K );

 // the bound of the search [see publish_bound()]: the best among those of
 // the nodes of the pool and of the subtrees that the workers have claimed,
 // each of which is no worse than that of the node at its root; it is
 // published by the ramp-up before each extraction, and by the workers
 // when they claim a subtree and when they end its exploration, under the
 // lock of the pool
 const double noBound = minimizing ? Inf< double >() : - Inf< double >();
 std::vector< double > claimed( std::size_t( K ) , noBound );
 const auto publishPool = [ & ]( void ) {
  double best = noBound;
  for( auto n : pool )
   best = minimizing ? std::min( best , n->get_bound() )
                     : std::max( best , n->get_bound() );
  for( auto b : claimed )
   best = minimizing ? std::min( best , b ) : std::max( best , b );
  if( best != noBound )
   publish_bound( best , minimizing );
  };

 while( ( ! pool.empty() ) && ( pool.size() < target ) &&
        ( --nodeBdg > 0 ) &&
        ( std::chrono::steady_clock::now() < deadline ) ) {
  publishPool();
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
    f_evaluated_nodes.fetch_add( 1 , std::memory_order_relaxed );
    res = evaluateNode( new_node , solvers , &f_RelaxationSolvers ,
                        &f_HeuristicSolvers , minimizing , bestBound ,
                        bestSolution , f_incumbentCell , f_deadline ,
                        toPrune , nodeBranchSolver ,
                        1 , relTol , absTol , globalMutex );
    if( inconclusive( res ) ) {    // set aside, as if fenced
     set_aside( new_node->get_bound() , minimizing );
     toPrune = true;
     res = Solver::kOK;
     }
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
  // an idle worker waits on this for the pool to be fed or for the last
  // busy worker to finish, both of which happen under poolMutex
  std::condition_variable poolChanged;

  auto workerLoop = [ & ]( int w ) {
   auto & ws = v_workerSolvers[ w ];
   std::list< ChangeSolver * > wsolvers;
   for( auto & pr : ws.relaxation )
    wsolvers.push_back( pr.first );
   for( auto & pr : ws.heuristic )
    wsolvers.push_back( pr.first );
   int wNameCounter = ( w + 1 ) * 10000000;  // disjoint per-worker names
   std::vector< Change * > path;  // the Changes from a node to the root

   while( result.load() == Solver::kOK ) {
    ExploringNode * node;
    {
     std::unique_lock< std::mutex > guard( poolMutex );
     // the pool may be re-fed by the work-sharing of a busy worker: only an
     // empty pool with NO busy worker means there is no work left
     poolChanged.wait( guard , [ & ] {
      return( ( ! pool.empty() ) || ( busy.load() == 0 ) ||
              ( result.load() != Solver::kOK ) );
      } );
     if( pool.empty() || ( result.load() != Solver::kOK ) )
      break;
     node = pool.front();
     pool.pop_front();
     --poolSize;
     ++busy;
     claimed[ std::size_t( w ) ] = node->get_bound();
     publishPool();
     }

    // claim the subtree: re-apply the path of Changes from the root
    path.clear();
    for( auto n = node ; n ; n = n->get_parent() )
     if( n->get_f_change() )
      path.push_back( n->get_f_change() );
    for( auto it = path.rbegin() ; it != path.rend() ; ++it )
     for( auto s : wsolvers )
      s->apply( *it , false );

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
       f_evaluated_nodes.fetch_add( 1 , std::memory_order_relaxed );
       int RV = computeRelaxations( &ws.relaxation , new_node , minimizing ,
                                    bestBound , bestSolution ,
                                    f_incumbentCell , f_deadline , toPrune ,
                                    nodeBranchSolver , relTol , absTol ,
                                    1 , &incumbentMutex );
       if( RV == Solver::kOK && ! toPrune )
        RV = computeHeuristic( &ws.heuristic , new_node , minimizing ,
                               bestBound , bestSolution ,
                               f_incumbentCell , f_deadline , toPrune ,
                               1 , &incumbentMutex );
       if( inconclusive( RV ) ) {  // set aside, as if fenced
        {
         std::lock_guard< std::mutex > guard( incumbentMutex );
         set_aside( new_node->get_bound() , minimizing );
         }
        toPrune = true;
        RV = Solver::kOK;
        }
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
        {
         std::lock_guard< std::mutex > guard( poolMutex );
         pool.push_back( new_node );
         ++poolSize;
         }
        poolChanged.notify_one();
        }
       else                      // pruned or fenced: nothing to donate
        delete new_node;
       continue;
       }

      // the DFSNode takes the ownership of the branching Change: null the
      // entry so that ~Node does not double-delete it on the final cleanup
      DFSNode new_node( br , ++wNameCounter );
      new_node.set_bound( node->get_bound() );  // inherited
      br = nullptr;
      first = false;
      int RV = workerDFS( & new_node , wsolvers , ws.relaxation ,
                          ws.heuristic , minimizing , incumbentMutex ,
                          nodeBdg , deadline , wNameCounter );
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

    {
     std::lock_guard< std::mutex > guard( poolMutex );
     --busy;
     claimed[ std::size_t( w ) ] = noBound;
     if( result.load() == Solver::kOK )
      publishPool();
     }
    poolChanged.notify_all();
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
 // all explored but the nodes set aside unresolved, which alone bound what
 // is left [see add_unresolved_bounds()]
 if( result.load() == Solver::kOK ) {
  f_open_bound = minimizing ? Inf< double >() : - Inf< double >();
  if( add_unresolved_bounds( minimizing ) )
   return( out_of_time() ? Solver::kStopTime : Solver::kLowPrecision );
  }
 return( result.load() );

 }  // end( BranchAndXSolver::ParallelDFSSolve )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::PBeFSSolve( std::mutex & globalMutex , int K )
{
 K = std::max( 1 , K );
 // the copies of the Block and their Solver last as long as this solve
 struct DropWorkers {
  BranchAndXSolver * bxs;
  ~DropWorkers() { bxs->dropWorkerSolvers(); }
  } dropWorkers{ this };
 createWorkerSolvers( Index( K ) );
 const bool minimizing =
  ( f_Block->get_objective_sense() != Objective::eMax );
 BXContext c( relTol , absTol , minimizing , boundingProtocol == Eager ,
              ( solveType == DFS ) || ( solveType == BestFSDive ) ,
              globalMutex , f_deadline , f_probeBound );
 for( auto & ws : v_workerSolvers ) {
  c.relaxation.push_back( & ws.relaxation );
  c.heuristic.push_back( & ws.heuristic );
  }
 c.set_aside = [ this , minimizing ]( double b ) {
  set_aside( b , minimizing );
  };
 // the frontier is within PBeFS, which gives its bound (in the sense of
 // the minimization) by on_bound when it improves and at the end, and may
 // give none before the end: until it does, the bound is that of the root;
 // both are given by the workers, which publish_bound() allows
 c.root_evaluated = [ this , minimizing ]( double b ) {
  publish_bound( b , minimizing );
  };

 // the strategy as a selection of PBeFS: the dive of BestFSDive goes as
 // far as dblPlungeQuot says, by default down to a leaf whatever the bound
 // of the child
 PB::ParallelBBOptions o;
 const auto select = [ this ]( auto & b , int m ) {
  switch( m ) {
   case( DFS ):    b.selection = PB::Selection::depth_first; break;
   case( BFS ):    b.selection = PB::Selection::breadth_first; break;
   case( BestFS ): b.selection = PB::Selection::best_first; break;
   case( BestEstimate ): b.selection = PB::Selection::best_estimate; break;
   default:
    b.selection = PB::Selection::plunging;
    b.plunge_gap = plungeQuot;
   }
  };
 select( o , solveType );
 // the racers, at most one per worker [see vintPBeFSRace]
 for( std::size_t k = 0 ; ( k < pbefsRace.size() ) && ( k < Index( K ) ) ;
      ++k ) {
  o.race.emplace_back();
  select( o.race.back() , pbefsRace[ k ] );
  }
 if( pbefsRaceNodes > 0 )
  o.race_nodes = std::size_t( pbefsRaceNodes );
 o.rel_gap = relTol;
 o.abs_gap = absTol == Inf< double >() ? 0 : absTol;
 if( nodeBudget < INT_MAX )
  o.node_limit = std::size_t( std::max( 1 , nodeBudget ) );
 if( timeBudget < Inf< double >() )
  o.time_limit = timeBudget;
 o.threads = unsigned( K );
 if( pbefsSpawnDepth > 0 )
  o.spawn_depth = unsigned( pbefsSpawnDepth );
 if( pbefsBudget > 0 )
  o.budget = std::uint64_t( pbefsBudget );
 o.split_shared = pbefsSplitShared;
 if( f_incumbentCell->load() != no_incumbent() )
  o.incumbent = c.sign() * f_incumbentCell->load();
 o.on_bound = [ this , & c , minimizing ]( double b ) {
  publish_bound( c.sign() * b , minimizing );
  };

 const BXTree t{ c };
 const auto s = PB::parallel_branch_and_bound( t ,
                                               PB::bb_scheme( pbefsScheme ) ,
                                               o );

 if( c.error.load() != kOK )
  return( c.error.load() );
 // the solution is that of the incumbent node
 if( s.node && ( *s.node )->solution ) {
  bestBound = ( *s.node )->found;
  f_incumbentCell->store( bestBound );
  delete bestSolution;
  bestSolution = ( *s.node )->solution.release();
  }
 // the nodes set aside unresolved [see set_aside()] bound what is left,
 // as the open ones do [see add_unresolved_bounds()]
 switch( s.status ) {
  case( PB::BBStatus::node_limit ):
   f_open_bound = c.sign() * s.bound;
   add_unresolved_bounds( minimizing );
   return( kStopIter );
  case( PB::BBStatus::time_limit ):
   f_open_bound = c.sign() * s.bound;
   add_unresolved_bounds( minimizing );
   return( kStopTime );
  default:
   f_open_bound = minimizing ? Inf< double >() : - Inf< double >();
   if( add_unresolved_bounds( minimizing ) )
    return( out_of_time() ? kStopTime : kLowPrecision );
   return( kOK );
  }
 }  // end( BranchAndXSolver::PBeFSSolve )

/*--------------------------------------------------------------------------*/

void BranchAndXSolver::checkPBeFSScheme( const std::string & name )
{
 try {
  ( void ) PB::bb_scheme( name );
  }
 catch( const std::invalid_argument & ) {
  throw( std::invalid_argument( "BranchAndXSolver::set_par: PBeFS has no "
                                "scheme " + name ) );
  }
 }

/*--------------------------------------------------------------------------*/

void BranchAndXSolver::set_aside( double bound , bool minimizing )
{
 v_unresolved.push_back( bound );
 // the best of the bounds set aside is lowered lock-free, the caller's lock
 // being on v_unresolved only
 double cur = f_unresolvedBound.load( std::memory_order_relaxed );
 while( ( minimizing ? bound < cur : bound > cur ) &&
        ( ! f_unresolvedBound.compare_exchange_weak( cur , bound ) ) ) {}
 }

/*--------------------------------------------------------------------------*/

void BranchAndXSolver::publish_bound( double bound , bool minimizing )
{
 const double u = f_unresolvedBound.load( std::memory_order_relaxed );
 const double z = f_incumbentCell->load( std::memory_order_relaxed );
 bound = minimizing ? std::min( { bound , u , z } )
                    : std::max( { bound , u , z } );
 // only improved, by whichever thread finds a better one
 double cur = f_bestBoundCell->load( std::memory_order_relaxed );
 while( ( minimizing ? bound > cur : bound < cur ) &&
        ( ! f_bestBoundCell->compare_exchange_weak( cur , bound ) ) ) {}
 }

/*--------------------------------------------------------------------------*/

bool BranchAndXSolver::add_unresolved_bounds( bool minimizing )
{
 bool open = false;
 for( const auto b : v_unresolved ) {
  f_open_bound = minimizing ? std::min( f_open_bound , b )
                            : std::max( f_open_bound , b );
  open |= ! cannot_improve( b , f_incumbentCell->load() , minimizing ,
                            relTol , absTol );
  }
 return( open );
 }

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

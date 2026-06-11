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
 * Solver. WARNING: WORK IN PROGRESS - see the file-level comment of
 * BranchAndXSolver.h.
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

// register GroupChange to the Change factory

SMSpp_insert_in_factory_cpp_0( GroupChange );

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
 * dual bound, the branching solver, and - when a true solution improves it
 * - the incumbent bestBound / bestSol; sets \p toPrune when the node can
 * be discarded.
 *  @return the sol_type [see Solver.h] of the computation */

static int computeRelaxations(
                  std::vector< RelaxationSolver * > * f_RelaxationSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol ,
                  bool & toPrune , RelaxationSolver * & branchSolver ,
                  double relAcc = 0 , double absAcc = 0 )
{
 for( auto s : *f_RelaxationSolvers ) {
  auto z = s->compute();
  if( z == Solver::kInfeasible ) {
   if( ! currentNode->get_toFather() )    // the root is infeasible
    return( Solver::kInfeasible );
   toPrune = true;
   break;
   }
  if( z != ThinComputeInterface::kOK )
   return( z );

  // see if the primal bound improves thanks to a true solution
  if( s->has_true_var_solution() ) {
   double primal_bound = minimizing ? s->get_true_ub() : s->get_true_lb();
   // TODO: evaluate primal_bound - 1e-6 to ward off numerical issues
   if( minimizing ? primal_bound < bestBound : primal_bound > bestBound ) {
    bestBound = primal_bound;
    delete bestSol;
    bestSol = s->get_true_solution();
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
                  double & bestBound , Solution * & bestSol , bool & toPrune )
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
   bestBound = primal_bound;
   delete bestSol;
   bestSol = s->get_Solution();
   }
  }
 return( Solver::kOK );
 }

/*--------------------------------------------------------------------------*/
/// parallel version of computeRelaxations(): NOT implemented yet

static int computeRelaxationsParallel(
                  std::vector< RelaxationSolver * > * f_RelaxationSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol , bool & toPrune ,
                  RelaxationSolver * & branchSolver ,
                  std::mutex & globalMutex , int maxThreads )
{
 // TODO: implement the parallel computation of the relaxations
 return( computeRelaxations( f_RelaxationSolvers , currentNode , minimizing ,
                             bestBound , bestSol , toPrune , branchSolver )
         );
 }

/*--------------------------------------------------------------------------*/
/// parallel version of computeHeuristic(): NOT implemented yet

static int computeHeuristicParallel(
                  std::vector< ChangeSolver * > * f_HeuristicSolvers ,
                  Node * currentNode , const bool minimizing ,
                  double & bestBound , Solution * & bestSol , bool & toPrune ,
                  std::mutex & globalMutex , int maxThreads )
{
 // TODO: implement the parallel computation of the heuristics
 return( computeHeuristic( f_HeuristicSolvers , currentNode , minimizing ,
                           bestBound , bestSol , toPrune ) );
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
/*----------------- METHODS OF BranchAndXSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int BranchAndXSolver::compute( bool changedvars )
{
 // TODO: understand whether this has to be done in parallel, too
 lock();
 process_outstanding_Modification();
 int old_state = f_state;
 if( f_state == kStillRunning )
  return( kError );
 f_state = kStillRunning;
 std::mutex globalMutex;

 // per-solve budgets and tolerances from the inherited standard parameters
 nodeBudget = get_int_par( intMaxIter );
 timeBudget = get_dbl_par( dblMaxTime );
 relTol = get_dbl_par( dblRelAcc );
 absTol = get_dbl_par( dblAbsAcc );

 if( changes == 0 )         // nothing changed since the last solve
  f_state = old_state;
 else {
  // TODO: cases changes == 1 / 2 / 3 (objective and/or r.h.s. changes only)
  // should reoptimize the tree out of the subOptimalNodes / infeasibleNodes
  // / integerNodes lists rather than solving from scratch; until that is
  // implemented every change triggers a full solve
  bestBound = ( f_Block->get_objective_sense() == Objective::eMax )
              ? - Inf< double >() : Inf< double >();
  switch( solveType ) {
   case( DFS ): {
    int counter = 0;
    DFSNode * root_node = new DFSNode( nullptr );
    f_state = DFSSolve( globalMutex , root_node , counter );
    delete root_node;
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

int BranchAndXSolver::DFSSolve( std::mutex & globalMutex ,
                                    DFSNode * currentNode ,
                                    int & nameCounter )
{
 nodeBudget--;
 std::list< ChangeSolver * > solvers;
 bool minimizing;
 initializeVariables( &solvers , &f_RelaxationSolvers , &f_HeuristicSolvers ,
                      minimizing );
 auto start = std::chrono::high_resolution_clock::now();
 if( nodeBudget <= 0 )
  return( Solver::kStopIter );
 if( timeBudget <= 0 )
  return( Solver::kStopTime );

 // move all the :ChangeSolver to the current node
 // TODO: do it in parallel when maxThreadForSolvers > 1
 moveSolverToSon( currentNode , &solvers );

 currentNode->initializeBound( minimizing );
 RelaxationSolver * branchSolver = nullptr;
 // TODO: evaluate computing relaxations and heuristics in parallel

 bool toPrune = false;
 int res;
 if( maxThreadForSolvers == 1 )
  res = computeRelaxations( &f_RelaxationSolvers , currentNode , minimizing ,
                            bestBound , bestSolution , toPrune ,
                            branchSolver , relTol , absTol );
 else
  res = computeRelaxationsParallel( &f_RelaxationSolvers , currentNode ,
                                    minimizing , bestBound , bestSolution ,
                                    toPrune , branchSolver , globalMutex ,
                                    maxThreadForSolvers );
 if( res != ThinComputeInterface::kOK )
  return( res );
 if( toPrune ) {
  for( const auto s : solvers )
   s->apply( currentNode->get_toFather() , false );
  timeBudget -= std::chrono::duration< double >(
              std::chrono::high_resolution_clock::now() - start ).count();
  if( f_log )
   *f_log << "Node " << currentNode->get_name()
          << " pruned by relaxation." << std::endl;
  return( ThinComputeInterface::kOK );
  }

 // heuristic solvers
 if( maxThreadForSolvers == 1 )
  res = computeHeuristic( &f_HeuristicSolvers , currentNode , minimizing ,
                          bestBound , bestSolution , toPrune );
 else
  res = computeHeuristicParallel( &f_HeuristicSolvers , currentNode ,
                                  minimizing , bestBound , bestSolution ,
                                  toPrune , globalMutex ,
                                  maxThreadForSolvers );
 if( res != ThinComputeInterface::kOK )
  return( res );
 if( toPrune ) {
  // TODO: do it in parallel when maxThreadForSolvers > 1
  for( const auto s : solvers )
   s->apply( currentNode->get_toFather() , false );
  timeBudget -= std::chrono::duration< double >(
              std::chrono::high_resolution_clock::now() - start ).count();
  if( f_log )
   *f_log << "Node " << currentNode->get_name()
          << " pruned by heuristic." << std::endl;
  return( ThinComputeInterface::kOK );
  }

 std::string log = "Exploring node " +
  std::to_string( currentNode->get_name() ) + " with dual bound " +
  std::to_string( currentNode->get_dual_bound() ) + ", global upper bound " +
  std::to_string( bestBound ) + "\n";
 log += "Compute time " + std::to_string( std::chrono::duration< double >(
         std::chrono::high_resolution_clock::now() - start ).count() ) +
        " seconds.\n";

 if( ! cannot_improve( currentNode->get_dual_bound() , bestBound ,
                       minimizing , relTol , absTol ) ) {
  auto branches = branchSolver->branch();
  for( auto br : branches ) {
   DFSNode * new_node = new DFSNode( br , ++nameCounter );
   int RV = DFSSolve( globalMutex , new_node , nameCounter );
   delete new_node;
   if( RV != Solver::kOK )
    return( RV );
   }
  }

 if( currentNode->get_toFather() ) {
  for( const auto s : solvers )
   s->apply( currentNode->get_toFather() , false );
  timeBudget -= std::chrono::duration< double >(
              std::chrono::high_resolution_clock::now() - start ).count();
  }

 log += "Time to evaluate his subtree: " +
  std::to_string( std::chrono::duration< double >(
   std::chrono::high_resolution_clock::now() - start ).count() ) +
  " seconds.\n";
 if( f_log )
  *f_log << log;
 return( Solver::kOK );

 }  // end( BranchAndXSolver::DFSSolve )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::BestFirstSolve( std::mutex & globalMutex )
{
 auto * solvers = new std::list< ChangeSolver * >();
 bool minimizing;
 initializeVariables( solvers , &f_RelaxationSolvers , &f_HeuristicSolvers ,
                      minimizing );

 // initialize the priority queue (by best dual bound first)
 auto start = std::chrono::high_resolution_clock::now();
 int counter = 0;
 auto cmp = [ minimizing ]( ExploringNode * a , ExploringNode * b ) {
  return( minimizing ? a->get_dual_bound() > b->get_dual_bound()
                     : a->get_dual_bound() < b->get_dual_bound() );
  };
 auto pq = std::priority_queue< ExploringNode * ,
                                std::vector< ExploringNode * > ,
                                decltype( cmp ) >( cmp );
 ExploringNode * rootNode = new ExploringNode( nullptr , nullptr , 0 );
 RelaxationSolver * branchSolver = nullptr;

 // cleanup of the whole tree, for the error exits
 auto cleanupAll = [ & ]( ExploringNode * root ,
                          std::priority_queue< ExploringNode * ,
                           std::vector< ExploringNode * > ,
                           decltype( cmp ) > & queue ,
                          std::list< ChangeSolver * > * solversPtr ) {
  while( ! queue.empty() ) {
   ExploringNode * node = queue.top();
   queue.pop();
   delete node;
   }
  std::function< void( ExploringNode * ) > deleteTree =
   [ & ]( ExploringNode * node ) {
    if( ! node )
     return;
    for( auto child : node->get_children() )
     deleteTree( child );
    if( node != root )
     delete node;
    };
  deleteTree( root );
  delete root;
  delete solversPtr;
  };

 int res = initializeRoot( &f_RelaxationSolvers , &f_HeuristicSolvers ,
                           rootNode , minimizing , bestBound , bestSolution ,
                           branchSolver , globalMutex );
 if( res != Solver::kOK ) {
  cleanupAll( rootNode , pq , solvers );
  return( res );
  }
 pq.push( rootNode );
 ExploringNode * currentNode = rootNode;
 ExploringNode * oldNode = nullptr;

 // TODO: understand how to make this parallel, given that the tree is
 // continuously modified by the exploration
 while( ( ! pq.empty() ) && ( nodeBudget > 1 ) &&
        ( timeBudget > std::chrono::duration< double >(
           std::chrono::high_resolution_clock::now() - start ).count() ) ) {
  nodeBudget--;
  oldNode = currentNode;
  currentNode = pq.top();
  pq.pop();
  // TODO: evaluate whether this check is needed
  if( currentNode->get_parent() )
   ExploringNode::moveBetweenNodes( oldNode , currentNode , solvers );

  // branching and evaluation of the new children
  if( ! cannot_improve( currentNode->get_dual_bound() , bestBound ,
                        minimizing , relTol , absTol ) ) {
   auto & branches = currentNode->getBranches();
   // TODO: make this a parallel loop
   for( auto br : branches ) {
    ExploringNode * new_node = new ExploringNode( br , currentNode ,
                                                  currentNode->get_level()
                                                  + 1 , ++counter );
    moveSolverToSon( new_node , solvers );
    new_node->initializeBound( minimizing );
    bool toPrune = false;
    res = computeRelaxations( &f_RelaxationSolvers , new_node , minimizing ,
                              bestBound , bestSolution , toPrune ,
                              branchSolver , relTol , absTol );
    if( res != Solver::kOK ) {
     cleanupAll( rootNode , pq , solvers );
     return( res );
     }
    if( ! toPrune ) {
     res = computeHeuristic( &f_HeuristicSolvers , new_node , minimizing ,
                             bestBound , bestSolution , toPrune );
     if( res != Solver::kOK ) {
      cleanupAll( rootNode , pq , solvers );
      return( res );
      }
     }
    // keep the node only if its dual bound can improve the incumbent
    if( ( ! toPrune ) &&
        ( ! cannot_improve( new_node->get_dual_bound() , bestBound ,
                            minimizing , relTol , absTol ) ) ) {
     new_node->obtainBranchList( branchSolver );
     pq.push( new_node );
     currentNode->get_children().push_back( new_node );
     moveSolverToFather( new_node , solvers );
     }
    else {
     moveSolverToFather( new_node , solvers );
     *( std::find( branches.begin() , branches.end() ,
                   new_node->get_f_change() ) ) = nullptr;
     delete new_node;
     }
    }
   branches.erase( std::remove( branches.begin() , branches.end() ,
                                nullptr ) , branches.end() );
   if( branches.empty() && currentNode->get_toFather() )
    currentNode = ExploringNode::prune( currentNode , solvers );
   }
  else if( currentNode->get_toFather() )
   currentNode = ExploringNode::prune( currentNode , solvers );
  }

 // move the :ChangeSolver back to the root
 while( currentNode->get_toFather() ) {
  for( const auto s : *solvers )
   s->apply( currentNode->get_toFather() , false );
  currentNode = currentNode->get_parent();
  }

 // cleanup: delete all the nodes left in the priority queue ...
 std::vector< ExploringNode * > nodesToDelete;
 while( ! pq.empty() ) {
  nodesToDelete.push_back( pq.top() );
  pq.pop();
  }
 for( auto node : nodesToDelete )
  delete node;

 // ... then, recursively, all the remaining children
 std::function< void( ExploringNode * ) > deleteTree =
  [ & ]( ExploringNode * node ) {
   if( ! node )
    return;
   for( auto child : node->get_children() )
    deleteTree( child );
   if( node != rootNode )
    delete node;
   };
 deleteTree( rootNode );
 delete rootNode;
 delete solvers;

 if( nodeBudget <= 0 )
  return( Solver::kStopIter );
 if( timeBudget <= std::chrono::duration< double >(
                 std::chrono::high_resolution_clock::now() - start ).count()
     )
  return( Solver::kStopTime );
 return( kOK );

 }  // end( BranchAndXSolver::BestFirstSolve )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::BFSSolve( std::mutex & globalMutex )
{
 auto * solvers = new std::list< ChangeSolver * >();
 bool minimizing;
 initializeVariables( solvers , &f_RelaxationSolvers , &f_HeuristicSolvers ,
                      minimizing );
 int counter = 0;
 auto start = std::chrono::high_resolution_clock::now();
 std::queue< ExploringNode * > nodesQueue;
 ExploringNode * rootNode = new ExploringNode( nullptr , nullptr , 0 );
 RelaxationSolver * branchSolver = nullptr;
 int res = initializeRoot( &f_RelaxationSolvers , &f_HeuristicSolvers ,
                           rootNode , minimizing , bestBound , bestSolution ,
                           branchSolver , globalMutex );
 if( res != Solver::kOK ) {
  delete rootNode;
  delete solvers;
  return( res );
  }
 nodesQueue.push( rootNode );
 ExploringNode * currentNode = rootNode;
 ExploringNode * oldNode = nullptr;

 while( ( ! nodesQueue.empty() ) && ( nodeBudget > 1 ) &&
        ( timeBudget > std::chrono::duration< double >(
           std::chrono::high_resolution_clock::now() - start ).count() ) ) {
  nodeBudget--;
  oldNode = currentNode;
  currentNode = nodesQueue.front();
  nodesQueue.pop();
  if( currentNode->get_parent() )
   ExploringNode::moveBetweenNodes( oldNode , currentNode , solvers );

  // branching: only the best node can branch
  if( ! cannot_improve( currentNode->get_dual_bound() , bestBound ,
                        minimizing , relTol , absTol ) ) {
   auto & branches = currentNode->getBranches();
   for( Change * br : branches ) {
    ExploringNode * new_node = new ExploringNode( br , currentNode ,
                                                  currentNode->get_level()
                                                  + 1 , ++counter );
    moveSolverToSon( new_node , solvers );
    new_node->initializeBound( minimizing );
    bool toPrune = false;
    RelaxationSolver * nodeBranchSolver = nullptr;

    int res = computeRelaxations( &f_RelaxationSolvers , new_node ,
                                  minimizing , bestBound , bestSolution ,
                                  toPrune , nodeBranchSolver , relTol ,
                                  absTol );
    if( res != Solver::kOK )
     return( res );
    if( ! toPrune ) {
     res = computeHeuristic( &f_HeuristicSolvers , new_node , minimizing ,
                             bestBound , bestSolution , toPrune );
     if( res != Solver::kOK )
      return( res );
     }
    if( ( ! toPrune ) &&
        ( ! cannot_improve( new_node->get_dual_bound() , bestBound ,
                            minimizing , relTol , absTol ) ) ) {
     new_node->obtainBranchList( nodeBranchSolver );
     nodesQueue.push( new_node );
     currentNode->get_children().push_back( new_node );
     moveSolverToFather( new_node , solvers );
     }
    else {
     moveSolverToFather( new_node , solvers );
     *( std::find( branches.begin() , branches.end() ,
                   new_node->get_f_change() ) ) = nullptr;
     delete new_node;
     }
    }
   branches.erase( std::remove( branches.begin() , branches.end() ,
                                nullptr ) , branches.end() );
   if( branches.empty() && currentNode->get_toFather() )
    currentNode = ExploringNode::prune( currentNode , solvers );
   }
  else if( currentNode->get_toFather() )
   currentNode = ExploringNode::prune( currentNode , solvers );
  }

 // move the :ChangeSolver back to the root
 while( currentNode->get_toFather() ) {
  for( const auto s : *solvers )
   s->apply( currentNode->get_toFather() , false );
  currentNode = currentNode->get_parent();
  }

 if( nodeBudget <= 0 )
  return( Solver::kStopIter );
 if( timeBudget <= std::chrono::duration< double >(
                 std::chrono::high_resolution_clock::now() - start ).count()
     )
  return( Solver::kStopTime );
 return( kOK );

 }  // end( BranchAndXSolver::BFSSolve )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

// TODO: classify the Modification (objective only / r.h.s. only / both /
// general) so that compute() can reoptimize rather than solve from scratch

void BranchAndXSolver::process_outstanding_Modification( void )
{
 Lst_sp_Mod v_mod_tmp;              // temporary list of modifications

 // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) );

 for( auto mod : v_mod )
  v_mod_tmp.push_back( mod );       // copy v_mod in v_mod_tmp

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );  // release lock

 if( ! v_mod_tmp.empty() )
  changes = 4;                      // for now: any change = solve again

 v_mod_tmp.clear();

 }  // end( BranchAndXSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/
/*------------------ End File BranchAndXSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/

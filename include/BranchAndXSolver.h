/*--------------------------------------------------------------------------*/
/*---------------------- File BranchAndXSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class BranchAndXSolver, which implements the
 * Solver concept [see Solver.h] for a Relaxation-Agnostic Branch-and-X (RABaX)
 * Solver: the problem encoded by the Block is solved with an enumerative
 * algorithm built on top of the ChangeSolver / RelaxationSolver concepts [see
 * ChangeSolver.h]. The nodes of the enumeration tree are identified by the
 * Change that generates them from their parent, the attached
 * RelaxationSolver(s) provide relaxation bounds, true (primal) solutions and
 * the branching Changes, and the attached heuristic ChangeSolver(s), if any,
 * provide further primal bounds. The tree can be explored depth-first,
 * breadth-first, best-first or by the estimates that the RelaxationSolver
 * gives of its children (see intSolveMethod); nodes are pruned by the
 * standard inherited tolerances (dblRelAcc / dblAbsAcc), and intMaxNodes / the
 * inherited dblMaxTime bound the number of explored nodes and the total time.
 *
 * The solver is generic ("X" is for Bound / Cut / Price): nothing in here
 * depends on the specific :Block being solved, all the problem-specific
 * knowledge lives in the attached :ChangeSolver, which provides the relaxation
 * bounds and the branching Changes (any tightening it does --- preprocessing,
 * reduced-cost fixing, cuts --- it does on its own terms, reading the search-
 * global data through a GlobalInformation [see ChangeSolver.h]). The tree is
 * explored serially or in parallel, both at the tree level (see
 * ParallelDFSSolve()) and across the solvers of a single node (see
 * intThreadForDifferentSolvers), and can be retained across re-solves for
 * reoptimization (see intReoptimize).
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
 * \copyright &copy; by Antonio Frangioni, Filippo Magi, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BranchAndXSolver
 #define __BranchAndXSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <limits>
#include <list>
#include <mutex>
#include <queue>
#include <vector>

#include "BlockSolverConfig.h"

#include "BranchingRule.h"

#include "ChangeSolver.h"

#include "Solution.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

 class Node;            // forward declaration of Node
 class DFSNode;         // forward declaration of DFSNode
 class ExploringNode;   // forward declaration of ExploringNode
 class OpenList;        // forward declaration of OpenList

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS BranchAndXSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a generic Branch-and-Bound Solver driven by ChangeSolver
/** The BranchAndXSolver solves the problem
 *  \f[
 *    z^* = \max \{ f(x) : x \in X \}
 *  \f]
 * encoded by its Block (minimization being the mirror image, handled by the
 * sense of the Objective) by recursively partitioning \f$ X \f$ and bounding
 * each part through a relaxation \f$ \hat{X} \supseteq X \f$ over which
 * maximizing \f$ f \f$ is "easy", the relaxation being solved by the
 * attached :RelaxationSolver [see ChangeSolver.h]. In the following
 * \f$ \bar{z} \f$ is the value of the best feasible solution found so far
 * (the incumbent), and, for a node \f$ \nu \f$ of the enumeration tree,
 * \f$ X_\nu \subseteq X \f$ is the subproblem it represents and
 * \f$ D( \nu ) \geq \max \{ f(x) : x \in X_\nu \} \f$ the (dual) bound
 * given by its relaxation.
 *
 * \par The enumeration tree.
 * A node \f$ \nu \f$ is identified by the Change \f$ c_\nu \f$ that
 * generates it from its father \f$ p( \nu ) \f$, and it owns \f$ c_\nu \f$,
 * the undo \f$ c_\nu^{-1} \f$ that goes back to the father, and its bound;
 * hence \f$ X_\nu \f$ is the result of applying to \f$ X \f$ the Changes
 * along the path from the root, i.e.,
 *  \f[
 *    X_\nu = ( c_\nu \circ c_{p( \nu )} \circ \dots \circ c_{r'} ) X
 *  \f]
 * where \f$ r' \f$ is the child of the root on that path. A GroupChange
 * \f$ g = \langle c_1 , \dots , c_k \rangle \f$ is applied as the single
 * Change \f$ c_k \circ \dots \circ c_1 \f$, and its undo is the composition
 * of the undos in reverse order,
 *  \f[
 *    g^{-1} = c_1^{-1} \circ \dots \circ c_k^{-1} \;,
 *  \f]
 * which is what makes a branching decision touching several things at once
 * a Change like any other, and the round trip down and up a subtree exact.
 * Since the Changes of the children of a node partition its subproblem
 * (the RelaxationSolver must ensure that \f$ X_\nu = \bigcup_i
 * X_{\nu c_i} \f$ for those returned by branch()), the leaves of any cut
 * through the tree partition \f$ X \f$, and no feasible solution is lost.
 *
 * \par Moving along the tree.
 * One :ChangeSolver visits many nodes without being rebuilt: it is moved
 * from a node \f$ s \f$ to a node \f$ t \f$ along the unique path through
 * their lowest common ancestor \f$ a \f$, applying the undos from
 * \f$ s \f$ up to \f$ a \f$ and then the Changes from \f$ a \f$ down to
 * \f$ t \f$. A :ChangeSolver may intercept the Changes it understands and
 * apply them to its own internal state rather than to the Block, which is
 * where the efficiency of the enumeration (possibly) comes from, since
 * moving along the tree then never locks nor reads the Block. The
 * depth-first exploration needs no explicit links (the recursion is the
 * path to the root), while the others keep the father of each node, so as
 * to jump between any two open nodes.
 *
 * \par Bounding and pruning.
 * Since \f$ X_\nu \subseteq \hat{X}_\nu \f$, no solution in \f$ X_\nu \f$
 * is better than \f$ D( \nu ) \f$; hence, \f$ \nu \f$ is pruned when its
 * bound cannot improve the incumbent within the inherited tolerances
 * \f$ \varepsilon_a \f$ = dblAbsAcc and \f$ \varepsilon_r \f$ = dblRelAcc,
 * i.e., when
 *  \f[
 *    D( \nu ) \leq \bar{z} + \max \{ \varepsilon_a \,,\, \varepsilon_r
 *                                    \max \{ | \bar{z} | , 1 \} \}
 *  \f]
 * (the inequality reversed for minimization), as well as when its
 * relaxation is empty. Note that the default dblAbsAcc = +INF is taken as
 * "not active", i.e., as \f$ \varepsilon_a = 0 \f$: plugged in literally it
 * would prune every node as soon as there is an incumbent.
 *
 * \par The Solver of a node.
 * A node is evaluated by a set of inner Solver, all described by the
 * BlockSolverConfig [see strNameOfBlockSolverConfigurationFile]: one or more
 * :RelaxationSolver \f$ \mathcal{R} \f$, each giving a bound
 * \f$ D_R( \nu ) \f$ and possibly a solution feasible for the true problem
 * of value \f$ \underline{z}_R \f$, and zero or more plain :ChangeSolver
 * \f$ \mathcal{H} \f$, the primal heuristics, each giving only a feasible
 * solution \f$ x_H \f$ and never a bound. The bound of the node is the
 * tightest of those of the relaxations,
 *  \f[
 *    D( \nu ) = \min \{ D_R( \nu ) : R \in \mathcal{R} \}
 *  \f]
 * (the max for minimization), and the relaxation attaining it is the one
 * that branches the node; the incumbent is updated with the best solution
 * found by any Solver of the node,
 *  \f[
 *    \bar{z} \leftarrow \max \{ \bar{z} \,,\,
 *                     \max \{ \underline{z}_R : R \in \mathcal{R} \} \,,\,
 *                     \max \{ f( x_H ) : H \in \mathcal{H} \} \} \;.
 *  \f]
 * A heuristic therefore changes neither the partition nor the bounds, only
 * \f$ \bar{z} \f$; yet, a good incumbent found early is what makes the
 * pruning test effective, and it is often decisive in practice.
 *
 * \par Exploration strategies.
 * If \f$ \mathcal{O} \f$ is the set of the open nodes (generated but not
 * yet branched), a strategy [see intSolveMethod] is only the rule choosing
 * the next node in \f$ \mathcal{O} \f$: the deepest one for depth-first, the
 * shallowest one for breadth-first, the one with the best \f$ D( \nu ) \f$
 * for best-first (diving to a leaf below it in BestFSDive), and the one with
 * the best estimate of the value under it for best-estimate. Best-first
 * explores the fewest nodes for a given relaxation, since no node with
 * \f$ D( \nu ) \leq z^* \f$ is ever branched, but it stores the whole of
 * \f$ \mathcal{O} \f$; depth-first only stores the current path, hence its
 * memory is linear in the depth, which makes it the safer choice on deep
 * trees and the one run by the workers of the parallel exploration.
 *
 * \par Parallel exploration.
 * Two orthogonal kinds of parallelism are available. Across the tree [see
 * intMaxThread and ParallelDFSSolve()], after an ordered ramp-up \f$ K \f$
 * workers, each driving its own copy of the inner Solver (so that applying
 * the Changes needs no synchronization), repeatedly claim the oldest open
 * subtree and explore it depth-first, sharing only the incumbent and the
 * pool of open subtrees, and a worker donates a branch of its subtree when
 * the pool runs low; claiming in the order of discovery is what keeps the
 * runs replicable, since a random distribution of the work is known to make
 * them erratic on knapsack-like problems. Across the Solver of a node [see
 * intThreadForDifferentSolvers], the compute() of the several relaxations
 * and heuristics are run concurrently while the reduction of their results
 * (incumbent, tightest bound, pruning) stays serial, and therefore the
 * outcome is the same as that of the serial code. Of course, the first kind
 * is subject to the usual anomalies of parallel Branch-and-Bound: a worker
 * finding a good incumbent early may let the others prune subtrees that the
 * sequential search would enter (and vice versa), hence speedups may be
 * superlinear as well as below one.
 *
 * \par Tightening.
 * Any tightening of the relaxation (preprocessing, reduced-cost fixing,
 * cuts) is done by the RelaxationSolver on its own terms, which reads what
 * is global to the search (the incumbent, the pool of the globally valid
 * cuts) in a GlobalInformation [see str_GlobalCuts]; a tightening is valid
 * if it does not remove any solution of \f$ X_\nu \f$ better than the
 * incumbent, and one that is only locally valid is undone when the subtree
 * is left. Reoptimization across compute() calls is described in
 * intReoptimize.
 *
 * \par References.
 * E. L. Lawler, D. E. Wood, "Branch-and-bound methods: a survey",
 * Operations Research 14(4):699-719, 1966. \n
 * D. R. Morrison, S. H. Jacobson, J. J. Sauppe, E. C. Sewell,
 * "Branch-and-bound algorithms: a survey of recent advances in searching,
 * branching, and pruning", Discrete Optimization 19:79-102, 2016. \n
 * M. Fischetti, A. Lodi, "Heuristics in mixed integer programming", in
 * Wiley Encyclopedia of Operations Research and Management Science, 2011. \n
 * T. Berthold, "Primal heuristics for mixed integer programs", Diploma
 * thesis, Technische Universitaet Berlin, 2006. \n
 * B. Gendron, T. G. Crainic, "Parallel branch-and-bound algorithms: survey
 * and synthesis", Operations Research 42(6):1042-1066, 1994. \n
 * D. A. Bader, W. E. Hart, C. A. Phillips, "Parallel algorithm design for
 * branch and bound", in H. J. Greenberg (ed.), Tutorials on Emerging
 * Methodologies and Applications in Operations Research, Springer, 2005. \n
 * B. Archibald, P. Maier, C. McCreesh, R. Stewart, P. Trinder, "Replicable
 * parallel branch and bound search", Journal of Parallel and Distributed
 * Computing 113:92-114, 2018. */

class BranchAndXSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 @{ */

 using Index = Block::Index;

 /// public enum for the int algorithmic parameters

 enum int_par_type_BXS {
  intSolveMethod = intLastAlgPar ,  ///< how to explore the tree
  intThreadForDifferentSolvers ,    ///< max threads for solvers at each node
  intBoundingProtocol ,             /**< when a node is evaluated:
   * with the eager protocol (the default) a node is evaluated as soon as it is
   * created, i.e., inside the branching of its parent, so that it is only put
   * in the open set if its relaxation bound can improve the incumbent; with
   * the lazy one it is put there unevaluated, carrying the relaxation bound of
   * its parent, and it is evaluated when it is extracted. Evaluating late can
   * save the whole evaluation of the nodes that a better incumbent, found in
   * the meantime, prunes without ever looking at them; evaluating early gives
   * the open set the true bound of each node, which is what a best-first
   * discipline orders by [see intSolveMethod]. Either way a node is evaluated
   * at most once [see Node::is_evaluated()]. */
  intReoptimize ,                   /**< retain the tree to reoptimize:
   * with it any serial exploration keeps the tree and its fenced frontier
   * alive across compute() calls, and a re-solve under changes that leave
   * the tree valid [see keeps_tree()] re-seeds the search from the
   * re-evaluated frontier instead of re-deriving the whole tree [see
   * reseedFrontier()]; the strategy only decides the order in which the
   * frontier is re-explored, which is that of its own open set. The
   * parallel depth-first exploration [see ParallelDFSSolve()] does not
   * retain anything, and always solves from scratch. Note
   * that in a binary enumeration the fenced frontier is about as large as
   * the interior, so this pays off when evaluating a node is expensive
   * (say, LP-like relaxations) and/or few nodes re-open, while for very
   * cheap relaxations re-evaluating the frontier may cost as much as
   * solving from scratch; it is off by default, also because retaining
   * the tree can use a lot of memory. */
  intMaxNodes ,                     ///< max number of nodes to explore
  intLastBXSPar                     ///< first allowed new int parameter
  };

 /// public enum for the string algorithmic parameters

 enum str_par_type_BXS {
  /// name of the BlockSolverConfig file providing the inner Solver
  strNameOfBlockSolverConfigurationFile = strLastAlgPar ,
  strLastBXSPar                    ///< first allowed new string parameter
  };

 /// public enum for when the inner Solver are run on a node

 enum BoundingProtocol {
  Eager = 0 ,                       ///< when the node is created
  Lazy = 1                          ///< when the node is extracted
  };

/*--------------------------------------------------------------------------*/
 /// public enum for the possible tree exploration strategies

 enum SolveMethod {
  DFS = 0 ,                         ///< depth-first search
  BFS = 1 ,                         ///< breadth-first search
  BestFS = 2 ,                      ///< best-first search
  BestFSDive = 3 ,                  ///< best-first search with dives
  BestEstimate = 4                  ///< best-estimate search
  };

/*--------------------------------------------------------------------------*/
 /// GlobalInformation names / keys defined by the BranchAndXSolver
 /** The names and keys, beyond the reserved ones [see GlobalInformation.h],
  * that the BranchAndXSolver defines in the GlobalInformation it hands to
  * its relaxations:
  *
  * - str_GlobalCuts ("GlobalCuts") names a
  *   Collection< std::vector< std::shared_ptr< Change > > > holding, under
  *   the single key str_Cuts ("Cuts"), the pool of the globally valid cuts
  *   (valid in every node of the tree) contributed and consulted by the
  *   relaxations: append-only (entries are only ever push_back()-ed, under
  *   write_with()), with shared ownership of the Change so that the same
  *   cut can be consumed by several solvers at once; a consumer remembers
  *   how many entries it has already seen and only reads the new tail
  *   (under read_with()). In perspective a symmetric pool serves the
  *   globally valid columns.
  *
  * The BranchAndXSolver does not create this Collection: it is a channel
  * between relaxations, which it neither fills nor reads, so it is created
  * (add_to_Universe()) by the first relaxation that contributes a cut. Only
  * the name and the key are fixed here, so that producers and consumers can
  * agree on them. */

 static constexpr const char * str_GlobalCuts = "GlobalCuts";

 /// key (in str_GlobalCuts) of the pool of the globally valid cuts
 static constexpr const char * str_Cuts = "Cuts";

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor: the inner Solver are provided later via the configuration

 BranchAndXSolver() : Solver() , f_RelaxationSolvers() ,
                      f_HeuristicSolvers() , f_state( kUnEval ) ,
                      bestBound( 0 ) , f_open_bound( - Inf< double >() ) ,
                      bestSolution( nullptr ) ,
                      solveType( BestFS ) , maxThreadForSolvers( 1 ) ,
                      boundingProtocol( Eager ) ,
                      maxThread( 0 ) , reoptimize( 0 ) ,
                      maxNodes( INT_MAX ) , f_treeRoot( nullptr ) ,
                      nodeBudget( INT_MAX ) , timeBudget( Inf< double >() ) ,
                      relTol( 0 ) , absTol( 0 ) ,
                      configurationRS( nullptr ) ,
                      changes( Modification::eModAll ) ,
                      subOptimalNodes() , infeasibleNodes() ,
                      integerNodes() ,
                      f_incumbentCell( nullptr ) , f_lfaCell( nullptr ) { }

/*--------------------------------------------------------------------------*/
 /// destructor

 ~BranchAndXSolver() override {
  discardRetainedTree();
  delete bestSolution;
  delete configurationRS;
  for( auto slvr : v_created ) {
   slvr->set_Block( nullptr );
   delete slvr;
   }
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve
 /** Attaching to a Block different from the current one (nullptr, i.e.,
  * detaching, included) drops everything that belongs to the old one: the
  * inner Solver created out of the BlockSolverConfig, the tree retained for
  * reoptimization and the best solution found; the inner Solver are then
  * created anew on the new Block, if a BlockSolverConfig is there, and the
  * next compute() solves from scratch. */

 void set_Block( Block * block ) override {
  if( block == f_Block )
   return;
  if( f_Block )
   dropInnerSolvers();
  Solver::set_Block( block );
  f_state = kUnEval;
  changes = Modification::eModAll;
  if( block ) {
   initializeGlobalInformation();
   if( configurationRS )
    applyConfigurationToSolvers();
   }
  }

/*--------------------------------------------------------------------------*/

 using Solver::set_par;  // restore the hidden overloaded methods

/*--------------------------------------------------------------------------*/
 /// set the int parameters of BranchAndXSolver
 /** Set the int parameters specific of BranchAndXSolver:
  *
  * - intSolveMethod [BestFS]: how to explore the tree (see SolveMethod);
  *   BestEstimate extracts first the node whose estimate of the best value
  *   under it is the best, the estimate being the one that the
  *   RelaxationSolver that branched gives of the child, if it derives from
  *   BranchEstimator [see BranchingRule.h], and the bound of the parent
  *   otherwise, and never better than the bound of the node; the nodes are
  *   still pruned by their bounds, so that the estimates only change the
  *   order of the exploration.
  *
  * - intBoundingProtocol [Eager]: when a node is evaluated (see
  *   BoundingProtocol); it changes nothing of what is explored, only when
  *   the inner Solver are run on a node.
  *
  * - intThreadForDifferentSolvers [1]: maximum number of threads used to
  *   evaluate the (multiple) inner Solver at each node in parallel [see
  *   computeRelaxationsParallel() / computeHeuristicParallel()]; 1 = serial.
  *   It only pays when several (expensive) relaxations / heuristics are
  *   attached: their compute() is run concurrently and the results reduced
  *   exactly as in the serial path. Active in all the (serial-tree)
  *   explorations: depth-, breadth- and best-first.
  *
  * - intMaxNodes [no limit]: the budget of tree nodes a solve may explore;
  *   <= 0 means no limit. It is the enumerative counterpart of the inherited
  *   intMaxIter (a node is not an iteration: a node may take several inner
  *   iterations), which is left to the inner Solvers.
  *
  * The inherited intMaxThread [0] sets the number of workers of the
  * parallel tree exploration [see ParallelDFSSolve(); <= 1 = serial; the
  * base Solver classes do not store it, so it is stored here]. The inherited
  * dblMaxTime / dblRelAcc / dblAbsAcc respectively bound the total solve time
  * and set the optimality tolerances of the pruning. */

 void set_par( idx_type par , int value ) override {
  switch( par ) {
   case( intSolveMethod ):
    if( ( value < DFS ) || ( value > BestEstimate ) )
     throw( std::invalid_argument( "BranchAndXSolver::set_par: "
            "invalid intSolveMethod" ) );
    solveType = static_cast< SolveMethod >( value );
    break;
   case( intMaxThread ):
    maxThread = std::max( 0 , value );
    break;
   case( intThreadForDifferentSolvers ):
    if( value <= 0 )
     throw( std::invalid_argument( "BranchAndXSolver::set_par: "
            "intThreadForDifferentSolvers must be positive" ) );
    maxThreadForSolvers = value;
    break;
   case( intBoundingProtocol ):
    if( ( value < Eager ) || ( value > Lazy ) )
     throw( std::invalid_argument( "BranchAndXSolver::set_par: "
            "invalid intBoundingProtocol" ) );
    boundingProtocol = static_cast< BoundingProtocol >( value );
    break;
   case( intReoptimize ):
    reoptimize = value;
    break;
   case( intMaxNodes ):
    maxNodes = ( value <= 0 ) ? INT_MAX : value;
    break;
   default:
    Solver::set_par( par , value );
   }
  }

/*--------------------------------------------------------------------------*/
 /// set the double parameters of BranchAndXSolver
 /** The base Solver classes do not store the double parameters, so the
  * inherited ones that BranchAndXSolver uses are stored here:
  *
  * - dblMaxTime [+INF]: the time budget of a solve, in seconds;
  *
  * - dblRelAcc [1e-6] / dblAbsAcc [+INF, i.e., not active]: the relative and
  *   absolute tolerances of the pruning, a node being pruned when its
  *   relaxation bound cannot improve the incumbent by more than the largest
  *   of the two (the relative one taken w.r.t. max( | incumbent | , 1 )). */

 void set_par( idx_type par , double value ) override {
  switch( par ) {
   case( dblMaxTime ):
    maxTime = value;
    break;
   case( dblRelAcc ):
    relAcc = value;
    break;
   case( dblAbsAcc ):
    absAcc = value;
    break;
   default:
    Solver::set_par( par , value );
   }
  }

/*--------------------------------------------------------------------------*/
 /// set the string parameters of BranchAndXSolver
 /** Set the string parameters specific of BranchAndXSolver:
  *
  * - strNameOfBlockSolverConfigurationFile [""]: name of the file out of
  *   which the BlockSolverConfig providing the inner Solver (the
  *   RelaxationSolver(s) and heuristic ChangeSolver(s)) is deserialized;
  *   it is applied to the Block as soon as both are available, the inner
  *   Solver it describes replacing those of any previous one. The empty
  *   name (the default) means no BlockSolverConfig: the inner Solver
  *   already there, if any, are kept.
  *
  * The string parameters are taken by the rvalue version, which the one
  * taking a const reference forwards to [see ThinComputeInterface]. */

 void set_par( idx_type par , std::string && value ) override {
  if( par != strNameOfBlockSolverConfigurationFile ) {
   Solver::set_par( par , std::move( value ) );
   return;
   }
  BlockSolverConfig * cfg = nullptr;
  if( ! value.empty() ) {
   auto c = Configuration::deserialize( value );
   cfg = dynamic_cast< BlockSolverConfig * >( c );
   if( ! cfg ) {
    delete c;
    throw( std::invalid_argument( "BranchAndXSolver::set_par: unable to "
           "create a BlockSolverConfig out of " + value ) );
    }
   }
  nameRS = value;
  delete configurationRS;
  configurationRS = cfg;
  if( f_Block && configurationRS ) {
   dropInnerSolvers();
   changes = Modification::eModAll;
   applyConfigurationToSolvers();
   }
  }

/*--------------------------------------------------------------------------*/
 /// receive a Modification and forward it to the (private) inner Solver

 void add_Modification( sp_Mod & mod ) override {
  Solver::add_Modification( mod );
  for( auto slvr : v_created )
   slvr->add_Modification( mod );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

 /// solve the problem encoded in the Block with Branch-and-Bound

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// write the best solution found into the Variable of the Block

 void get_var_solution( Configuration * solc = nullptr ) override {
  if( bestSolution )
   bestSolution->write( f_Block );
  }

/*--------------------------------------------------------------------------*/
 /// return the value of the best solution found

 OFValue get_var_value( void ) override { return( bestBound ); }

/*--------------------------------------------------------------------------*/
 /// return a lower bound on the optimal objective function value
 /** For a maximization problem the incumbent is a valid lower bound; for a
  * minimization one it is the global relaxation bound of the exploration,
  * i.e., the smallest among the relaxation bounds of the nodes left open and
  * the incumbent, which is the incumbent itself once the enumeration has been
  * completed; a node not evaluated yet has the relaxation bound of its father,
  * hence the bound is finite as soon as the root has been. An exploration
  * stopped by a budget with the parallel depth-first search [see
  * intSolveMethod] claims no bound on the side of the open nodes. */

 OFValue get_lb( void ) override;

/*--------------------------------------------------------------------------*/
 /// return an upper bound on the optimal objective function value
 /** The symmetric of get_lb(): the incumbent for a minimization problem, the
  * global relaxation bound (the largest among the relaxation bounds of the
  * nodes left open and the incumbent) for a maximization one. */

 OFValue get_ub( void ) override;

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the BranchAndXSolver
 *  @{ */

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastBXSPar ) );
  }

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return( idx_type( strLastBXSPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  switch( par ) {
   case( intSolveMethod ):               return( int( BestFS ) );
   case( intThreadForDifferentSolvers ): return( 1 );
   case( intBoundingProtocol ):          return( int( Eager ) );
   case( intReoptimize ):                return( 0 );
   case( intMaxNodes ):                  return( INT_MAX );
   }
  return( Solver::get_dflt_int_par( par ) );
  }

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override {
  if( par == strNameOfBlockSolverConfigurationFile ) {
   static const std::string dfltRS;
   return( dfltRS );
   }
  return( Solver::get_dflt_str_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  switch( par ) {
   case( intSolveMethod ):               return( int( solveType ) );
   case( intMaxThread ):                 return( maxThread );
   case( intThreadForDifferentSolvers ): return( maxThreadForSolvers );
   case( intBoundingProtocol ):          return( int( boundingProtocol ) );
   case( intReoptimize ):                return( reoptimize );
   case( intMaxNodes ):                  return( maxNodes );
   }
  return( Solver::get_int_par( par ) );
  }

 [[nodiscard]] double get_dbl_par( idx_type par ) const override {
  switch( par ) {
   case( dblMaxTime ): return( maxTime );
   case( dblRelAcc ):  return( relAcc );
   case( dblAbsAcc ):  return( absAcc );
   }
  return( Solver::get_dbl_par( par ) );
  }

 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override {
  if( par == strNameOfBlockSolverConfigurationFile )
   return( nameRS );
  return( Solver::get_str_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "intSolveMethod" )
   return( intSolveMethod );
  if( name == "intThreadForDifferentSolvers" )
   return( intThreadForDifferentSolvers );
  if( name == "intBoundingProtocol" )
   return( intBoundingProtocol );
  if( name == "intReoptimize" )
   return( intReoptimize );
  if( name == "intMaxNodes" )
   return( intMaxNodes );
  return( Solver::int_par_str2idx( name ) );
  }

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override {
  if( name == "strNameOfBlockSolverConfigurationFile" )
   return( strNameOfBlockSolverConfigurationFile );
  return( Solver::str_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::string pars[] = { "intSolveMethod" ,
                                      "intThreadForDifferentSolvers" ,
                                      "intBoundingProtocol" ,
                                      "intReoptimize" , "intMaxNodes" };
  if( ( idx >= intSolveMethod ) && ( idx < intLastBXSPar ) )
   return( pars[ idx - intSolveMethod ] );
  return( Solver::int_par_idx2str( idx ) );
  }

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override {
  if( idx == strNameOfBlockSolverConfigurationFile ) {
   static const std::string par = "strNameOfBlockSolverConfigurationFile";
   return( par );
   }
  return( Solver::str_par_idx2str( idx ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// setter for the tree exploration strategy

 void setSolveType( SolveMethod st ) { solveType = st; }

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// the inner solvers of the relaxations: ( trait , Solver personality )
 /** Each inner solver is one object with two personalities [see
  * ChangeSolver.h]: the RelaxationSolver trait (branch(), apply(),
  * get_true_*()) and the Solver proper (compute(), get_lb() / get_ub()).
  * Since neither derives from the other, both pointers to the same object
  * are kept, obtained by cross-cast at creation. */
 std::vector< std::pair< RelaxationSolver * , Solver * > >
  f_RelaxationSolvers;

 /// the inner heuristic solvers: ( trait , Solver personality )
 std::vector< std::pair< ChangeSolver * , Solver * > > f_HeuristicSolvers;

 int f_state;             ///< the (current) state of the compute() process

 double bestBound;        ///< the best primal (hence feasible) bound found
 double f_open_bound;     ///< global relaxation bound of the last compute()

 Solution * bestSolution; ///< the best solution found, valued bestBound

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// create the inner Solver out of configurationRS
 /** Creates and ComputeConfig-ures the Solver described by the
  * BlockSolverConfig and attaches them to the Block, WITHOUT registering
  * them: the inner Solver are private to the BranchAndXSolver, which
  * owns them, drives them, and forwards them the Modification it receives
  * [see add_Modification()]; this way they never appear in the Block's
  * registered-Solver list (whose order and content belong to the user). */

 /// declare in the GlobalInformation what this search shares
 /** Declares in the GlobalInformation the Collection that the search shares
  * with its inner Solver: the reserved "hot" scalars and flags of the whole
  * cooperation [see GlobalInformation.h], i.e., the incumbent and the
  * local-fixing-allowed flag. The references to the two individual atomics
  * are cached once here, and used lock-free ever after; the incumbent starts
  * at "no incumbent", which is why this needs the Block: which infinity that
  * is depends on the sense of the Objective. Called by set_Block(), i.e.,
  * before any inner Solver is handed the GlobalInformation. */

 void initializeGlobalInformation( void ) {
  if( ! f_incumbentCell ) {   // the first Block: declare the Collection
   f_globalInfo.add_to_Universe< std::atomic< double > >(
				GlobalInformation::str_AtomicScalars );
   f_incumbentCell = &( ( * f_globalInfo.get_from_Universe<
			  std::atomic< double > >(
			   GlobalInformation::str_AtomicScalars ) )[
			    GlobalInformation::str_Incumbent ] );

   f_globalInfo.add_to_Universe< std::atomic< bool > >(
				GlobalInformation::str_AtomicFlags );
   f_lfaCell = &( ( * f_globalInfo.get_from_Universe< std::atomic< bool > >(
		     GlobalInformation::str_AtomicFlags ) )[
		      GlobalInformation::str_LocalFixingAllowed ] );
   }
  f_incumbentCell->store( no_incumbent() );
  f_lfaCell->store( true );
  }

/*--------------------------------------------------------------------------*/
 /// drop what belongs to the current Block, before attaching to another
 /** Deletes the inner Solver, the tree retained for reoptimization and the
  * best Solution, all of which refer to the current Block, and the pool of
  * the globally valid cuts that the relaxations may have created in the
  * GlobalInformation [see str_GlobalCuts]. */

 void dropInnerSolvers( void ) {
  discardRetainedTree();
  delete bestSolution;
  bestSolution = nullptr;
  bestBound = 0;
  for( auto slvr : v_created ) {
   slvr->set_Block( nullptr );
   delete slvr;
   }
  v_created.clear();
  f_RelaxationSolvers.clear();
  f_HeuristicSolvers.clear();
  v_workerSolvers.clear();
  f_globalInfo.remove_from_Universe( str_GlobalCuts );
  }

/*--------------------------------------------------------------------------*/
 /// the value of the incumbent meaning that there is none
 /** The value that the incumbent has while no feasible solution has been
  * found yet: the infinity that no solution can improve upon, i.e., +INF
  * for a minimization problem and -INF for a maximization one. Any bound is
  * therefore "better" than it without any special case, and a consumer that
  * rather wants to check whether there is an incumbent at all just tests it
  * for finiteness. */

 [[nodiscard]] double no_incumbent( void ) const;

/*--------------------------------------------------------------------------*/
 /// create the per-worker private Solver sets for the parallel exploration
 /** Creates \p K additional sets of inner Solver out of configurationRS,
  * one per worker of the parallel exploration (see v_workerSolvers): each
  * worker drives its own clones, so that the (cheap, internal) application
  * of the branching Changes never needs synchronization; the incumbent and
  * the pool of open subtrees are the only shared state. The clones are
  * attached to the Block (set_Block, no registration) and owned. */

 void createWorkerSolvers( Index K ) {
  if( ! ( configurationRS && f_Block ) )
   throw( std::logic_error( "BranchAndXSolver::createWorkerSolvers: "
          "configurationRS or f_Block not set" ) );
  v_workerSolvers.resize( K );
  for( Index w = 0 ; w < K ; ++w ) {
   auto & ws = v_workerSolvers[ w ];
   if( ! ws.relaxation.empty() )    // already created (re-solve)
    continue;
   for( Index i = 0 ; i < configurationRS->num_ComputeConfig() ; ++i ) {
    auto slvr = Solver::new_Solver( configurationRS->get_SolverName( i ) );
    if( auto cfg = configurationRS->get_SolverConfig( i ) )
     slvr->set_ComputeConfig( cfg );
    slvr->set_Block( f_Block );
    v_created.push_back( slvr );    // owned like the serial ones
    // the traits are discovered on the Solver by cross-cast [see
    // ChangeSolver.h] and kept alongside the Solver personality
    if( auto rs = dynamic_cast< RelaxationSolver * >( slvr ) ) {
     rs->set_global_information( &f_globalInfo );
     ws.relaxation.push_back( { rs , slvr } );
     }
    else if( auto hs = dynamic_cast< ChangeSolver * >( slvr ) )
     ws.heuristic.push_back( { hs , slvr } );
    else
     throw( std::invalid_argument( "BranchAndXSolver::"
            "createWorkerSolvers: the BlockSolverConfig must only "
            "describe :ChangeSolver" ) );
    }
   }
  }

/*--------------------------------------------------------------------------*/

 void applyConfigurationToSolvers( void ) {
  if( ! ( configurationRS && f_Block ) )
   throw( std::logic_error( "BranchAndXSolver::"
          "applyConfigurationToSolvers: configurationRS or f_Block not set"
          ) );
  for( Index i = 0 ; i < configurationRS->num_ComputeConfig() ; ++i ) {
   auto slvr = Solver::new_Solver( configurationRS->get_SolverName( i ) );
   if( auto cfg = configurationRS->get_SolverConfig( i ) )
    slvr->set_ComputeConfig( cfg );
   slvr->set_Block( f_Block );
   v_created.push_back( slvr );
   // the traits are discovered on the Solver by cross-cast [see
   // ChangeSolver.h] and kept alongside the Solver personality
   if( auto rs = dynamic_cast< RelaxationSolver * >( slvr ) ) {
    rs->set_global_information( &f_globalInfo );
    f_RelaxationSolvers.push_back( { rs , slvr } );
    }
   else if( auto hs = dynamic_cast< ChangeSolver * >( slvr ) )
    f_HeuristicSolvers.push_back( { hs , slvr } );
   else
    throw( std::invalid_argument( "BranchAndXSolver::"
           "applyConfigurationToSolvers: the BlockSolverConfig must only "
           "describe :ChangeSolver" ) );
   }
  }

/*--------------------------------------------------------------------------*/
 /// parallel depth-first exploration with \p K workers
 /** Parallel depth-first exploration, used when the inherited intMaxThread
  * parameter is > 1: a serial, ordered (FIFO, i.e., discovery order)
  * ramp-up expands the tree BestFS-style until enough open subtrees exist,
  * then \p K workers, each driving its own private set of inner Solver
  * [see createWorkerSolvers()], repeatedly claim the OLDEST open subtree
  * (preserving the sequential search order, which is what keeps parallel
  * performance replicable) and explore it depth-first, sharing only the
  * incumbent (under mutex) and the pool of open subtrees.
  *
  * The distribution of work is dynamic: when the pool runs low a worker
  * DONATES one of the (eagerly evaluated, hence claimable) top-level
  * branches of its current subtree instead of exploring it, so that
  * deep-and-narrow trees keep all the workers busy; donated children hang
  * off the claimed node, which always outlives them, keeping their path
  * replayable and the final cleanup single-owner. Plain std::thread is
  * used on purpose: at this granularity (one whole subtree per claim) a
  * task framework would add nothing. */

 int ParallelDFSSolve( std::mutex & globalMutex , int K );

/*--------------------------------------------------------------------------*/
 /// delete the tree retained for reoptimization (if any) and its frontier

 void discardRetainedTree( void );

/*--------------------------------------------------------------------------*/
 /// the depth-first subtree exploration of one worker
 /** The recursive core of one worker of ParallelDFSSolve(): a depth-first
  * exploration on the worker's own Solver set, with the shared
  * incumbent updated under \p incumbentMutex, a shared atomic node budget
  * and a wall-clock deadline. */

 int workerDFS( Node * currentNode , std::list< ChangeSolver * > & solvers ,
                std::vector< std::pair< RelaxationSolver * ,
                             Solver * > > & relaxation ,
                std::vector< std::pair< ChangeSolver * ,
                             Solver * > > & heuristic ,
                bool minimizing , std::mutex & incumbentMutex ,
                std::atomic< int > & nodeBdg ,
                std::chrono::steady_clock::time_point deadline ,
                int & nameCounter );

/*--------------------------------------------------------------------------*/
 /// solve the tree with the exploration strategy that is set
 /** The one serial exploration: it builds the OpenList that the chosen
  * strategy asks for [see intSolveMethod and OpenList] - a LIFO stack for
  * depth-first, a FIFO queue for breadth-first, a relaxation-bound priority
  * queue for best-first, its diving variant for BestFSDive - seeds it with the
  * root (or with the frontier of the retained tree, see reseedFrontier()) and
  * hands it to the shared loop [see explore()]. The strategies differ in
  * nothing else.
  *  @param globalMutex mutex protecting the shared state of this class
  *  @return the sol_type [see Solver.h] of the computation */

 int treeSolve( std::mutex & globalMutex );

/*--------------------------------------------------------------------------*/
 /// the single serial tree exploration, driven by the open-list discipline
 /** The one exploration loop shared by every (serial) strategy: it
  * repeatedly takes the next node from \p open [see OpenList], moves the
  * :ChangeSolver to it, evaluates and branches it, and pushes the surviving
  * children back into \p open; the discipline of \p open (stack / queue /
  * priority) is what makes it a depth- / breadth- / best-first search.
  * \p currentNode is where the :ChangeSolver currently sits (the root, or
  * the last re-seeded node in a reoptimization), \p retain keeps the tree
  * and its fenced frontier alive for a later reoptimization [see
  * intReoptimize]. */

 int explore( OpenList & open , std::mutex & globalMutex ,
              std::list< ChangeSolver * > * solvers , bool minimizing ,
              int & counter , ExploringNode * rootNode ,
              ExploringNode * currentNode , bool retain ,
              std::chrono::high_resolution_clock::time_point start );

/*--------------------------------------------------------------------------*/
 /// re-seed the exploration from the frontier of the retained tree
 /** Reoptimization [see intReoptimize]: instead of building a new tree, the
  * one retained by the previous solve is taken back (\p rootNode is set to
  * it) and its fenced frontier is re-evaluated under the new data, each node
  * being either re-fenced or re-opened into \p open; the interior of the
  * tree is NOT re-derived. When the changes cannot make the feasible
  * region grow [see Modification::may_grow_region()], e.g., when only the
  * Objective changes, no infeasible node can be un-fenced, and that part of
  * the frontier is left untouched.
  *
  * The re-opened nodes are handed to \p open respecting its discipline (a
  * LIFO stack receives them in reverse), exactly as explore() does with the
  * children of a node: the frontier is sorted by the relaxation bound of the
  * previous solve, so that whatever the strategy the most promising nodes
  * are looked at first and the incumbent warms up immediately. This is why
  * the mechanism is independent of the exploration strategy: the order in
  * which the frontier is re-explored is entirely a matter of \p open.
  *
  * \p currentNode is left where the :ChangeSolver actually sit, which on
  * success is the last re-evaluated node and on failure is \p rootNode.
  *  @return the sol_type [see Solver.h] of the re-evaluation */

 int reseedFrontier( OpenList & open , std::list< ChangeSolver * > * solvers ,
                     bool minimizing , ExploringNode * & rootNode ,
                     ExploringNode * & currentNode );

/*--------------------------------------------------------------------------*/
 /// tear down a retained / in-progress tree and its open set and frontier
 /** Besides, the global relaxation bound is set to the best among those of the
  * nodes left in the open set and of \p current, the node whose evaluation
  * (or that of a child of which) has stopped the exploration, if any: this
  * is what get_lb() / get_ub() return when the stop is due to a budget. */

 void discardTree( ExploringNode * root , OpenList & open ,
                   std::list< ChangeSolver * > * solvers ,
                   const ExploringNode * current = nullptr );

/*--------------------------------------------------------------------------*/
 /// process the outstanding Modification
 /** Processes the queued Modification and accumulates what they change
  * [see Modification::changes()] in the changes field, which tells a
  * re-solve whether the retained tree can be reused [see keeps_tree()] and
  * which part of its frontier is to be re-evaluated [see
  * reseedFrontier()]. */

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
 /// true if the retained tree is still valid under the given changes
 /** True if a tree retained across compute() calls can be reused under the
  * changes \p c, which are those of a single Modification [see
  * process_outstanding_Modification()]. It cannot when the data of some
  * Variable change, when Variable or Constraint are removed or when some
  * Variable becomes continuous, since the branching Changes of the tree,
  * which act on the Variable and on the Constraint, may then undo the
  * change, act on something that is no longer there or no longer cover the
  * region (x <= 0 and x >= 1 leave out ( 0 , 1 )). It can under anything
  * else, which is a relaxation (added Variable, enlarged sides) or a
  * restriction (added Constraint, Variable becoming integer, shrunk sides)
  * of the problem that the re-evaluation of the fenced frontier takes care
  * of [see reseedFrontier()]: the infeasible nodes are re-evaluated only
  * when the region may grow. A physical Modification that does not say
  * what it changes of the model counts as changing anything. */

 static bool keeps_tree( Modification::ModConcern c );

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SolveMethod solveType;       ///< chosen method to explore the tree

 int maxThreadForSolvers;     ///< max threads for solvers at each node

 BoundingProtocol boundingProtocol;  ///< when a node is evaluated (see
                                     ///< intBoundingProtocol)

 int maxThread;               ///< workers of the parallel tree exploration

 int reoptimize;              ///< retain the tree to reoptimize (see
                              ///< intReoptimize / treeSolve())

 int maxNodes;                ///< node budget of a solve (see intMaxNodes)

 /// the inherited dblMaxTime [see set_par( double )]
 double maxTime = Inf< double >();

 /// the inherited dblRelAcc [see set_par( double )]
 double relAcc = 1e-6;

 /// the inherited dblAbsAcc [see set_par( double )]
 double absAcc = Inf< double >();

 /// the global information shared with the relaxations
 /** The GlobalInformation handed to every relaxation: the reserved hot
  * scalars / flags (incumbent, local-fixing-allowed), declared by
  * initializeGlobalInformation(), plus whatever the relaxations put there
  * themselves (say, the pool of the globally valid cuts, see
  * str_GlobalCuts). The individual atomics are reached through the
  * references cached below. */
 GlobalInformation f_globalInfo;

 /// the cached reference to the atomic incumbent cell
 /** Mirrors bestBound: every improvement of the incumbent is store()-d
  * here too, which is how the relaxations see it lock-free [see
  * GlobalInformation::str_Incumbent]; not finite when no feasible solution
  * has been found yet [see no_incumbent()]. nullptr until set_Block(). */
 std::atomic< double > * f_incumbentCell;

 /// the cached reference to the atomic local-fixing-allowed flag
 /** Set by each compute(), which is where it is known whether the tree is
  * being retained across re-solves [see intReoptimize]. nullptr until
  * set_Block(). */
 std::atomic< bool > * f_lfaCell;

 /// the root of the tree retained for reoptimization, nullptr if none
 ExploringNode * f_treeRoot;

 /// residual node budget of the current solve (from intMaxNodes), consumed
 /// by the exploration
 int nodeBudget;

 /// residual time budget of the current solve (from the inherited
 /// dblMaxTime), consumed by the exploration
 double timeBudget;

 double relTol;   ///< the inherited dblRelAcc, cached at compute() start

 double absTol;   ///< the inherited dblAbsAcc, cached at compute() start

 std::string nameRS;          ///< name of the BlockSolverConfig file

 BlockSolverConfig * configurationRS;  ///< the inner BlockSolverConfig

 /// the (private) inner Solver created out of configurationRS, owned
 std::vector< Solver * > v_created;

 /// the private Solver set of one worker of the parallel exploration
 struct WorkerSolvers {
  std::vector< std::pair< RelaxationSolver * , Solver * > > relaxation;
  std::vector< std::pair< ChangeSolver * , Solver * > > heuristic;
  };

 /// the per-worker Solver sets [see createWorkerSolvers()]
 std::vector< WorkerSolvers > v_workerSolvers;

 /// what the outstanding Modification change [see Modification::changes()],
 /// 0 if there are none
 Modification::ModConcern changes;

 /// nodes pruned as sub-optimal, to be revisited on an objective change
 std::list< ExploringNode * > subOptimalNodes;

 /// nodes pruned as infeasible, to be revisited on a r.h.s. change
 std::list< ExploringNode * > infeasibleNodes;

 /// nodes pruned as integer-feasible, to be revisited on reoptimization
 std::list< ExploringNode * > integerNodes;

 SMSpp_insert_in_factory_h;   // insert BranchAndXSolver in the factory

/*--------------------------------------------------------------------------*/

 };  // end( class( BranchAndXSolver ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS Node --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a node of the Branch-and-Bound tree
/** A Node of the enumeration tree is identified by the Change that leads to
 * it from its parent; it also holds the (undo) Change to return to the
 * parent, the Changes generating its children (see obtainBranchList()) and
 * its relaxation bound. The Node owns all these Change. */

class Node {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
 /// constructor
 /** @param change the Change to apply to reach this node from its parent
  *  @param nodeName name of the node (debugging purposes only) */

 Node( Change * change , int nodeName = 0 )
  : bound( - Inf< double >() ) ,
    estimate( std::numeric_limits< double >::quiet_NaN() ) ,
    f_change( change ) , toFather( nullptr ) , branches() ,
    f_infeasible( false ) , f_evaluated( false ) , name( nodeName ) {}

/*--------------------------------------------------------------------------*/
 /// destructor: deletes the owned Changes

 virtual ~Node() {
  delete f_change;
  delete toFather;
  for( auto ch : branches )
   delete ch;
  }

/*--------------------------------------------------------------------------*/

 /// return the value of the relaxation bound
 double get_bound( void ) const { return( bound ); }

 /// set the value of the relaxation bound
 void set_bound( double db ) { bound = db; }

 /// return the estimate of the best value under the node, NaN if none
 /** The estimate that the RelaxationSolver that branched the parent gave of
  * this node [see BranchEstimator]. */
 double get_estimate( void ) const { return( estimate ); }

 /// set the estimate of the best value under the node
 void set_estimate( double e ) { estimate = e; }

 /// the key of the node in a best-estimate search
 /** Its estimate, but never better than its bound (the bound itself if it
  * has no estimate), in the sense \p minimizing. */
 double estimate_key( bool minimizing ) const {
  if( std::isnan( estimate ) )
   return( bound );
  return( minimizing ? std::max( estimate , bound )
                     : std::min( estimate , bound ) );
  }

 /// return the Change that leads to this node
 Change * get_f_change( void ) const { return( f_change ); }

 /// return the Change to undo to return to the parent node
 Change * get_toFather( void ) const { return( toFather ); }

 /// set the Change to undo to return to the parent node
 void set_toFather( Change * tf ) { toFather = tf; }

 /// return the name of the node (debugging purposes only)
 int get_name( void ) const { return( name ); }

 /// initialize the relaxation bound to the appropriate infinity
 void initializeBound( bool minimizing ) {
  bound = minimizing ? - Inf< double >() : Inf< double >();
  f_infeasible = false;
  f_evaluated = false;
  }

 /// whether the node was found infeasible (its relaxation has no solution)
 bool is_infeasible( void ) const { return( f_infeasible ); }

 /// record that the node is infeasible
 void set_infeasible( bool i ) { f_infeasible = i; }

 /// whether the inner Solver have already been run on this node
 /** Tells whether the node has already been evaluated, i.e., whether its
  * relaxation bound and branching Changes are those of the inner Solver rather
  * than the ones it inherited when it was created; this is what makes the
  * lazy bounding protocol [see BranchAndXSolver::intBoundingProtocol] never
  * evaluate a node twice. */

 bool is_evaluated( void ) const { return( f_evaluated ); }

 /// record that the inner Solver have been run on this node
 void set_evaluated( bool e ) { f_evaluated = e; }

 /// have the given RelaxationSolver produce the branching Changes
 /** The Changes produced by a previous evaluation, those not yet handed to
  * a child (whose entries are nullptr), are deleted; the estimates of the
  * children are taken too, if the RelaxationSolver gives them [see
  * BranchEstimator], and kept only if there is one per Change. */
 void obtainBranchList( RelaxationSolver * solver ) {
  for( auto ch : branches )
   delete ch;
  branches = solver->branch();
  estimates.clear();
  if( auto be = dynamic_cast< BranchEstimator * >( solver ) ) {
   estimates = be->branch_estimates();
   if( estimates.size() != branches.size() )
    estimates.clear();
   }
  }

 /// return the Changes generating the children of this node
 std::vector< Change * > & getBranches( void ) { return( branches ); }

 /// forget the estimates of the children, once they have been handed to
 /// them or when the branching Changes no longer match them
 void clear_estimates( void ) { estimates.clear(); }

 /// the estimate of the child generated by the i-th branching Change, NaN
 /// if none [see obtainBranchList()]
 double child_estimate( std::size_t i ) const {
  return( i < estimates.size() ? estimates[ i ]
                               : std::numeric_limits< double >::quiet_NaN() );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 double bound;     ///< value of the relaxation bound at this node

 double estimate;  ///< estimate of the best value under it, NaN if none

 Change * f_change;     ///< the Change to reach this node from the parent

 Change * toFather;     ///< the (undo) Change to return to the parent

 /// the Changes to reach each child from this node
 std::vector< Change * > branches;

 /// the estimates of the children, one per Change of branches, or none
 std::vector< double > estimates;

 bool f_infeasible;     ///< whether the node's relaxation has no solution

 bool f_evaluated;      ///< whether the inner Solver have been run here

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 int name;              ///< name of the node, debugging purposes only

/*--------------------------------------------------------------------------*/

 };  // end( class( Node ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS DFSNode ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Node used in the depth-first exploration of the tree
/** In the recursive depth-first exploration no explicit tree structure is
 * needed (the recursion stack is the path to the root), so a DFSNode is
 * just a plain Node. */

class DFSNode : public Node {

 public:

 /// constructor, see Node::Node()

 DFSNode( Change * change , int nodeName = 0 ) : Node( change , nodeName ) {}

 ~DFSNode() override = default;   ///< destructor

 };  // end( class( DFSNode ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ExploringNode ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the Node used in the breadth-first / best-first explorations
/** An ExploringNode also knows its parent and children, so that the
 * :ChangeSolver can be moved between two arbitrary nodes of the tree (see
 * moveBetweenNodes()) by climbing up and down along the joining path. */

class ExploringNode : public Node {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
 /// constructor
 /** @param change the Change to apply to reach this node from its parent
  *  @param parent pointer to the parent node
  *  @param depth depth of the node in the tree
  *  @param nodeName name of the node (debugging purposes only) */

 ExploringNode( Change * change , ExploringNode * parent , int depth ,
                int nodeName = 0 )
  : Node( change , nodeName ) , f_children() , level( depth ) ,
    parent( parent ) {}

 ~ExploringNode() override = default;   ///< destructor

/*--------------------------------------------------------------------------*/

 /// getter for the level of the node
 int get_level( void ) const { return( level ); }

 /// getter for the parent node
 ExploringNode * get_parent( void ) const { return( parent ); }

 /// setter for the parent node
 void set_parent( ExploringNode * p ) { parent = p; }

 /// getter for the children nodes
 std::vector< ExploringNode * > & get_children( void ) {
  return( f_children );
  }

/*--------------------------------------------------------------------------*/
 /// move the given :ChangeSolver between two nodes of the tree
 /** Moves the given :ChangeSolver from \p sourceNode to \p destNode by
  * applying, along the path joining them, the undo Changes (climbing up
  * from the source to the common ancestor) and then, top-down, the node
  * Changes (descending from the common ancestor to the destination). The
  * nodes to descend are collected bottom-up in a buffer of the thread, which
  * is taken for the call, so that a move costs no allocation. */

 static void moveBetweenNodes( ExploringNode * sourceNode ,
                               ExploringNode * destNode ,
                               std::list< ChangeSolver * > * solvers ) {
  thread_local std::vector< ExploringNode * > buffer;
  std::vector< ExploringNode * > toDescend;
  toDescend.swap( buffer );
  toDescend.clear();
  while( sourceNode->level != destNode->level ) {
   if( sourceNode->level > destNode->level ) {
    for( const auto s : *solvers )
     s->apply( sourceNode->toFather , false );
    sourceNode = sourceNode->parent;
    }
   else {
    toDescend.push_back( destNode );
    destNode = destNode->parent;
    }
   }
  while( sourceNode != destNode ) {
   for( const auto s : *solvers )
    s->apply( sourceNode->toFather , false );
   sourceNode = sourceNode->parent;
   toDescend.push_back( destNode );
   destNode = destNode->parent;
   }
  // descending, the first solver produces the undo Change of any node that
  // does not have one yet, which is how a node entered for the first time
  // gets it [see BranchAndXSolver::intBoundingProtocol: a node evaluated
  // lazily is created without the solvers ever descending into it]
  for( auto it = toDescend.rbegin() ; it != toDescend.rend() ; ++it ) {
   auto node = *it;
   bool modified = false;
   for( const auto s : *solvers ) {
    if( ( ! modified ) && ( ! node->toFather ) )
     node->set_toFather( s->apply( node->f_change , true ) );
    else
     s->apply( node->f_change , false );
    modified = true;
    }
   }
  toDescend.swap( buffer );
  }

/*--------------------------------------------------------------------------*/
 /// prune the node, recursively pruning any parent left childless
 /** Prunes \p nodeToPrune, moving the :ChangeSolver to its parent; if the
  * parent is left with no children it is recursively pruned as well.
  *  @return the node the :ChangeSolver is left at, i.e., the parent of the
  *          last pruned node */

 static ExploringNode * prune( ExploringNode * nodeToPrune ,
                               std::list< ChangeSolver * > * solvers ) {
  ExploringNode * parentNode = nodeToPrune->get_parent();
  if( nodeToPrune->get_toFather() )
   for( const auto s : *solvers )
    s->apply( nodeToPrune->get_toFather() , false );
  auto & children = parentNode->get_children();
  children.erase( std::remove( children.begin() , children.end() ,
                               nodeToPrune ) , children.end() );
  auto & branches = parentNode->getBranches();
  branches.erase( std::remove( branches.begin() , branches.end() ,
                               nodeToPrune->get_f_change() ) ,
                  branches.end() );
  delete nodeToPrune;
  if( parentNode->get_children().empty() && parentNode->get_parent() )
   return( prune( parentNode , solvers ) );
  return( parentNode );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 std::vector< ExploringNode * > f_children;  ///< the children nodes

 int level;                                ///< level of the node in the tree

 ExploringNode * parent;                   ///< pointer to the parent node

/*--------------------------------------------------------------------------*/

 };  // end( class( ExploringNode ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS OpenList ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the set of open nodes of an exploration, abstracted over its discipline
/** The frontier of open nodes that a tree exploration keeps and picks the
 * next node from. Abstracting the discipline (which node comes next) behind
 * this interface lets a single exploration loop [see
 * BranchAndXSolver::explore()] serve every strategy: a LIFO stack gives
 * depth-first, a FIFO queue gives breadth-first, a priority queue ordered by
 * relaxation bound gives best-first. A new strategy is just a new OpenList. */

class OpenList {

 public:

 virtual ~OpenList() = default;

 /// whether there are no open nodes left
 [[nodiscard]] virtual bool empty( void ) const = 0;

 /// add a node to the open set
 virtual void push( ExploringNode * node ) = 0;

 /// remove and return the next node to explore
 virtual ExploringNode * pop( void ) = 0;

 /// whether the discipline is last-in first-out (a stack)
 /** Tells the exploration how to order a node's freshly generated children
  * in the open set so that the most promising one [the first returned by
  * RelaxationSolver::branch()] is explored first: a LIFO stack must receive
  * them in reverse, the other disciplines in branching order. */

 [[nodiscard]] virtual bool isLIFO( void ) const { return( false ); }

 };  // end( class( OpenList ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* BranchAndXSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File BranchAndXSolver.h --------------------------*/
/*--------------------------------------------------------------------------*/

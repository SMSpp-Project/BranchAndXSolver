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
#include <memory>
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
 * A node whose evaluation is inconclusive, i.e., whose relaxation returns
 * kStopTime, kStopIter or kLowPrecision (say, a Lagrangian Dual not solved
 * to the required accuracy whose solution gives nothing to branch upon), is
 * set aside: it is neither branched upon nor fenced, the exploration goes on
 * with the other nodes, and its bound enters that of the whole solve as the
 * bound of an open node does. If, at the end, some node so set aside can
 * still improve the incumbent, compute() returns kLowPrecision rather than
 * kOK; a relaxation returning kError still stops the whole solve.
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
 * \par Bounds of the children.
 * A child starts from the bound of its parent, which is valid for it; if
 * the RelaxationSolver that branched derives from BranchBounder [see
 * BranchingRule.h] and gives a better valid bound of a child, which it may
 * know without evaluating it, the child starts from that one instead. In
 * the serial exploration, the children of a node are then taken best bound
 * first, and, unless the tree is retained [see intReoptimize], a child
 * whose bound cannot improve the incumbent is pruned without ever being
 * created, and therefore evaluated, whichever the bounding protocol.
 *
 * \par Parallel exploration.
 * Two orthogonal kinds of parallelism are available. Across the tree [see
 * intMaxThread and ParallelDFSSolve()], after an ordered ramp-up \f$ K \f$
 * workers, each driving its own copy of the inner Solver attached to its own
 * copy of the Block (so that applying the Changes needs no synchronization,
 * which requires the Block to give an identical copy of itself, see
 * Block::get_R3_Block()), repeatedly claim the oldest open subtree and
 * explore it depth-first, sharing only the incumbent and the pool of open
 * subtrees, and a worker donates a branch of its subtree when the pool runs
 * low; claiming in the order of discovery is what keeps the
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
 * \par Cutoff probing.
 * With intCutOffProbing the solve is a sequence of decision problems, as in
 * Posta, Ferland and Michelon (2012): rather than waiting for the inner
 * Solver to find an incumbent good enough for the pruning, and for the
 * tightenings that read it (say, fixing by reduced cost), an artificial
 * one, the target, is set close to the bound of the root, and each
 * exploration either finds a solution better than the target or proves that
 * there is none, in which case the target is moved away from the bound and
 * the exploration is repeated. Let \f$ D \f$ be the bound of the root, given
 * by its relaxations before the first exploration, and
 * \f$ \varepsilon( v ) = \max \{ \varepsilon_a \,,\, \varepsilon_r
 * \max \{ | v | , 1 \} \} \f$ the tolerance of the pruning against
 * \f$ v \f$. The \f$ k \f$-th target (\f$ k = 0, 1, \ldots \f$) is
 *  \f[
 *    t_k = D - \Delta_0 \sum_{h = 0}^{k} \varphi^h \;,
 *  \f]
 * with \f$ \Delta_0 \f$ = dblCutOffProbingGap, multiplied by
 * \f$ \max \{ | D | , 1 \} \f$ if intCutOffProbing is RelativeGap, and
 * \f$ \varphi \f$ = dblCutOffProbingFactor \f$ \geq 1 \f$: with
 * \f$ \varphi = 1 \f$ the target moves by the same step every time, with
 * \f$ \varphi = 2 \f$ the step doubles. The \f$ k \f$-th exploration
 * starts with \f$ t_k \f$ as the incumbent, i.e., it is the one with the
 * cutoff dblLwCutOff = \f$ t_k \f$ [see set_par( double )], and since a
 * node is fenced when \f$ D( \nu ) \leq t_k + \varepsilon( t_k ) \f$ it
 * ends in one of two ways. Either it finds a solution better than
 * \f$ t_k \f$, or it proves that the target is refuted, i.e., that
 * \f$ z^* \leq u_k = t_k + \varepsilon( t_k ) \f$. In the second case
 * the next target is tried, and every node of the next exploration, the
 * root included, starts from the bound \f$ u_k \f$ [see f_probeBound], so
 * that a node is fenced also when
 * \f$ u_k \leq \bar{z} + \varepsilon( \bar{z} ) \f$, whatever its own
 * \f$ D( \nu ) \f$ (\f$ u_{-1} = D \f$).
 *
 * The solve ends, and what it gives is optimal within the tolerances, in
 * one of three ways:
 *
 * - The first exploration that finds a solution better than its target
 *   \f$ t_k \f$ goes on as an ordinary one against the incumbent
 *   \f$ \bar{z} \f$ thus found, with the bound of each node no larger than
 *   \f$ u_{k-1} \f$; when it ends every node is fenced, hence
 *   \f$ z^* \leq \bar{z} + \varepsilon( \bar{z} ) \f$, the same guarantee
 *   as without probing. In particular it stops as soon as
 *   \f$ \bar{z} \geq u_{k-1} - \varepsilon( \bar{z} ) \f$: with integer
 *   objective values, \f$ \Delta_0 = \varphi = 1 \f$ and tolerances
 *   below 1 this happens at the first solution found, which is the method of
 *   Posta et al., where the target decreases by 1 at each refutation.
 *
 * - A target is never set beyond the true incumbent, the best of the
 *   cutoff dblLwCutOff and of the solution that the relaxations may have
 *   found at the root: when \f$ t_k \f$ would be, the last exploration is
 *   the ordinary one against it, with the bound \f$ u_{k-1} \f$, and its
 *   outcome (kCutOff or kInfeasible if it finds nothing better) is that of
 *   the solve.
 *
 * - The same last exploration follows the intCutOffProbingRuns-th refuted
 *   target, which bounds the number of explorations also on a problem with
 *   no solution better than the cutoff, where every target is refuted.
 *
 * An exploration stopped by a budget (dblMaxTime bounds the whole solve)
 * ends the solve, with a bound no smaller than the target, those of the
 * nodes fenced against it being so. The inner Solver keep their state
 * across the explorations (say, the multipliers of a Lagrangian Dual at the
 * root, which make its next evaluation cheap), and so does the solution
 * found at the root, if any; the tree is not retained [see intReoptimize].
 * The target is the cutoff that each relaxation is given at every node and
 * the incumbent of GlobalInformation::str_Incumbent, hence any tightening
 * that reads them, as fixing by reduced cost, uses the target: this is
 * the point of the method, since with an incumbent close to the bound such
 * a tightening fixes much more, and every exploration is that much smaller.
 * Of course, a refuted target costs a whole exploration, which the next one
 * largely repeats: the probing pays when \f$ D - z^* \f$ is small with
 * respect to the gap between \f$ D \f$ and the incumbents that the inner
 * Solver find.
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
 * Computing 113:92-114, 2018. \n
 * M. Posta, J. A. Ferland, P. Michelon, "An exact method with variable
 * fixing for solving the generalized assignment problem", Computational
 * Optimization and Applications 52(3):629-644, 2012. */

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
  intTreeExplorer ,                 ///< who explores the tree
  intPBeFSSpawnDepth ,              ///< the spawn depth of PBeFS
  intPBeFSBudget ,                  ///< the spawn budget of PBeFS
  intPBeFSSplitShared ,             ///< the SelfSplit of PBeFS shares
  intPBeFSRaceNodes ,               ///< the nodes of a racer of PBeFS
  intCutOffProbing ,                ///< the scheme of the cutoff probing
  intCutOffProbingRuns ,            ///< max refuted targets of the probing
  intLastBXSPar                     ///< first allowed new int parameter
  };

 /// public enum for the double algorithmic parameters

 enum dbl_par_type_BXS {
  dblCutOffProbingGap = dblLastAlgPar ,  ///< first step of the probing
  dblCutOffProbingFactor ,               ///< growth of the step
  dblLastBXSPar                          ///< first allowed new dbl parameter
  };

 /// public enum for the vector-of-int algorithmic parameters

 enum vint_par_type_BXS {
  vintPBeFSRace = vintLastAlgPar ,  ///< the racers of PBeFS
  vintLastBXSPar                    ///< first allowed new vint parameter
  };

 /// public enum for the string algorithmic parameters

 enum str_par_type_BXS {
  /// name of the BlockSolverConfig file providing the inner Solver
  strNameOfBlockSolverConfigurationFile = strLastAlgPar ,
  strPBeFSScheme ,                 ///< the scheme of PBeFS
  strLastBXSPar                    ///< first allowed new string parameter
  };

 /// public enum for the vector-of-string algorithmic parameters

 enum vstr_par_type_BXS {
  /// the BlockSolverConfig files of the workers but the first
  vstrWorkerBlockSolverConfigurationFile = vstrLastAlgPar ,
  vstrLastBXSPar                   ///< first allowed new vstr parameter
  };

 /// public enum for when the inner Solver are run on a node

 enum BoundingProtocol {
  Eager = 0 ,                       ///< when the node is created
  Lazy = 1                          ///< when the node is extracted
  };

/*--------------------------------------------------------------------------*/
 /// public enum for who explores the tree [see intTreeExplorer]

 enum TreeExplorer {
  OwnExplorer = 0 ,                 ///< the explorations of this class
  PBeFS = 1                         ///< the branch-and-bound of PBeFS
  };

/*--------------------------------------------------------------------------*/
 /// public enum for the schemes of the cutoff probing [see intCutOffProbing]

 enum CutOffProbing {
  NoProbing = 0 ,                   ///< one exploration against the cutoff
  AbsoluteGap = 1 ,                 ///< first step dblCutOffProbingGap
  RelativeGap = 2                   ///< the same times max( | D | , 1 )
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
                      maxNodes( INT_MAX ) , treeExplorer( OwnExplorer ) ,
                      f_incumbentCell( nullptr ) ,
                      f_bestBoundCell( nullptr ) , f_lfaCell( nullptr ) ,
                      f_treeRoot( nullptr ) ,
                      nodeBudget( INT_MAX ) , timeBudget( Inf< double >() ) ,
                      relTol( 0 ) , absTol( 0 ) ,
                      configurationRS( nullptr ) ,
                      changes( Modification::eModAll ) ,
                      subOptimalNodes() , infeasibleNodes() ,
                      integerNodes() { }

/*--------------------------------------------------------------------------*/
 /// destructor

 ~BranchAndXSolver() override {
  dropWorkerSolvers();
  discardRetainedTree();
  delete bestSolution;
  delete configurationRS;
  for( auto cfg : v_workerRS )
   delete cfg;
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
  * - intTreeExplorer [OwnExplorer]: who explores the tree (see
  *   TreeExplorer): the explorations of this class, or the branch-and-bound
  *   of PBeFS [see PBeFSSolve()] by the scheme that strPBeFSScheme names,
  *   on max( intMaxThread , 1 ) workers; the strategy, the bounding
  *   protocol, the tolerances and the budgets are the same, and the tree
  *   is never retained [see intReoptimize].
  *
  * - intPBeFSSpawnDepth [0]: the spawn depth of the schemes of PBeFS that
  *   make tasks of the nodes above a depth (depthbounded) or at a depth
  *   (ordered); a binary tree has at most \f$ 2^d \f$ nodes of depth
  *   \f$ d \f$, so that it has to grow with the workers; <= 0 means the
  *   default of PBeFS.
  *
  * - intPBeFSBudget [0]: the backtracks of a worker of the scheme budget of
  *   PBeFS after which the children left on its stack become tasks; <= 0
  *   means the default of PBeFS.
  *
  * - intPBeFSSplitShared [0]: if not 0, the workers of the scheme selfsplit
  *   of PBeFS share their incumbent, which prunes more but makes the nodes
  *   evaluated depend on the timing; with 0 each worker prunes only by its
  *   own incumbent, and the nodes evaluated depend only on the workers.
  *
  * - intPBeFSRaceNodes [0]: the nodes that each racer of PBeFS evaluates at
  *   most [see vintPBeFSRace]; <= 0 means the default of PBeFS.
  *
  * - intCutOffProbing [NoProbing]: the scheme of the cutoff probing (see
  *   CutOffProbing), by which the solve is a sequence of explorations, each
  *   against a target incumbent moved away from the bound of the root after
  *   each one that proves it too good [see the class description].
  *
  * - intCutOffProbingRuns [100]: the targets that the cutoff probing may
  *   find refuted before the last exploration, which is against the true
  *   cutoff; <= 0 means no limit, with which a problem with no solution
  *   better than the cutoff is explored with ever worse targets forever.
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
   case( intTreeExplorer ):
    if( ( value < OwnExplorer ) || ( value > PBeFS ) )
     throw( std::invalid_argument( "BranchAndXSolver::set_par: "
            "invalid intTreeExplorer" ) );
    treeExplorer = static_cast< TreeExplorer >( value );
    break;
   case( intPBeFSSpawnDepth ):
    pbefsSpawnDepth = std::max( 0 , value );
    break;
   case( intPBeFSBudget ):
    pbefsBudget = std::max( 0 , value );
    break;
   case( intPBeFSSplitShared ):
    pbefsSplitShared = ( value != 0 );
    break;
   case( intPBeFSRaceNodes ):
    pbefsRaceNodes = std::max( 0 , value );
    break;
   case( intCutOffProbing ):
    if( ( value < NoProbing ) || ( value > RelativeGap ) )
     throw( std::invalid_argument( "BranchAndXSolver::set_par: "
            "invalid intCutOffProbing" ) );
    probing = static_cast< CutOffProbing >( value );
    break;
   case( intCutOffProbingRuns ):
    probingRuns = std::max( 0 , value );
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
  * - dblMaxTime [+INF]: the time budget of a solve, in seconds; each
  *   compute() of an inner Solver has the dblMaxTime of its own or the time
  *   left, whichever is smaller, and a node that it leaves unsolved is set
  *   aside with the bound it has [see v_unresolved];
  *
  * - dblRelAcc [1e-6] / dblAbsAcc [+INF, i.e., not active]: the relative and
  *   absolute tolerances of the pruning, a node being pruned when its
  *   relaxation bound cannot improve the incumbent by more than the largest
  *   of the two (the relative one taken w.r.t. max( | incumbent | , 1 )).
  *
  * - dblUpCutOff [+INF] when minimizing, dblLwCutOff [-INF] when
  *   maximizing: a value that no solution of interest is worse than, which
  *   the exploration starts from as the incumbent, so that a node is pruned
  *   against it until a better solution is found; if none is, compute()
  *   returns kCutOff with no solution, and get_lb() (get_ub() when
  *   maximizing) gives the cutoff. A change of it discards the tree retained
  *   for reoptimization, which may have been pruned against the old one.
  *
  * The double parameters specific of BranchAndXSolver are those of the
  * cutoff probing [see intCutOffProbing and the class description]:
  *
  * - dblCutOffProbingGap [1]: \f$ \Delta_0 > 0 \f$, the distance of the
  *   first target from the bound of the root, absolute or relative to it
  *   as intCutOffProbing says; 1 is the step of the integer objective
  *   values.
  *
  * - dblCutOffProbingFactor [1]: \f$ \varphi \geq 1 \f$, the factor by
  *   which the step of the target grows at each refuted target. */

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
   case( dblUpCutOff ):
    if( value != upCutOff )
     changes = Modification::eModAll;
    upCutOff = value;
    break;
   case( dblLwCutOff ):
    if( value != lwCutOff )
     changes = Modification::eModAll;
    lwCutOff = value;
    break;
   case( dblCutOffProbingGap ):
    if( ! ( value > 0 ) )
     throw( std::invalid_argument( "BranchAndXSolver::set_par: "
            "dblCutOffProbingGap must be positive" ) );
    probingGap = value;
    break;
   case( dblCutOffProbingFactor ):
    if( ! ( value >= 1 ) )
     throw( std::invalid_argument( "BranchAndXSolver::set_par: "
            "dblCutOffProbingFactor must be at least 1" ) );
    probingFactor = value;
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
  * - strPBeFSScheme ["shared"]: the name of the scheme of the
  *   branch-and-bound of PBeFS that explores the tree when intTreeExplorer
  *   is PBeFS, which is given to PBeFS as it is, so that every scheme of
  *   PBeFS can be named [see bb_scheme() in PBeFS]; an unknown name throws
  *   std::invalid_argument.
  *
  * The string parameters are taken by the rvalue version, which the one
  * taking a const reference forwards to [see ThinComputeInterface]. */

 void set_par( idx_type par , std::string && value ) override {
  if( par == strPBeFSScheme ) {
   checkPBeFSScheme( value );
   pbefsScheme = std::move( value );
   return;
   }
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
 /// set the vector-of-int parameters of BranchAndXSolver
 /** - vintPBeFSRace [empty]: the strategies [see intSolveMethod] of the
  *   racers that the branch-and-bound of PBeFS [see intTreeExplorer] starts
  *   with: the racer k explores the tree from the root by the k-th of them
  *   with the inner Solver of the worker k [see
  *   vstrWorkerBlockSolverConfigurationFile], all of them at once and
  *   pruning by the best incumbent of all, for at most intPBeFSRaceNodes
  *   nodes, and the scheme of strPBeFSScheme goes on from the open nodes of
  *   the racer with the best bound; the racers beyond the workers are left
  *   out, and the empty vector (the default) means no race. A strategy out
  *   of SolveMethod throws std::invalid_argument. */

 void set_par( idx_type par , std::vector< int > && value ) override {
  if( par != vintPBeFSRace ) {
   Solver::set_par( par , std::move( value ) );
   return;
   }
  for( const auto m : value )
   if( ( m < DFS ) || ( m > BestEstimate ) )
    throw( std::invalid_argument( "BranchAndXSolver::set_par: invalid "
                                  "strategy in vintPBeFSRace" ) );
  pbefsRace = std::move( value );
  }

/*--------------------------------------------------------------------------*/
 /// set the vector-of-string parameters of BranchAndXSolver
 /** - vstrWorkerBlockSolverConfigurationFile [empty]: the names of the
  *   files of the BlockSolverConfig of the workers of the parallel
  *   explorations but the first, the i-th entry for the worker i + 1, whose
  *   inner Solver are created out of it rather than out of the one of
  *   strNameOfBlockSolverConfigurationFile, which is that of the first
  *   worker, of the empty entries and of the workers beyond the vector;
  *   every BlockSolverConfig must describe the same kinds of inner Solver
  *   (:RelaxationSolver and heuristic :ChangeSolver), possibly differently
  *   configured, since each Change of the tree is applied to the inner
  *   Solver of whichever worker explores the node. A file out of which no
  *   BlockSolverConfig can be created throws std::invalid_argument. */

 void set_par( idx_type par , std::vector< std::string > && value )
  override {
  if( par != vstrWorkerBlockSolverConfigurationFile ) {
   Solver::set_par( par , std::move( value ) );
   return;
   }
  std::vector< BlockSolverConfig * > cfgs;
  for( const auto & name : value ) {
   BlockSolverConfig * cfg = nullptr;
   if( ! name.empty() ) {
    auto c = Configuration::deserialize( name );
    cfg = dynamic_cast< BlockSolverConfig * >( c );
    if( ! cfg ) {
     delete c;
     for( auto d : cfgs )
      delete d;
     throw( std::invalid_argument( "BranchAndXSolver::set_par: unable to "
            "create a BlockSolverConfig out of " + name ) );
     }
    }
   cfgs.push_back( cfg );
   }
  for( auto cfg : v_workerRS )
   delete cfg;
  v_workerRS = std::move( cfgs );
  v_nameWorkerRS = std::move( value );
  }

/*--------------------------------------------------------------------------*/
 /// the kinds of Modification the BranchAndXSolver reads
 /** The inner Solver are not registered to the Block: they have the
  * Modification from the BranchAndXSolver [see add_Modification()], which
  * therefore reads what any of them reads [see Solver::concerned_by()], and
  * all the kinds before they are created. What the tree needs to know of a
  * change [see keeps_tree()] is what its RelaxationSolver read of it: a
  * change that none of them reads does not change the bounds of the
  * nodes. */

 [[nodiscard]] Modification::ModConcern concerned_by( void ) const override {
  if( v_created.empty() )
   return( Modification::eModAnything );
  Modification::ModConcern c = 0;
  for( auto slvr : v_created )
   c |= slvr->concerned_by();
  return( c );
  }

/*--------------------------------------------------------------------------*/
 /// receive a Modification and forward it to the (private) inner Solver
 /** The Modification is passed to each of the inner Solver that reads its
  * kind [see Modification::is_of_concern_to()], as the Block would do if they
  * were registered to it. */

 void add_Modification( sp_Mod & mod ) override {
  Solver::add_Modification( mod );
  for( auto slvr : v_created )
   if( mod->is_of_concern_to( slvr->concerned_by() ) )
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

 /// whether an incumbent solution is actually available
 /** A solve stopped before finding an incumbent has no solution to write.
  * bestSolution alone can still belong to an earlier solve, whose bound
  * has been reset to no_incumbent() at the start of a new exploration. */

 [[nodiscard]] bool has_var_solution( void ) override {
  return( bestSolution && ( bestBound != no_incumbent() ) );
  }

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
 /// number of node evaluations started by the last compute(), root included
 /** Counts each node once per evaluation, independently of the number of
  * inner solvers; reoptimization counts the frontier evaluations of this
  * solve. Nodes fenced without evaluation are excluded. Zero before the
  * first compute(), or when an unchanged problem reuses its answer. */

 [[nodiscard]] long get_evaluated_nodes( void ) const {
  return( f_evaluated_nodes.load( std::memory_order_relaxed ) );
  }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the BranchAndXSolver
 *  @{ */

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastBXSPar ) );
  }

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override {
  return( idx_type( dblLastBXSPar ) );
  }

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return( idx_type( strLastBXSPar ) );
  }

 [[nodiscard]] idx_type get_num_vint_par( void ) const override {
  return( idx_type( vintLastBXSPar ) );
  }

 [[nodiscard]] idx_type get_num_vstr_par( void ) const override {
  return( idx_type( vstrLastBXSPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  switch( par ) {
   case( intSolveMethod ):               return( int( BestFS ) );
   case( intThreadForDifferentSolvers ): return( 1 );
   case( intBoundingProtocol ):          return( int( Eager ) );
   case( intReoptimize ):                return( 0 );
   case( intMaxNodes ):                  return( INT_MAX );
   case( intTreeExplorer ):              return( int( OwnExplorer ) );
   case( intPBeFSSpawnDepth ):           return( 0 );
   case( intPBeFSBudget ):               return( 0 );
   case( intPBeFSSplitShared ):          return( 0 );
   case( intPBeFSRaceNodes ):            return( 0 );
   case( intCutOffProbing ):             return( int( NoProbing ) );
   case( intCutOffProbingRuns ):         return( 100 );
   }
  return( Solver::get_dflt_int_par( par ) );
  }

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override {
  switch( par ) {
   case( dblCutOffProbingGap ):    return( 1 );
   case( dblCutOffProbingFactor ): return( 1 );
   }
  return( Solver::get_dflt_dbl_par( par ) );
  }

 [[nodiscard]] const std::vector< int > & get_dflt_vint_par( idx_type par )
  const override {
  static const std::vector< int > none;
  if( par == vintPBeFSRace )
   return( none );
  return( Solver::get_dflt_vint_par( par ) );
  }

 [[nodiscard]] const std::vector< std::string > & get_dflt_vstr_par(
                                         idx_type par ) const override {
  static const std::vector< std::string > none;
  if( par == vstrWorkerBlockSolverConfigurationFile )
   return( none );
  return( Solver::get_dflt_vstr_par( par ) );
  }

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override {
  if( par == strNameOfBlockSolverConfigurationFile ) {
   static const std::string dfltRS;
   return( dfltRS );
   }
  if( par == strPBeFSScheme ) {
   static const std::string dfltScheme = "shared";
   return( dfltScheme );
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
   case( intTreeExplorer ):              return( int( treeExplorer ) );
   case( intPBeFSSpawnDepth ):           return( pbefsSpawnDepth );
   case( intPBeFSBudget ):               return( pbefsBudget );
   case( intPBeFSSplitShared ):          return( pbefsSplitShared ? 1 : 0 );
   case( intPBeFSRaceNodes ):            return( pbefsRaceNodes );
   case( intCutOffProbing ):             return( int( probing ) );
   case( intCutOffProbingRuns ):         return( probingRuns );
   }
  return( Solver::get_int_par( par ) );
  }

 [[nodiscard]] const std::vector< int > & get_vint_par( idx_type par )
  const override {
  if( par == vintPBeFSRace )
   return( pbefsRace );
  return( Solver::get_vint_par( par ) );
  }

 [[nodiscard]] const std::vector< std::string > & get_vstr_par(
                                         idx_type par ) const override {
  if( par == vstrWorkerBlockSolverConfigurationFile )
   return( v_nameWorkerRS );
  return( Solver::get_vstr_par( par ) );
  }

 [[nodiscard]] double get_dbl_par( idx_type par ) const override {
  switch( par ) {
   case( dblMaxTime ): return( maxTime );
   case( dblRelAcc ):  return( relAcc );
   case( dblAbsAcc ):  return( absAcc );
   case( dblUpCutOff ): return( upCutOff );
   case( dblLwCutOff ): return( lwCutOff );
   case( dblCutOffProbingGap ):    return( probingGap );
   case( dblCutOffProbingFactor ): return( probingFactor );
   }
  return( Solver::get_dbl_par( par ) );
  }

 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override {
  if( par == strNameOfBlockSolverConfigurationFile )
   return( nameRS );
  if( par == strPBeFSScheme )
   return( pbefsScheme );
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
  if( name == "intTreeExplorer" )
   return( intTreeExplorer );
  if( name == "intPBeFSSpawnDepth" )
   return( intPBeFSSpawnDepth );
  if( name == "intPBeFSBudget" )
   return( intPBeFSBudget );
  if( name == "intPBeFSSplitShared" )
   return( intPBeFSSplitShared );
  if( name == "intPBeFSRaceNodes" )
   return( intPBeFSRaceNodes );
  if( name == "intCutOffProbing" )
   return( intCutOffProbing );
  if( name == "intCutOffProbingRuns" )
   return( intCutOffProbingRuns );
  return( Solver::int_par_str2idx( name ) );
  }

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override {
  if( name == "dblCutOffProbingGap" )
   return( dblCutOffProbingGap );
  if( name == "dblCutOffProbingFactor" )
   return( dblCutOffProbingFactor );
  return( Solver::dbl_par_str2idx( name ) );
  }

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override {
  if( name == "strNameOfBlockSolverConfigurationFile" )
   return( strNameOfBlockSolverConfigurationFile );
  if( name == "strPBeFSScheme" )
   return( strPBeFSScheme );
  return( Solver::str_par_str2idx( name ) );
  }

 [[nodiscard]] idx_type vint_par_str2idx( const std::string & name )
  const override {
  if( name == "vintPBeFSRace" )
   return( vintPBeFSRace );
  return( Solver::vint_par_str2idx( name ) );
  }

 [[nodiscard]] idx_type vstr_par_str2idx( const std::string & name )
  const override {
  if( name == "vstrWorkerBlockSolverConfigurationFile" )
   return( vstrWorkerBlockSolverConfigurationFile );
  return( Solver::vstr_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::string pars[] = { "intSolveMethod" ,
                                      "intThreadForDifferentSolvers" ,
                                      "intBoundingProtocol" ,
                                      "intReoptimize" , "intMaxNodes" ,
                                      "intTreeExplorer" ,
                                      "intPBeFSSpawnDepth" ,
                                      "intPBeFSBudget" ,
                                      "intPBeFSSplitShared" ,
                                      "intPBeFSRaceNodes" ,
                                      "intCutOffProbing" ,
                                      "intCutOffProbingRuns" };
  if( ( idx >= intSolveMethod ) && ( idx < intLastBXSPar ) )
   return( pars[ idx - intSolveMethod ] );
  return( Solver::int_par_idx2str( idx ) );
  }

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override {
  static const std::string pars[] = { "dblCutOffProbingGap" ,
                                      "dblCutOffProbingFactor" };
  if( ( idx >= dblCutOffProbingGap ) && ( idx < dblLastBXSPar ) )
   return( pars[ idx - dblCutOffProbingGap ] );
  return( Solver::dbl_par_idx2str( idx ) );
  }

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override {
  if( idx == strNameOfBlockSolverConfigurationFile ) {
   static const std::string par = "strNameOfBlockSolverConfigurationFile";
   return( par );
   }
  if( idx == strPBeFSScheme ) {
   static const std::string par = "strPBeFSScheme";
   return( par );
   }
  return( Solver::str_par_idx2str( idx ) );
  }

 [[nodiscard]] const std::string & vint_par_idx2str( idx_type idx )
  const override {
  if( idx == vintPBeFSRace ) {
   static const std::string par = "vintPBeFSRace";
   return( par );
   }
  return( Solver::vint_par_idx2str( idx ) );
  }

 [[nodiscard]] const std::string & vstr_par_idx2str( idx_type idx )
  const override {
  if( idx == vstrWorkerBlockSolverConfigurationFile ) {
   static const std::string par = "vstrWorkerBlockSolverConfigurationFile";
   return( par );
   }
  return( Solver::vstr_par_idx2str( idx ) );
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
  * cooperation [see GlobalInformation.h], i.e., the incumbent, the best
  * bound and the local-fixing-allowed flag. The references to the three
  * individual atomics are cached once here, and used lock-free ever after;
  * the incumbent starts at "no incumbent" and the best bound at "no bound",
  * which is why this needs the Block: which infinity that is depends on the
  * sense of the Objective. Called by set_Block(), i.e., before any inner
  * Solver is handed the GlobalInformation. */

 void initializeGlobalInformation( void ) {
  if( ! f_incumbentCell ) {   // the first Block: declare the Collection
   f_globalInfo.add_to_Universe< std::atomic< double > >(
				GlobalInformation::str_AtomicScalars );
   auto scalars = f_globalInfo.get_from_Universe< std::atomic< double > >(
					   GlobalInformation::str_AtomicScalars );
   f_incumbentCell = &( ( * scalars )[ GlobalInformation::str_Incumbent ] );
   f_bestBoundCell = &( ( * scalars )[ GlobalInformation::str_BestBound ] );

   f_globalInfo.add_to_Universe< std::atomic< bool > >(
				GlobalInformation::str_AtomicFlags );
   f_lfaCell = &( ( * f_globalInfo.get_from_Universe< std::atomic< bool > >(
		     GlobalInformation::str_AtomicFlags ) )[
		      GlobalInformation::str_LocalFixingAllowed ] );
   }
  f_incumbentCell->store( no_incumbent() );
  f_bestBoundCell->store( - no_incumbent() );
  f_lfaCell->store( true );
  }

/*--------------------------------------------------------------------------*/
 /// drop what belongs to the current Block, before attaching to another
 /** Deletes the inner Solver, the tree retained for reoptimization and the
  * best Solution, all of which refer to the current Block, and the pool of
  * the globally valid cuts that the relaxations may have created in the
  * GlobalInformation [see str_GlobalCuts]. */

 void dropInnerSolvers( void ) {
  dropWorkerSolvers();
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
 /// the value that the exploration starts from as the incumbent
 /** dblUpCutOff for a minimization problem, dblLwCutOff for a maximization
  * one [see set_par( double )]: no_incumbent() when it is left at its
  * default. */

 [[nodiscard]] double start_incumbent( void ) const;

/*--------------------------------------------------------------------------*/
 /// create the per-worker private Solver sets for the parallel exploration
 /** Creates the \p K Solver sets of the workers of the parallel exploration
  * (see v_workerSolvers), so that each worker applies the branching Changes
  * to nothing but its own Solver and its own Block, and the incumbent and
  * the pool of open subtrees are the only shared state. The first worker
  * drives the inner Solver of f_Block; each of the others gets an identical
  * copy of f_Block [see copy_of_Block()] and a set of inner Solver created
  * out of configurationRS and registered to it, so that they have the
  * Modification of the copy as those of f_Block have the ones of f_Block
  * [see add_Modification()]. The Changes are applied to any of them, since
  * they identify what they change by its AbstractPath relative to the
  * Block. The copies and their Solver are owned, and they live until
  * dropWorkerSolvers(), as they do not receive the Modification of
  * f_Block. */

 void createWorkerSolvers( Index K ) {
  if( ! ( configurationRS && f_Block ) )
   throw( std::logic_error( "BranchAndXSolver::createWorkerSolvers: "
          "configurationRS or f_Block not set" ) );
  dropWorkerSolvers();
  v_workerSolvers.resize( K );
  if( ! K )
   return;
  v_workerSolvers[ 0 ].relaxation = f_RelaxationSolvers;
  v_workerSolvers[ 0 ].heuristic = f_HeuristicSolvers;
  for( Index w = 1 ; w < K ; ++w ) {
   auto & ws = v_workerSolvers[ w ];
   ws.block = copy_of_Block();
   // the BlockSolverConfig of the worker, if it has one of its own [see
   // vstrWorkerBlockSolverConfigurationFile]
   const auto rs = ( w - 1 < v_workerRS.size() ) && v_workerRS[ w - 1 ] ?
                   v_workerRS[ w - 1 ] : configurationRS;
   for( Index i = 0 ; i < rs->num_ComputeConfig() ; ++i ) {
    const auto name = rs->get_SolverName( i );
    std::unique_ptr< Solver > s( Solver::new_Solver( name ) );
    if( auto cfg = rs->get_SolverConfig( i ) )
     s->set_ComputeConfig( cfg );
    auto slvr = s.release();
    ws.block->register_Solver( slvr );
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
 /// delete the per-worker Solver sets, and the copies of the Block
 /** The Solver of the first worker are those of f_Block, and they stay. */

 void dropWorkerSolvers( void ) {
  for( auto & ws : v_workerSolvers )
   if( ws.block ) {
    ws.block->unregister_Solvers( true );
    delete ws.block;
    }
  v_workerSolvers.clear();
  }

/*--------------------------------------------------------------------------*/
 /// an identical copy of f_Block, for a worker of the parallel exploration
 /** The copy is the R3 Block that f_Block gives with no Configuration [see
  * Block::get_R3_Block()], with the weights of its Objective and of its
  * Variable; a Block that gives none cannot be explored in parallel, and
  * this throws. */

 [[nodiscard]] Block * copy_of_Block( void ) const {
  static const std::string _prfx = "BranchAndXSolver::copy_of_Block: ";
  static const std::string why = ", and every worker of the parallel "
   "exploration but the first needs one";
  Block * copy;
  try {
   copy = f_Block->get_R3_Block( nullptr );
   }
  catch( const std::exception & e ) {
   throw( std::invalid_argument( _prfx + "the Block cannot be copied ("
                                 + e.what() + ")" + why ) );
   }
  if( ! copy )
   throw( std::invalid_argument( _prfx + "the Block gives no copy of "
                                 "itself [see Block::get_R3_Block()]" +
                                 why ) );
  copy->set_objective_weight( f_Block->get_objective_weight() , eNoMod );
  copy->set_variable_weight( f_Block->get_variable_weight() , eNoMod );
  return( copy );
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
 /// the exploration of the tree by the branch-and-bound of PBeFS
 /** The tree explored by the branch-and-bound of PBeFS, by the scheme that
  * strPBeFSScheme names [see intTreeExplorer], with \p K workers, each
  * driving its own set of inner Solver [see createWorkerSolvers()]. A node
  * is the Change that generates it from its parent, linked to the node of
  * the parent, so that any worker can reach it: the worker moves its Solver
  * from its current node to the target one by undoing the Changes up to
  * their deepest common ancestor (with undo Changes of its own) and
  * applying those down from it. The evaluation of a node is that of the
  * other explorations [see evaluateNode()], against the incumbent that
  * PBeFS gives the worker before each evaluation (by cutoff()), which the
  * worker keeps on its own together with the Solution it finds; hence, the
  * outcome is that of the scheme, also when the scheme makes it the same
  * at every run. The Solution stays in the node, and that of the incumbent
  * node is the one of the solve. With the eager bounding protocol the
  * children are evaluated when they are generated, with the lazy one when
  * they are extracted; with the eager one, a child whose evaluation has
  * found a feasible solution is given to the search even when it cannot
  * improve on that solution, and the solution is lost only if the
  * incumbent is within the tolerances of it when the child is extracted.
  * The bound of the search [see publish_bound()] is the one that PBeFS
  * gives each time it improves, and that of the root until it gives one;
  * a node set aside unresolved, which PBeFS takes as infeasible, bounds it
  * through set_aside(). */

 int PBeFSSolve( std::mutex & globalMutex , int K );

/*--------------------------------------------------------------------------*/
 /// throws std::invalid_argument if PBeFS has no scheme of the given name
 /// [see strPBeFSScheme]

 static void checkPBeFSScheme( const std::string & name );

/*--------------------------------------------------------------------------*/
 /// delete the tree retained for reoptimization (if any) and its frontier

 void discardRetainedTree( void );

/*--------------------------------------------------------------------------*/
 /// whether the time budget of the current solve is over

 [[nodiscard]] bool out_of_time( void ) const {
  return( std::chrono::steady_clock::now() >= f_deadline );
  }

/*--------------------------------------------------------------------------*/
 /// set aside a node whose evaluation is inconclusive, of the given bound
 /** Adds \p bound to v_unresolved, under the lock of the caller if any, and
  * improves f_unresolvedBound with it [see publish_bound()]. */

 void set_aside( double bound , bool minimizing );

/*--------------------------------------------------------------------------*/
 /// improve the best bound of the GlobalInformation
 /** Takes \p bound, the best among the relaxation bounds of the open nodes
  * and of the node being evaluated as the caller knows them, makes it no
  * better than those of the nodes set aside unresolved [see
  * f_unresolvedBound] and than the incumbent, and stores it in
  * f_bestBoundCell if it improves it (i.e., if it is larger for a
  * minimization problem, smaller for a maximization one), so that the cell
  * is monotone and a call that does not change it costs only the loads.
  * The incumbent may be a cutoff, or the target of the cutoff probing
  * [see intCutOffProbing], with which the value is still a bound, since a
  * solution better than it can only lie under an open node. */

 void publish_bound( double bound , bool minimizing );

/*--------------------------------------------------------------------------*/
 /// account for the nodes set aside unresolved in the global bound
 /** Improves f_open_bound with the relaxation bounds of the nodes set aside
  * unresolved in this solve [see v_unresolved], and returns true if some of
  * them can still improve the incumbent, i.e., if the solve, all the rest
  * being explored, is still not conclusive. */

 bool add_unresolved_bounds( bool minimizing );

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
 /// one exploration of the tree, by the explorer and the strategy asked
 /** The branch-and-bound of PBeFS [see PBeFSSolve()] if intTreeExplorer
  * says so, the parallel depth-first exploration [see ParallelDFSSolve()]
  * if the strategy is DFS on more than one thread, and the serial one [see
  * treeSolve()] otherwise; it starts from the incumbent in f_incumbentCell
  * and from the bound f_probeBound. */

 int runExploration( std::mutex & globalMutex );

/*--------------------------------------------------------------------------*/
 /// the sequence of explorations of the cutoff probing
 /** Evaluates the relaxations of the root, out of which the first target,
  * then runs an exploration [see runExploration()] per target, until one
  * finds a solution better than its target or is the last one, as the
  * class description details; returns the sol_type of that exploration,
  * and that of the only one if the root is infeasible, fenced or
  * inconclusive, or has no finite bound. */

 int probeCutOff( std::mutex & globalMutex );

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

 TreeExplorer treeExplorer;   ///< who explores the tree (see
                              ///< intTreeExplorer)

 /// the inherited dblMaxTime [see set_par( double )]
 double maxTime = Inf< double >();

 /// the inherited dblRelAcc [see set_par( double )]
 double relAcc = 1e-6;

 /// the inherited dblAbsAcc [see set_par( double )]
 double absAcc = Inf< double >();

 /// the inherited dblUpCutOff [see set_par( double )]
 double upCutOff = Inf< double >();

 /// the inherited dblLwCutOff [see set_par( double )]
 double lwCutOff = - Inf< double >();

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

 /// the cached reference to the atomic best-bound cell
 /** The global relaxation bound of the exploration under way, i.e., the
  * best among the relaxation bounds of the open nodes, of the node being
  * evaluated, of the nodes set aside unresolved, and the incumbent [see
  * get_lb()], which is how the relaxations see it lock-free [see
  * GlobalInformation::str_BestBound and publish_bound()]: it starts at "no
  * bound" (- no_incumbent()) at each compute(), and is only ever improved
  * during it. nullptr until set_Block(). */
 std::atomic< double > * f_bestBoundCell;

 /// the best of the bounds of the nodes set aside unresolved in this solve
 /** The best (smallest for a minimization problem) of v_unresolved, kept
  * lock-free by set_aside() so that publish_bound() can account for them
  * at no cost; no_incumbent() when there are none. */
 std::atomic< double > f_unresolvedBound{ 0 };

 /// the cached reference to the atomic local-fixing-allowed flag
 /** Set by each compute(), which is where it is known whether the tree is
  * being retained across re-solves [see intReoptimize]. nullptr until
  * set_Block(). */
 std::atomic< bool > * f_lfaCell;

 /// the root of the tree retained for reoptimization, nullptr if none
 ExploringNode * f_treeRoot;

 /// node evaluations in the current compute(), also shared by DFS workers
 std::atomic< long > f_evaluated_nodes{ 0 };

 /// children pruned before their evaluation in the current compute(), by
 /// the bounds that the RelaxationSolver gave of them [see BranchBounder]
 long f_bounded_children = 0;

 /// residual node budget of the current solve (from intMaxNodes), consumed
 /// by the exploration
 int nodeBudget;

 /// residual time budget of the current solve (from the inherited
 /// dblMaxTime), consumed by the exploration
 double timeBudget;

 /// the time by which the current solve has to end, max() if none
 std::chrono::steady_clock::time_point f_deadline =
  std::chrono::steady_clock::time_point::max();

 double relTol;   ///< the inherited dblRelAcc, cached at compute() start

 double absTol;   ///< the inherited dblAbsAcc, cached at compute() start

 std::string nameRS;          ///< name of the BlockSolverConfig file

 std::string pbefsScheme = "shared";  ///< the scheme of PBeFS (see
                                      ///< strPBeFSScheme)

 int pbefsSpawnDepth = 0;  ///< see intPBeFSSpawnDepth

 int pbefsBudget = 0;      ///< see intPBeFSBudget

 bool pbefsSplitShared = false;  ///< see intPBeFSSplitShared

 int pbefsRaceNodes = 0;   ///< see intPBeFSRaceNodes

 std::vector< int > pbefsRace;  ///< see vintPBeFSRace

 CutOffProbing probing = NoProbing;  ///< see intCutOffProbing

 int probingRuns = 100;              ///< see intCutOffProbingRuns

 double probingGap = 1;              ///< see dblCutOffProbingGap

 double probingFactor = 1;           ///< see dblCutOffProbingFactor

 /// the bound that every root of the current exploration starts from
 /** In the sense of the Block, the infinity that any bound improves upon
  * but in the explorations of the cutoff probing, where it is
  * \f$ u_{k-1} \f$, the bound of the root or that proven by the last
  * refuted target [see the class description]: the children inherit it,
  * so that it is the bound of a node whose relaxation gives a worse one. */
 double f_probeBound = 0;

 BlockSolverConfig * configurationRS;  ///< the inner BlockSolverConfig

 /// the names of vstrWorkerBlockSolverConfigurationFile
 std::vector< std::string > v_nameWorkerRS;

 /// the BlockSolverConfig of vstrWorkerBlockSolverConfigurationFile, owned,
 /// nullptr for an empty name
 std::vector< BlockSolverConfig * > v_workerRS;

 /// the (private) inner Solver created out of configurationRS, owned
 std::vector< Solver * > v_created;

 /// the private Solver set of one worker of the parallel exploration
 /** The Block the Solver are registered to, a copy of f_Block, is owned
  * together with them; none for the first worker, which drives the inner
  * Solver of f_Block [see createWorkerSolvers()]. */
 struct WorkerSolvers {
  std::vector< std::pair< RelaxationSolver * , Solver * > > relaxation;
  std::vector< std::pair< ChangeSolver * , Solver * > > heuristic;
  Block * block = nullptr;           ///< the copy of f_Block, owned
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

 /// the relaxation bounds of the nodes set aside unresolved in this solve
 /** A node whose evaluation is inconclusive (its relaxation returns
  * kStopTime, kStopIter or kLowPrecision) can be neither branched upon nor
  * fenced: it is set aside, and its relaxation bound (that of its father,
  * if its own evaluation has not improved it) is kept here, so that the
  * bound of the whole solve accounts for the part of the region that it
  * covers [see add_unresolved_bounds()]. Written under the incumbent mutex
  * by the workers of the parallel exploration. */
 std::vector< double > v_unresolved;

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
  * a child (whose entries are nullptr), are deleted; the estimates and the
  * bounds of the children are taken too, if the RelaxationSolver gives
  * them [see BranchEstimator and BranchBounder], and kept only if there is
  * one per Change. */
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
  bounds.clear();
  if( auto bb = dynamic_cast< BranchBounder * >( solver ) ) {
   bounds = bb->branch_bounds();
   if( bounds.size() != branches.size() )
    bounds.clear();
   }
  }

 /// return the Changes generating the children of this node
 std::vector< Change * > & getBranches( void ) { return( branches ); }

 /// forget the estimates and the bounds of the children, once they have
 /// been handed to them or when the branching Changes no longer match them
 void clear_estimates( void ) {
  estimates.clear();
  bounds.clear();
  }

 /// whether the RelaxationSolver that branched gave the bounds of the
 /// children [see obtainBranchList()]
 bool has_child_bounds( void ) const { return( ! bounds.empty() ); }

 /// the bound of the child generated by the i-th branching Change
 /** That of this node, valid for the child, whose region is a part of that
  * of the node, or the better one that the RelaxationSolver that branched
  * gave of the child [see BranchBounder], in the sense \p minimizing. */
 double child_bound( std::size_t i , bool minimizing ) const {
  if( ( i >= bounds.size() ) || std::isnan( bounds[ i ] ) )
   return( bound );
  return( minimizing ? std::max( bound , bounds[ i ] )
                     : std::min( bound , bounds[ i ] ) );
  }

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

 /// the bounds of the children, one per Change of branches, or none
 std::vector< double > bounds;

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
 /// whether prune_away() can prune the node, the :ChangeSolver being on
 /// the node \p at
 /** The pruning of an open node deletes it and the ancestors it leaves
  * childless, the root excepted; it can be done without the :ChangeSolver
  * if none of them is the node they are on. */

 static bool can_prune_away( const ExploringNode * node ,
                             const ExploringNode * at ) {
  if( ( node == at ) || ( ! node->parent ) )
   return( false );
  for( auto p = node->parent ; p->parent && ( p->f_children.size() == 1 ) ;
       p = p->parent )
   if( p == at )
    return( false );
  return( true );
  }

/*--------------------------------------------------------------------------*/
 /// prune an open node that the :ChangeSolver are not on, not moving them
 /** Deletes \p node and, recursively, any ancestor left childless, the root
  * excepted, as prune() does, but applying no undo Change, the :ChangeSolver
  * being elsewhere [see can_prune_away()]. */

 static void prune_away( ExploringNode * node ) {
  for( ; ; ) {
   ExploringNode * parentNode = node->parent;
   auto & children = parentNode->get_children();
   children.erase( std::remove( children.begin() , children.end() , node ) ,
                   children.end() );
   auto & branches = parentNode->getBranches();
   branches.erase( std::remove( branches.begin() , branches.end() ,
                                node->get_f_change() ) , branches.end() );
   delete node;
   if( ( ! parentNode->get_children().empty() ) || ( ! parentNode->parent ) )
    return;
   node = parentNode;
   }
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

 /// the relaxation bound of the next node, if the nodes come by their bound
 /** Writes in \p bound that of the node that pop() gives next and returns
  * true if no open node has a better relaxation bound, i.e., if the list
  * gives the nodes best bound first; returns false otherwise, and for an
  * empty list. */

 [[nodiscard]] virtual bool next_bound( double & bound ) const {
  return( false );
  }

 /// the best relaxation bound of the open nodes, if the list knows it
 /** Writes in \p bound the best (the smallest for a minimization problem)
  * relaxation bound among those of the open nodes and returns true, which
  * the exploration publishes as the bound of the search [see
  * GlobalInformation::str_BestBound]; returns false for an empty list, and
  * for one that does not keep it. A node does not change its bound while it
  * is in the list. */

 [[nodiscard]] virtual bool best_bound( double & bound ) const {
  return( false );
  }

 };  // end( class( OpenList ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* BranchAndXSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File BranchAndXSolver.h --------------------------*/
/*--------------------------------------------------------------------------*/

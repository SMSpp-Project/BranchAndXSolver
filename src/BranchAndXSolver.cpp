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

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BranchAndXSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0(BranchAndXSolver);

/*--------------------------------------------------------------------------*/
/*-------------------------- INTERNAL FUNCTIONS ----------------------------*/
/*--------------------------------------------------------------------------*/

/// tells if a dual bound cannot improve the incumbent beyond the tolerances
/** True if \p dual cannot improve \p best by more than
 * max( absAcc , relAcc * max( | best | , 1 ) ), i.e., the node can be
 * pruned within the required optimality tolerances. */

static bool cannot_improve(double dual, double best, bool minimizing,
                           double relAcc, double absAcc)
{
    if (std::isinf(best)) // no incumbent yet: everything can improve
        return (false);
    // an absAcc at its default +Inf means "not active" [see Solver::dblAbsAcc]
    const double eps = std::max(absAcc == Inf<double>() ? 0.0 : absAcc,
                                relAcc * std::max(std::abs(best), 1.0));
    return (minimizing ? dual >= best - eps : dual <= best + eps);
}

/// the Solution corresponding to the current solution of \p slvr
/** Writes the current solution of \p slvr into its Block (under lock) and
 * returns the corresponding newly minted Solution object, whose ownership
 * is transferred to the caller. This is how a Solution is obtained from a
 * heuristic :ChangeSolver, which unlike a :RelaxationSolver [see
 * RelaxationSolver::get_true_solution()] has no way of its own to produce
 * one. */

static Solution *solution_of(Solver *slvr,
                             Configuration *solc = nullptr)
{
    auto blck = slvr->get_Block();
    // TODO: check that lock()-ing / unlock()-ing here is appropriate
    blck->lock(slvr);
    slvr->get_var_solution(solc);
    auto sol = blck->get_Solution(solc);
    blck->unlock(slvr);
    return (sol);
}

/*--------------------------------------------------------------------------*/
/// write on the log one line per explored node, in CSV form
/** When a log stream is set [see Solver::set_log()], each node is described
 * by one line as soon as its fate is decided (right after its evaluation,
 * whenever that happens [see BranchAndXSolver::intBoundingProtocol]), so that
 * the whole exploration can be reconstructed offline: which nodes were
 * generated, in which order they were evaluated, what their dual bound was,
 * how the incumbent moved, and why each node was fenced. \p evalTime is the
 * time spent evaluating this very node, \p elapsed the time since the
 * beginning of the solve. */

static void logNode(std::ostream *log, int iter, ExploringNode *node,
                    double bestBound, double elapsed, double evalTime,
                    bool boundPruned)
{
    if (!log)
        return;
    auto father = node->get_parent();
    (*log) << iter << ',' << node->get_name() << ','
           << (father ? father->get_name() : -1) << ','
           << node->get_level() << ',' << bestBound << ','
           << node->get_dual_bound() << ',' << elapsed << ',' << evalTime
           << ',' << boundPruned << ',' << node->is_infeasible() << std::endl;
}

/*--------------------------------------------------------------------------*/
/// the header of the per-node CSV log [see logNode()]

static void logHeader(std::ostream *log)
{
    if (log)
        (*log) << "iter,node,father,level,incumbent,dualBound,"
                  "elapsedTime,evaluationTime,boundPruned,infeasible"
               << std::endl;
}

/*--------------------------------------------------------------------------*/
/// initialize the common per-solve variables
/** @param solvers list to be filled with all the :ChangeSolver to drive
 *  @param f_RelaxationSolvers the relaxation solvers to insert in solvers
 *  @param f_HeuristicSolvers the heuristic solvers to insert in solvers
 *  @param minimizing set to true if the problem is a minimization one */

static void initializeVariables(
    std::list<ChangeSolver *> *solvers,
    std::vector<std::pair<RelaxationSolver *,
                          Solver *>> *f_RelaxationSolvers,
    std::vector<std::pair<ChangeSolver *,
                          Solver *>> *f_HeuristicSolvers,
    bool &minimizing)
{
    if (!f_HeuristicSolvers->empty())
        minimizing = f_HeuristicSolvers->front().second->get_Block()->get_objective_sense() == Objective::eMin;
    else if (!f_RelaxationSolvers->empty())
        minimizing = f_RelaxationSolvers->front().second->get_Block()->get_objective_sense() == Objective::eMin;
    else
        throw(std::logic_error("BranchAndXSolver::initializeVariables: "
                               "both the HeuristicSolvers and the RelaxationSolvers are empty"));

    // the solver set moved along the tree is made of the traits (the moves
    // only ever apply() Changes)
    for (auto &p : *f_RelaxationSolvers)
        solvers->push_back(p.first);
    for (auto &p : *f_HeuristicSolvers)
        solvers->push_back(p.first);
}

/*--------------------------------------------------------------------------*/
/// move the given :ChangeSolver down to the son node currentNode
/** The first Solver produces (and stores in the node) the undo Change, the
 * others just apply the node Change. */

static void moveSolverToSon(Node *currentNode,
                            std::list<ChangeSolver *> *solvers)
{
    if (!currentNode->get_f_change())
        return;
    bool modified = false;
    for (const auto s : *solvers)
    {
        if (!modified)
            currentNode->set_toFather(s->apply(currentNode->get_f_change(),
                                               true));
        else
            s->apply(currentNode->get_f_change(), false);
        modified = true;
    }
}

/*--------------------------------------------------------------------------*/
/// move the given :ChangeSolver back to the father of currentNode

static void moveSolverToFather(Node *currentNode,
                               std::list<ChangeSolver *> *solvers)
{
    if (currentNode->get_toFather())
        for (const auto s : *solvers)
            s->apply(currentNode->get_toFather(), false);
}

/*--------------------------------------------------------------------------*/
// runs every Solver's compute() concurrently [defined below]; used by the
// parallel evaluation path of computeRelaxations() / computeHeuristic()
template <typename SolverPtr>
static bool computeAllParallel(const std::vector<SolverPtr> &slvrs,
                               std::vector<int> &zs, int maxThreads);

/*--------------------------------------------------------------------------*/
/// compute the relaxations at the current node
/** Computes every RelaxationSolver at the current node, updating the node
 * dual bound, the branching solver, and, when a true solution improves it,
 * the incumbent bestBound / bestSol; sets \p toPrune when the node can
 * be discarded. With \p nThreads > 1 (and no shared incumbent lock) the
 * several relaxations are computed concurrently, bit-identically to the serial
 * reduction [see computeAllParallel()].
 *  @return the sol_type [see Solver.h] of the computation */

static int computeRelaxations(
    std::vector<std::pair<RelaxationSolver *,
                          Solver *>> *f_RelaxationSolvers,
    Node *currentNode, const bool minimizing,
    double &bestBound, Solution *&bestSol,
    std::atomic<double> *incumbentCell,
    bool &toPrune, RelaxationSolver *&branchSolver,
    double relAcc = 0, double absAcc = 0,
    int nThreads = 1, std::mutex *incumbentMutex = nullptr)
{
    // parallel evaluation of the (several) relaxations of this node: only when more
    // than one thread is asked, there are at least two solvers, and no shared
    // incumbent lock is in force (a parallel-tree worker evaluates its per-node
    // solvers serially); the reduction is bit-identical to the serial path below
    if ((nThreads > 1) && (f_RelaxationSolvers->size() > 1) &&
        (!incumbentMutex))
    {
        const std::size_t n = f_RelaxationSolvers->size();
        std::vector<int> zs(n, INT_MIN); // INT_MIN: not computed (early stop)
        if (computeAllParallel(*f_RelaxationSolvers, zs, nThreads))
        {
            if (!currentNode->get_toFather()) // the root is infeasible
                return (Solver::kInfeasible);
            toPrune = true;
            currentNode->set_infeasible(true);
            return (Solver::kOK);
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            auto &[rs, s] = (*f_RelaxationSolvers)[i];
            if (zs[i] != ThinComputeInterface::kOK)
                return (zs[i]);
            if (rs->has_true_var_solution())
            {
                double primal_bound = minimizing ? rs->get_true_ub() : rs->get_true_lb();
                if (minimizing ? primal_bound < bestBound : primal_bound > bestBound)
                {
                    bestBound = primal_bound;
                    incumbentCell->store(bestBound);
                    delete bestSol;
                    bestSol = rs->get_true_solution();
                }
            }
            auto dualBound = minimizing ? s->get_lb() : s->get_ub();
            if (cannot_improve(dualBound, bestBound, minimizing, relAcc, absAcc))
            {
                toPrune = true;
                return (Solver::kOK);
            }
            if (minimizing ? dualBound > currentNode->get_dual_bound()
                           : dualBound < currentNode->get_dual_bound())
            {
                currentNode->set_dual_bound(dualBound);
                branchSolver = rs;
            }
        }
        return (Solver::kOK);
    }

    // serial evaluation
    for (auto &[rs, s] : *f_RelaxationSolvers)
    {
        auto z = s->compute();
        if (z == Solver::kInfeasible)
        {
            if (!currentNode->get_toFather()) // the root is infeasible
                return (Solver::kInfeasible);
            toPrune = true;
            currentNode->set_infeasible(true);
            break;
        }
        if (z != ThinComputeInterface::kOK)
            return (z);

        // see if the primal bound improves thanks to a true solution
        if (rs->has_true_var_solution())
        {
            double primal_bound = minimizing ? rs->get_true_ub() : rs->get_true_lb();
            if (minimizing ? primal_bound < bestBound : primal_bound > bestBound)
            {
                // in the parallel exploration the incumbent is shared between the
                // workers: re-check the improvement under the mutex
                std::unique_lock<std::mutex> guard;
                if (incumbentMutex)
                    guard = std::unique_lock<std::mutex>(*incumbentMutex);
                if (minimizing ? primal_bound < bestBound : primal_bound > bestBound)
                {
                    bestBound = primal_bound;
                    incumbentCell->store(bestBound);
                    delete bestSol;
                    bestSol = rs->get_true_solution();
                }
            }
        }

        // update the dual bound of the node and select the branching solver
        auto dualBound = minimizing ? s->get_lb() : s->get_ub();
        if (cannot_improve(dualBound, bestBound, minimizing, relAcc,
                           absAcc))
        {
            toPrune = true;
            return (Solver::kOK);
        }
        if (minimizing ? dualBound > currentNode->get_dual_bound()
                       : dualBound < currentNode->get_dual_bound())
        {
            currentNode->set_dual_bound(dualBound);
            branchSolver = rs;
        }
    }
    return (Solver::kOK);
}

/*--------------------------------------------------------------------------*/
/// compute the heuristics at the current node
/** Computes every heuristic ChangeSolver at the current node, updating the
 * incumbent bestBound / bestSol when improved; sets \p toPrune when the
 * node can be discarded.
 *  @return the sol_type [see Solver.h] of the computation */

static int computeHeuristic(
    std::vector<std::pair<ChangeSolver *,
                          Solver *>> *f_HeuristicSolvers,
    Node *currentNode, const bool minimizing,
    double &bestBound, Solution *&bestSol,
    std::atomic<double> *incumbentCell, bool &toPrune,
    int nThreads = 1, std::mutex *incumbentMutex = nullptr)
{
    // parallel evaluation of the (several) heuristics of this node, bit-identical
    // to the serial path below [see computeRelaxations()]
    if ((nThreads > 1) && (f_HeuristicSolvers->size() > 1) &&
        (!incumbentMutex))
    {
        const std::size_t n = f_HeuristicSolvers->size();
        std::vector<int> zs(n, INT_MIN); // INT_MIN: not computed (early stop)
        if (computeAllParallel(*f_HeuristicSolvers, zs, nThreads))
        {
            if (!currentNode->get_toFather()) // the root is infeasible
                return (Solver::kInfeasible);
            toPrune = true;
            return (Solver::kOK);
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            auto &[hs, s] = (*f_HeuristicSolvers)[i];
            if (zs[i] != ThinComputeInterface::kOK)
                return (zs[i]);
            double primal_bound = minimizing ? s->get_ub() : s->get_lb();
            if (s->has_var_solution() && s->is_var_feasible() &&
                (minimizing ? primal_bound < bestBound
                            : primal_bound > bestBound))
            {
                bestBound = primal_bound;
                incumbentCell->store(bestBound);
                delete bestSol;
                bestSol = solution_of(s);
            }
        }
        return (Solver::kOK);
    }

    // serial evaluation
    for (auto &[hs, s] : *f_HeuristicSolvers)
    {
        auto z = s->compute();
        if (z == Solver::kInfeasible)
        {
            if (!currentNode->get_toFather()) // the root is infeasible
                return (Solver::kInfeasible);
            toPrune = true;
            break;
        }
        if (z != ThinComputeInterface::kOK)
            return (z);

        // see if the primal bound improves
        double primal_bound = minimizing ? s->get_ub() : s->get_lb();
        if (s->has_var_solution() && s->is_var_feasible() &&
            (minimizing ? primal_bound < bestBound
                        : primal_bound > bestBound))
        {
            std::unique_lock<std::mutex> guard;
            if (incumbentMutex)
                guard = std::unique_lock<std::mutex>(*incumbentMutex);
            if (minimizing ? primal_bound < bestBound : primal_bound > bestBound)
            {
                bestBound = primal_bound;
                incumbentCell->store(bestBound);
                delete bestSol;
                bestSol = solution_of(s);
            }
        }
    }
    return (Solver::kOK);
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

template <typename SolverPtr>
static bool computeAllParallel(const std::vector<SolverPtr> &slvrs,
                               std::vector<int> &zs, int maxThreads)
{
    const int n = int(slvrs.size());
    const int K = std::max(1, std::min(maxThreads, n));
    std::atomic<int> next(0);
    std::atomic<bool> infeasible(false);
    auto work = [&]()
    {
        for (int i; (!infeasible.load(std::memory_order_relaxed)) &&
                    ((i = next.fetch_add(1)) < n);)
        {
            zs[i] = slvrs[i].second->compute();
            if (zs[i] == Solver::kInfeasible)
                infeasible.store(true, std::memory_order_relaxed);
        }
    };
    std::vector<std::thread> pool;
    pool.reserve(K - 1);
    for (int t = 0; t < K - 1; ++t)
        pool.emplace_back(work);
    work(); // the calling thread is a worker too
    for (auto &th : pool)
        th.join();
    return (infeasible.load());
}

/*--------------------------------------------------------------------------*/
/// evaluate the root node and produce its branching list

static int initializeRoot(
    std::vector<std::pair<RelaxationSolver *,
                          Solver *>> *f_RelaxationSolvers,
    std::vector<std::pair<ChangeSolver *,
                          Solver *>> *f_HeuristicSolvers,
    Node *rootNode, const bool minimizing,
    double &bestBound, Solution *&bestSol,
    std::atomic<double> *incumbentCell,
    RelaxationSolver *&branchSolver,
    std::mutex &globalMutex)
{
    rootNode->initializeBound(minimizing);
    rootNode->set_evaluated(true);
    bool toPrune = false;
    int res = computeRelaxations(f_RelaxationSolvers, rootNode, minimizing,
                                 bestBound, bestSol, incumbentCell,
                                 toPrune, branchSolver);
    if (toPrune || (res != ThinComputeInterface::kOK))
        return (res);
    res = computeHeuristic(f_HeuristicSolvers, rootNode, minimizing,
                           bestBound, bestSol, incumbentCell, toPrune);
    if (toPrune || (res != ThinComputeInterface::kOK))
        return (res);
    rootNode->obtainBranchList(branchSolver);
    return (Solver::kOK);
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
    Node *node, std::list<ChangeSolver *> &solvers,
    std::vector<std::pair<RelaxationSolver *,
                          Solver *>> *relaxation,
    std::vector<std::pair<ChangeSolver *,
                          Solver *>> *heuristic,
    bool minimizing, double &bestBound, Solution *&bestSol,
    std::atomic<double> *incumbentCell,
    bool &toPrune, RelaxationSolver *&branchSolver,
    int nThreads,
    double relAcc, double absAcc, std::mutex &globalMutex,
    std::mutex *incumbentMutex = nullptr,
    bool moveToNode = true)
{
    if (moveToNode)                      // with the lazy protocol the solvers have
        moveSolverToSon(node, &solvers); // already been moved to the node
    node->initializeBound(minimizing);
    node->set_evaluated(true);

    // computeRelaxations() / computeHeuristic() dispatch internally on nThreads:
    // a worker of the parallel tree exploration (incumbentMutex set) evaluates its
    // per-node solvers serially, otherwise nThreads > 1 runs them in parallel
    int res = computeRelaxations(relaxation, node, minimizing, bestBound,
                                 bestSol, incumbentCell, toPrune,
                                 branchSolver, relAcc,
                                 absAcc, nThreads, incumbentMutex);
    if ((res != Solver::kOK) || toPrune)
        return (res);

    return (computeHeuristic(heuristic, node, minimizing, bestBound, bestSol,
                             incumbentCell, toPrune, nThreads,
                             incumbentMutex));
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OPEN-LIST DISCIPLINES -------------------------*/
/*--------------------------------------------------------------------------*/
// the three open-list disciplines, as thin adapters that give the standard
// container adapters a common interface [see OpenList]: a LIFO stack (depth-
// first), a FIFO queue (breadth-first), a dual-bound priority queue (best-
// first). The storage is entirely std::stack / std::queue / std::priority_-
// queue, there is no hand-rolled data structure here

namespace
{

    /// LIFO open list: a std::stack, giving depth-first exploration

    class StackOpenList final : public OpenList
    {

    public:
        [[nodiscard]] bool empty(void) const override { return (c.empty()); }

        void push(ExploringNode *node) override { c.push(node); }

        ExploringNode *pop(void) override
        {
            auto node = c.top();
            c.pop();
            return (node);
        }

        [[nodiscard]] bool isLIFO(void) const override { return (true); }

    private:
        std::stack<ExploringNode *> c;

    }; // end( class( StackOpenList ) )

    /*--------------------------------------------------------------------------*/
    /// FIFO open list: a std::queue, giving breadth-first exploration

    class QueueOpenList final : public OpenList
    {

    public:
        [[nodiscard]] bool empty(void) const override { return (c.empty()); }

        void push(ExploringNode *node) override { c.push(node); }

        ExploringNode *pop(void) override
        {
            auto node = c.front();
            c.pop();
            return (node);
        }

    private:
        std::queue<ExploringNode *> c;

    }; // end( class( QueueOpenList ) )

    /*--------------------------------------------------------------------------*/
    /// dual-bound priority open list: a std::priority_queue, giving best-first

    class PriorityOpenList final : public OpenList
    {

        using Cmp = std::function<bool(ExploringNode *, ExploringNode *)>;

    public:
        explicit PriorityOpenList(Cmp cmp) : c(std::move(cmp)) {}

        [[nodiscard]] bool empty(void) const override { return (c.empty()); }

        void push(ExploringNode *node) override { c.push(node); }

        ExploringNode *pop(void) override
        {
            auto node = c.top();
            c.pop();
            return (node);
        }

    private:
        std::priority_queue<ExploringNode *, std::vector<ExploringNode *>, Cmp>
            c;

    }; // end( class( PriorityOpenList ) )

    /*--------------------------------------------------------------------------*/
    /// best-first open list with depth-first dives
    /** The global frontier is a dual-bound priority queue, but once a node is taken
     * the search dives straight into its most promising child down to a leaf,
     * leaving the siblings to the global frontier. The dive reaches a complete
     * (feasible) solution quickly, so the incumbent tightens early and the pruning
     * bites sooner. The children of the just-expanded node arrive through push()
     * before the next pop(): pop() routes the most promising of them onward
     * (continuing the dive) and spills the rest to the priority queue; when no child
     * arrives the dive has bottomed out and the next global best is taken. */

    class DiveOpenList final : public OpenList
    {

        using Cmp = std::function<bool(ExploringNode *, ExploringNode *)>;

    public:
        explicit DiveOpenList(Cmp cmp) : best(std::move(cmp)) {}

        [[nodiscard]] bool empty(void) const override
        {
            return (best.empty() && children.empty());
        }

        void push(ExploringNode *node) override { children.push_back(node); }

        ExploringNode *pop(void) override
        {
            if (!children.empty())
            {
                // the loop is isLIFO(), so it pushed the children in reverse branching
                // order: back() is the first (most promising) branch, the dive follows it
                auto node = children.back();
                children.pop_back();
                for (auto sibling : children)
                    best.push(sibling);
                children.clear();
                return (node);
            }
            auto node = best.top();
            best.pop();
            return (node);
        }

        [[nodiscard]] bool isLIFO(void) const override { return (true); }

    private:
        std::priority_queue<ExploringNode *, std::vector<ExploringNode *>, Cmp>
            best;

        std::vector<ExploringNode *> children; // pushed since the last pop()

    }; // end( class( DiveOpenList ) )

} // anonymous namespace

/*--------------------------------------------------------------------------*/
/*----------------- METHODS OF BranchAndXSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

double BranchAndXSolver::no_incumbent(void) const
{
    return (f_Block->get_objective_sense() == Objective::eMax
                ? -Inf<double>()
                : Inf<double>());
}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue BranchAndXSolver::get_lb(void)
{
    if (f_Block->get_objective_sense() == Objective::eMax)
        return (bestBound); // the incumbent is a lower bound
    // the dual bound is only proven once the whole tree has been explored
    return (f_state == kOK ? bestBound : -Inf<OFValue>());
}

/*--------------------------------------------------------------------------*/

Solver::OFValue BranchAndXSolver::get_ub(void)
{
    if (f_Block->get_objective_sense() == Objective::eMin)
        return (bestBound); // the incumbent is an upper bound
    // the dual bound is only proven once the whole tree has been explored
    return (f_state == kOK ? bestBound : Inf<OFValue>());
}

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

int BranchAndXSolver::compute(bool changedvars)
{
    lock();
    process_outstanding_Modification();
    int old_state = f_state;
    if (f_state == kStillRunning)
        return (kError);
    f_state = kStillRunning;
    std::mutex globalMutex;

    // per-solve budgets and tolerances from the inherited standard parameters
    nodeBudget = get_int_par(intMaxNodes);
    timeBudget = get_dbl_par(dblMaxTime);
    relTol = get_dbl_par(dblRelAcc);
    absTol = get_dbl_par(dblAbsAcc);

    // the tree is retained by any serial exploration, but not by the parallel
    // depth-first one, whose worker forest is not kept [see ParallelDFSSolve()]
    const bool parallelDFS = (solveType == DFS) &&
                             (get_int_par(intMaxThread) > 1);

    // incumbent-dependent local fixing folded into the branching Changes is
    // unsafe when the tree is retained across re-solves with different
    // incumbents: forbid it in that case, allow it otherwise [see
    // GlobalInformation::local_fixing_allowed()]
    f_lfaCell->store(!((reoptimize > 0) && (!parallelDFS)));

    if (changes == 0) // nothing changed since the last solve
        f_state = old_state;
    else
    {
        // the outstanding changes are now properly classified (see
        // process_outstanding_Modification(): 1 = objective only, 2 = feasible
        // region only, 3 = both, 4 = anything); classes 1 / 2 / 3 reoptimize out
        // of the fenced frontier of the previous tree (e.g., an objective-only
        // change leaves every infeasibility certificate valid) rather than solving
        // from scratch, but the tree retained for reoptimization (if any) can only
        // be reused by a serial re-solve under class 1-3 changes when retention is
        // enabled [see intReoptimize / reseedFrontier()]; otherwise it is discarded
        if ((changes == 4) || (!reoptimize) || parallelDFS)
            discardRetainedTree();

        bestBound = no_incumbent();
        // no incumbent yet [see no_incumbent() and GlobalInformation::str_Incumbent]
        f_incumbentCell->store(bestBound);
        // one exploration, parameterized by the open set that the strategy asks
        // for [see treeSolve()]; the only separate one is the parallel depth-first,
        // which runs its own workers [see ParallelDFSSolve()]
        f_state = parallelDFS ? ParallelDFSSolve(globalMutex,
                                                 get_int_par(intMaxThread))
                              : treeSolve(globalMutex);

        // an exploration that has been completed without ever finding a feasible
        // solution proves that there is none
        if ((f_state == kOK) && (bestBound == no_incumbent()))
            f_state = kInfeasible;
    }

    changes = 0;
    unlock();
    return (f_state);

} // end( BranchAndXSolver::compute )

/*--------------------------------------------------------------------------*/

void BranchAndXSolver::discardTree(ExploringNode *root, OpenList &open,
                                   std::list<ChangeSolver *> *solvers)
{
    // the nodes still in the open set are also children in the tree, so the
    // recursive deletion of the tree covers them; emptying the open set first
    // (without deleting) avoids a double delete
    while (!open.empty())
        open.pop();
    std::function<void(ExploringNode *)> deleteTree =
        [&](ExploringNode *node)
    {
        if (!node)
            return;
        for (auto child : node->get_children())
            deleteTree(child);
        node->get_children().clear();
        delete node;
    };
    deleteTree(root);
    delete solvers;
    f_treeRoot = nullptr;
    subOptimalNodes.clear();
    infeasibleNodes.clear();
    integerNodes.clear();
}

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::reseedFrontier(OpenList &open,
                                     std::list<ChangeSolver *> *solvers,
                                     bool minimizing,
                                     ExploringNode *&rootNode,
                                     ExploringNode *&currentNode)
{
    // the retained tree is taken back: compute() guarantees it is only reused
    // for class 1-3 changes [see RelaxationSolver::classify()], so its interior
    // is NOT re-derived, only its fenced frontier is re-evaluated
    rootNode = f_treeRoot;
    f_treeRoot = nullptr; // ownership back to this solve
    currentNode = rootNode;

    std::list<ExploringNode *> frontier;
    frontier.swap(subOptimalNodes);
    if (changes != RelaxationSolver::eModObjective)
    {
        // an objective-only change cannot un-fence an infeasible node: in that
        // case that part of the frontier stays fenced, with no re-evaluation at
        // all; under any other change it goes back into the frontier
        frontier.splice(frontier.end(), infeasibleNodes);
        infeasibleNodes.clear();
    }

    // best (previous) bound first: the optimum almost surely lives in the first
    // few nodes, so the incumbent warms up immediately and the rest of the
    // frontier mostly just re-fences
    frontier.sort([minimizing](ExploringNode *a, ExploringNode *b)
                  { return (minimizing ? a->get_dual_bound() < b->get_dual_bound()
                                       : a->get_dual_bound() > b->get_dual_bound()); });

    int res = Solver::kOK;
    std::vector<ExploringNode *> reopened; // pushed after the loop
    for (auto F : frontier)
    {
        ExploringNode::moveBetweenNodes(currentNode, F, solvers);
        currentNode = F;
        F->initializeBound(minimizing);
        bool toPrune = false;
        RelaxationSolver *nodeBranchSolver = nullptr;
        res = computeRelaxations(&f_RelaxationSolvers, F, minimizing, bestBound,
                                 bestSolution, f_incumbentCell, toPrune,
                                 nodeBranchSolver,
                                 relTol, absTol, maxThreadForSolvers);
        if ((res == Solver::kOK) && (!toPrune))
            res = computeHeuristic(&f_HeuristicSolvers, F, minimizing, bestBound,
                                   bestSolution, f_incumbentCell, toPrune,
                                   maxThreadForSolvers);
        if (res != Solver::kOK)
            break;
        if ((!toPrune) &&
            (!cannot_improve(F->get_dual_bound(), bestBound, minimizing,
                             relTol, absTol)))
        {
            // re-opened: the branching Changes of the previous solve, if any, are
            // still valid (any branching is), so they are reused rather than leaked
            if (F->getBranches().empty())
                F->obtainBranchList(nodeBranchSolver);
            reopened.push_back(F);
        }
        else if (F->is_infeasible()) // fenced again, by reason
            infeasibleNodes.push_back(F);
        else
            subOptimalNodes.push_back(F);
    }

    if (res != Solver::kOK)
    { // leave the solvers at the root
        ExploringNode::moveBetweenNodes(currentNode, rootNode, solvers);
        currentNode = rootNode;
        return (res);
    }

    // hand the re-opened nodes to the open set so that the most promising is
    // explored first whatever the discipline: a LIFO stack receives them in
    // reverse, the others in bound order [see explore()]
    if (open.isLIFO())
        for (auto it = reopened.rbegin(); it != reopened.rend(); ++it)
            open.push(*it);
    else
        for (auto n : reopened)
            open.push(n);

    return (Solver::kOK);

} // end( BranchAndXSolver::reseedFrontier )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::explore(OpenList &open, std::mutex &globalMutex,
                              std::list<ChangeSolver *> *solvers,
                              bool minimizing, int &counter,
                              ExploringNode *rootNode,
                              ExploringNode *currentNode, bool retain,
                              std::chrono::high_resolution_clock::time_point
                                  start)
{
    RelaxationSolver *branchSolver = nullptr;
    int res = Solver::kOK;
    ExploringNode *oldNode = nullptr;
    int iterations = 0; // extractions from the open set, for the log

    while ((!open.empty()) && (nodeBudget > 1) &&
           (timeBudget > std::chrono::duration<double>(
                             std::chrono::high_resolution_clock::now() - start)
                             .count()))
    {
        nodeBudget--;
        ++iterations;
        oldNode = currentNode;
        currentNode = open.pop();
        if (currentNode->get_parent()) // the root: the solvers are already there
            ExploringNode::moveBetweenNodes(oldNode, currentNode, solvers);

        // with the lazy protocol the node is evaluated now, upon extraction,
        // rather than when it was created [see intBoundingProtocol]; a node is
        // evaluated at most once, whichever the protocol
        bool prunedHere = false;
        if (!currentNode->is_evaluated())
        {
            auto evalStart = std::chrono::high_resolution_clock::now();
            res = evaluateNode(currentNode, *solvers, &f_RelaxationSolvers,
                               &f_HeuristicSolvers, minimizing, bestBound,
                               bestSolution, f_incumbentCell, prunedHere,
                               branchSolver,
                               maxThreadForSolvers, relTol, absTol, globalMutex,
                               nullptr, false);
            if (res != Solver::kOK)
            { // an infeasible root ends up here
                discardTree(rootNode, open, solvers);
                return (res);
            }
            if (!prunedHere)
                currentNode->obtainBranchList(branchSolver);
            logNode(f_log, iterations, currentNode, bestBound,
                    std::chrono::duration<double>(
                        std::chrono::high_resolution_clock::now() - start)
                        .count(),
                    std::chrono::duration<double>(
                        std::chrono::high_resolution_clock::now() - evalStart)
                        .count(),
                    prunedHere || cannot_improve(currentNode->get_dual_bound(),
                                                 bestBound, minimizing, relTol,
                                                 absTol));
        }

        // branching and evaluation of the new children
        if ((!prunedHere) &&
            (!cannot_improve(currentNode->get_dual_bound(), bestBound,
                             minimizing, relTol, absTol)))
        {
            auto &branches = currentNode->getBranches();
            std::vector<ExploringNode *> kept; // survivors, pushed after the loop
            for (auto br : branches)
            {
                ExploringNode *new_node = new ExploringNode(br, currentNode,
                                                            currentNode->get_level() + 1, ++counter);
                if (boundingProtocol == Lazy)
                {
                    // the child is not evaluated now: it goes into the open set carrying
                    // the (valid, if weaker) dual bound of its parent, and will be
                    // evaluated if and when it is extracted [see intBoundingProtocol]
                    new_node->set_dual_bound(currentNode->get_dual_bound());
                    currentNode->get_children().push_back(new_node);
                    // the child now owns the branching Change as its f_change: null the
                    // entry so that ~Node does not double-delete it on teardown
                    *(std::find(branches.begin(), branches.end(),
                                new_node->get_f_change())) = nullptr;
                    kept.push_back(new_node);
                    continue;
                }
                bool toPrune = false;
                auto evalStart = std::chrono::high_resolution_clock::now();
                res = evaluateNode(new_node, *solvers, &f_RelaxationSolvers,
                                   &f_HeuristicSolvers, minimizing, bestBound,
                                   bestSolution, f_incumbentCell, toPrune,
                                   branchSolver,
                                   maxThreadForSolvers, relTol, absTol, globalMutex);
                if (res != Solver::kOK)
                {
                    discardTree(rootNode, open, solvers);
                    return (res);
                }
                logNode(f_log, iterations, new_node, bestBound,
                        std::chrono::duration<double>(
                            std::chrono::high_resolution_clock::now() - start)
                            .count(),
                        std::chrono::duration<double>(
                            std::chrono::high_resolution_clock::now() - evalStart)
                            .count(),
                        toPrune || cannot_improve(new_node->get_dual_bound(),
                                                  bestBound, minimizing, relTol,
                                                  absTol));
                // keep the node only if its dual bound can improve the incumbent
                if ((!toPrune) &&
                    (!cannot_improve(new_node->get_dual_bound(), bestBound,
                                     minimizing, relTol, absTol)))
                {
                    new_node->obtainBranchList(branchSolver);
                    currentNode->get_children().push_back(new_node);
                    // the kept child now owns the branching Change as its f_change: null
                    // the entry so that ~Node does not double-delete it on teardown
                    *(std::find(branches.begin(), branches.end(),
                                new_node->get_f_change())) = nullptr;
                    moveSolverToFather(new_node, solvers);
                    kept.push_back(new_node); // pushed to the open set after the loop
                }
                else
                { // fenced or pruned child
                    moveSolverToFather(new_node, solvers);
                    *(std::find(branches.begin(), branches.end(),
                                new_node->get_f_change())) = nullptr;
                    if (retain)
                    { // it belongs to the fenced frontier:
                        currentNode->get_children().push_back(new_node);
                        if (new_node->is_infeasible())           // infeasibility certificates survive
                            infeasibleNodes.push_back(new_node); //  objective-only changes
                        else
                            subOptimalNodes.push_back(new_node);
                    }
                    else
                        delete new_node;
                }
            }
            // hand the surviving children to the open set so that the most promising
            // (the first from branch()) is explored first whatever the discipline:
            // a LIFO stack receives them in reverse, the others in branching order
            if (open.isLIFO())
                for (auto it = kept.rbegin(); it != kept.rend(); ++it)
                    open.push(*it);
            else
                for (auto n : kept)
                    open.push(n);
            branches.erase(std::remove(branches.begin(), branches.end(),
                                       nullptr),
                           branches.end());
            if (currentNode->get_children().empty())
            {
                if (retain) // fathomed leaf: fenced frontier
                    subOptimalNodes.push_back(currentNode);
                if (currentNode->get_toFather())
                { // the root has nowhere to climb
                    if (retain)
                    {
                        for (const auto s : *solvers)
                            s->apply(currentNode->get_toFather(), false);
                        currentNode = currentNode->get_parent();
                    }
                    else
                        currentNode = ExploringNode::prune(currentNode, solvers);
                }
            }
        }
        else
        {
            // fenced at pop: it belongs to the fenced frontier, the root included -
            // a retained tree whose root is fathomed would otherwise have an empty
            // frontier, and the next re-solve would have nothing to re-explore
            if (retain)
                subOptimalNodes.push_back(currentNode);
            if (currentNode->get_toFather())
            { // the root has nowhere to climb
                if (retain)
                {
                    for (const auto s : *solvers)
                        s->apply(currentNode->get_toFather(), false);
                    currentNode = currentNode->get_parent();
                }
                else
                    currentNode = ExploringNode::prune(currentNode, solvers);
            }
        }
    }

    // move the :ChangeSolver back to the root
    while (currentNode->get_toFather())
    {
        for (const auto s : *solvers)
            s->apply(currentNode->get_toFather(), false);
        currentNode = currentNode->get_parent();
    }

    if (retain)
    {
        // retain the tree for future reoptimizations: the nodes still in the open
        // set (early stops) are open work, hence part of the frontier to re-seed
        while (!open.empty())
            subOptimalNodes.push_back(open.pop());
        f_treeRoot = rootNode;
    }
    else
    {
        // the nodes still in the open set are also children in the tree: the
        // recursive deletion of the tree covers them (deleting them from the open
        // set too would be a double delete)
        while (!open.empty())
            open.pop();
        std::function<void(ExploringNode *)> deleteTree =
            [&](ExploringNode *node)
        {
            if (!node)
                return;
            for (auto child : node->get_children())
                deleteTree(child);
            node->get_children().clear();
            delete node;
        };
        deleteTree(rootNode);
    }
    delete solvers;

    if (nodeBudget <= 0)
        return (Solver::kStopIter);
    if (timeBudget <= std::chrono::duration<double>(
                          std::chrono::high_resolution_clock::now() - start)
                          .count())
        return (Solver::kStopTime);
    return (kOK);

} // end( BranchAndXSolver::explore )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::treeSolve(std::mutex &globalMutex)
{
    auto *solvers = new std::list<ChangeSolver *>();
    bool minimizing;
    initializeVariables(solvers, &f_RelaxationSolvers, &f_HeuristicSolvers,
                        minimizing);

    auto start = std::chrono::high_resolution_clock::now();
    int counter = 0;

    // the open set *is* the exploration strategy: a LIFO stack explores
    // depth-first, a FIFO queue breadth-first, a dual-bound priority queue
    // best-first, and the dive variant wraps the latter so that each best node
    // is followed depth-first down to a leaf [see DiveOpenList]
    auto cmp = [minimizing](ExploringNode *a, ExploringNode *b)
    {
        return (minimizing ? a->get_dual_bound() > b->get_dual_bound()
                           : a->get_dual_bound() < b->get_dual_bound());
    };
    std::unique_ptr<OpenList> openPtr;
    switch (solveType)
    {
    case (DFS):
        openPtr.reset(new StackOpenList());
        break;
    case (BFS):
        openPtr.reset(new QueueOpenList());
        break;
    case (BestFS):
        openPtr.reset(new PriorityOpenList(cmp));
        break;
    case (BestFSDive):
        openPtr.reset(new DiveOpenList(cmp));
        break;
    default:
        delete solvers;
        throw(std::invalid_argument("BranchAndXSolver::treeSolve: invalid "
                                    "intSolveMethod"));
    }
    OpenList &open = *openPtr;

    logHeader(f_log);

    ExploringNode *rootNode = new ExploringNode(nullptr, nullptr, 0);
    RelaxationSolver *branchSolver = nullptr;
    // the tree is retained for reoptimization whatever the strategy [see
    // intReoptimize]: a re-solve re-seeds from the frontier of the retained
    // tree, in the order that the open set of this very strategy dictates
    const bool retain = (reoptimize > 0);
    int res;
    ExploringNode *currentNode;
    if (f_treeRoot)
    {
        delete rootNode; // the fresh root is not needed
        res = reseedFrontier(open, solvers, minimizing, rootNode, currentNode);
    }
    else
    {
        if (boundingProtocol == Lazy)
        {
            // like any other node, the root is evaluated when it is extracted
            rootNode->initializeBound(minimizing);
            res = Solver::kOK;
        }
        else
            res = initializeRoot(&f_RelaxationSolvers, &f_HeuristicSolvers,
                                 rootNode, minimizing, bestBound, bestSolution,
                                 f_incumbentCell, branchSolver, globalMutex);
        if (res == Solver::kOK)
            open.push(rootNode);
        currentNode = rootNode;
    }

    if (res != Solver::kOK)
    {
        discardTree(rootNode, open, solvers);
        return (res);
    }

    return (explore(open, globalMutex, solvers, minimizing, counter,
                    rootNode, currentNode, retain, start));

} // end( BranchAndXSolver::treeSolve )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

int BranchAndXSolver::workerDFS(Node *currentNode,
                                std::list<ChangeSolver *> &solvers,
                                std::vector<std::pair<RelaxationSolver *,
                                                      Solver *>> &relaxation,
                                std::vector<std::pair<ChangeSolver *,
                                                      Solver *>> &heuristic,
                                bool minimizing,
                                std::mutex &incumbentMutex,
                                std::atomic<int> &nodeBdg,
                                std::chrono::steady_clock::time_point
                                    deadline,
                                int &nameCounter)
{
    if (--nodeBdg <= 0)
        return (Solver::kStopIter);
    if (std::chrono::steady_clock::now() > deadline)
        return (Solver::kStopTime);

    RelaxationSolver *branchSolver = nullptr;
    bool toPrune = false;
    // out-of-mutex reads of the shared bestBound may be slightly stale, only
    // making the pruning marginally less aggressive; every incumbent update
    // happens under incumbentMutex with a double check. The solvers run serially
    // (the per-node-solvers parallelism does not lock the incumbent), so the
    // globalMutex argument is unused here
    int res = evaluateNode(currentNode, solvers, &relaxation, &heuristic,
                           minimizing, bestBound, bestSolution,
                           f_incumbentCell, toPrune,
                           branchSolver, 1,
                           relTol, absTol, incumbentMutex, &incumbentMutex);
    if (res != ThinComputeInterface::kOK)
        return (res);

    if ((!toPrune) &&
        (!cannot_improve(currentNode->get_dual_bound(), bestBound,
                         minimizing, relTol, absTol)))
    {
        auto branches = branchSolver->branch();
        for (auto br : branches)
        {
            DFSNode *new_node = new DFSNode(br, ++nameCounter);
            int RV = workerDFS(new_node, solvers, relaxation, heuristic,
                               minimizing, incumbentMutex, nodeBdg, deadline,
                               nameCounter);
            delete new_node;
            if (RV != Solver::kOK)
                return (RV);
        }
    }

    moveSolverToFather(currentNode, &solvers);
    return (Solver::kOK);

} // end( BranchAndXSolver::workerDFS )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::ParallelDFSSolve(std::mutex &globalMutex, int K)
{
    createWorkerSolvers(K);

    // the serial Solver set drives the ramp-up
    std::list<ChangeSolver *> solvers;
    bool minimizing;
    initializeVariables(&solvers, &f_RelaxationSolvers, &f_HeuristicSolvers,
                        minimizing);

    // an infinite time budget must not be duration_cast (it overflows into a
    // deadline in the past): it simply means no deadline at all
    const auto deadline = timeBudget == Inf<double>() ? std::chrono::steady_clock::time_point::max() : std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(timeBudget));
    std::atomic<int> nodeBdg(nodeBudget);

    // ramp-up- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // ordered (FIFO = discovery order) BestFS-style expansion with the serial
    // Solver set: every node put in the pool has its dual bound, branching
    // Changes (see obtainBranchList()) and undo (toFather) already computed,
    // so the workers can later claim it by just re-applying the Changes found
    // on its path, which are plain data usable by any Solver
    int nameCounter = 0;
    ExploringNode *rootNode = new ExploringNode(nullptr, nullptr, 0);
    RelaxationSolver *branchSolver = nullptr;
    int res = initializeRoot(&f_RelaxationSolvers, &f_HeuristicSolvers,
                             rootNode, minimizing, bestBound, bestSolution,
                             f_incumbentCell, branchSolver, globalMutex);
    if (res != Solver::kOK)
    {
        delete rootNode;
        return (res);
    }

    std::deque<ExploringNode *> pool;
    pool.push_back(rootNode);
    ExploringNode *currentNode = rootNode;
    const std::size_t target = std::size_t(4 * K);

    while ((!pool.empty()) && (pool.size() < target) &&
           (--nodeBdg > 0) &&
           (std::chrono::steady_clock::now() < deadline))
    {
        ExploringNode *oldNode = currentNode;
        currentNode = pool.front();
        pool.pop_front();
        if (currentNode->get_parent())
            ExploringNode::moveBetweenNodes(oldNode, currentNode, &solvers);

        if (!cannot_improve(currentNode->get_dual_bound(), bestBound,
                            minimizing, relTol, absTol))
        {
            auto &branches = currentNode->getBranches();
            for (Change *br : branches)
            {
                auto new_node = new ExploringNode(br, currentNode,
                                                  currentNode->get_level() + 1,
                                                  ++nameCounter);
                bool toPrune = false;
                RelaxationSolver *nodeBranchSolver = nullptr;
                res = evaluateNode(new_node, solvers, &f_RelaxationSolvers,
                                   &f_HeuristicSolvers, minimizing, bestBound,
                                   bestSolution, f_incumbentCell, toPrune,
                                   nodeBranchSolver,
                                   1, relTol, absTol, globalMutex);
                if (res != Solver::kOK)
                {
                    moveSolverToFather(new_node, &solvers);
                    delete new_node;
                    break;
                }
                // either way the ownership of br leaves the branches of currentNode
                // (the kept child owns it as its f_change, the discarded one died
                // with it): null the entry so that ~Node does not double-delete it
                *(std::find(branches.begin(), branches.end(),
                            new_node->get_f_change())) = nullptr;
                if ((!toPrune) &&
                    (!cannot_improve(new_node->get_dual_bound(), bestBound,
                                     minimizing, relTol, absTol)))
                {
                    new_node->obtainBranchList(nodeBranchSolver);
                    pool.push_back(new_node);
                    currentNode->get_children().push_back(new_node);
                    moveSolverToFather(new_node, &solvers);
                }
                else
                { // fenced or pruned: discard the child
                    moveSolverToFather(new_node, &solvers);
                    delete new_node;
                }
            }
            if (res != Solver::kOK)
                break;
        }
    }

    // back to the root, ready for the workers
    ExploringNode::moveBetweenNodes(currentNode, rootNode, &solvers);

    // workers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // each worker repeatedly claims the OLDEST open subtree of the pool
    // (preserving the sequential search order), positions its own Solver set
    // on it by re-applying the Changes of its path, explores it depth-first
    // and un-winds back to the root via the (already computed) toFather
    std::atomic<int> result(res);
    if ((res == Solver::kOK) && (!pool.empty()))
    {
        std::mutex poolMutex, incumbentMutex;
        std::atomic<std::size_t> poolSize(pool.size());
        std::atomic<int> busy(0);

        auto workerLoop = [&](int w)
        {
            auto &ws = v_workerSolvers[w];
            std::list<ChangeSolver *> wsolvers;
            for (auto &pr : ws.relaxation)
                wsolvers.push_back(pr.first);
            for (auto &pr : ws.heuristic)
                wsolvers.push_back(pr.first);
            int wNameCounter = (w + 1) * 10000000; // disjoint per-worker names

            while (result.load() == Solver::kOK)
            {
                ExploringNode *node;
                {
                    std::lock_guard<std::mutex> guard(poolMutex);
                    if (pool.empty())
                    {
                        // the pool may be re-fed by the work-sharing of a busy worker: only
                        // an empty pool with NO busy worker means there is no work left
                        if (busy.load() == 0)
                            break;
                        node = nullptr;
                    }
                    else
                    {
                        node = pool.front();
                        pool.pop_front();
                        --poolSize;
                        ++busy;
                    }
                }
                if (!node)
                { // pool empty but someone is still working
                    std::this_thread::sleep_for(std::chrono::microseconds(100));
                    continue;
                }

                // claim the subtree: re-apply the path of Changes from the root
                std::list<Change *> path;
                for (auto n = node; n; n = n->get_parent())
                    if (n->get_f_change())
                        path.push_front(n->get_f_change());
                for (auto chg : path)
                    for (auto s : wsolvers)
                        s->apply(chg, false);

                // explore the subtree depth-first, unless meanwhile fenced
                if (!cannot_improve(node->get_dual_bound(), bestBound, minimizing,
                                    relTol, absTol))
                {
                    bool first = true;
                    for (auto &br : node->getBranches())
                    {
                        if (!br) // child discarded during the ramp-up
                            continue;

                        // work-sharing: if the pool is REALLY starving (less than half the
                        // workers could claim something), donate this branch to it (eagerly
                        // evaluated like in the ramp-up, so it is claimable) rather than
                        // exploring it; the donated child hangs off the claimed node, which
                        // always outlives it, so the path stays replayable. The first branch
                        // is always explored locally, and the threshold is conservative: on
                        // small subtrees an aggressive donation policy costs more (one eager
                        // evaluation per donation) than it parallelizes
                        if ((!first) && (2 * poolSize.load() < std::size_t(K)))
                        {
                            auto new_node = new ExploringNode(br, node,
                                                              node->get_level() + 1,
                                                              ++wNameCounter);
                            br = nullptr;
                            moveSolverToSon(new_node, &wsolvers);
                            new_node->initializeBound(minimizing);
                            bool toPrune = false;
                            RelaxationSolver *nodeBranchSolver = nullptr;
                            int RV = computeRelaxations(&ws.relaxation, new_node, minimizing,
                                                        bestBound, bestSolution,
                                                        f_incumbentCell, toPrune,
                                                        nodeBranchSolver, relTol, absTol,
                                                        1, &incumbentMutex);
                            if (RV == Solver::kOK && !toPrune)
                                RV = computeHeuristic(&ws.heuristic, new_node, minimizing,
                                                      bestBound, bestSolution,
                                                      f_incumbentCell, toPrune,
                                                      1, &incumbentMutex);
                            moveSolverToFather(new_node, &wsolvers);
                            if (RV != Solver::kOK)
                            {
                                delete new_node;
                                int expected = Solver::kOK;
                                result.compare_exchange_strong(expected, RV);
                                break;
                            }
                            if ((!toPrune) &&
                                (!cannot_improve(new_node->get_dual_bound(), bestBound,
                                                 minimizing, relTol, absTol)))
                            {
                                new_node->obtainBranchList(nodeBranchSolver);
                                node->get_children().push_back(new_node);
                                std::lock_guard<std::mutex> guard(poolMutex);
                                pool.push_back(new_node);
                                ++poolSize;
                            }
                            else // pruned or fenced: nothing to donate
                                delete new_node;
                            continue;
                        }

                        // the DFSNode takes the ownership of the branching Change: null the
                        // entry so that ~Node does not double-delete it on the final cleanup
                        DFSNode *new_node = new DFSNode(br, ++wNameCounter);
                        br = nullptr;
                        first = false;
                        int RV = workerDFS(new_node, wsolvers, ws.relaxation,
                                           ws.heuristic, minimizing, incumbentMutex,
                                           nodeBdg, deadline, wNameCounter);
                        delete new_node;
                        if (RV != Solver::kOK)
                        {
                            int expected = Solver::kOK;
                            result.compare_exchange_strong(expected, RV);
                            break;
                        }
                    }
                }

                // un-wind back to the root via the toFather of the path
                for (auto n = node; n; n = n->get_parent())
                    if (n->get_toFather())
                        for (auto s : wsolvers)
                            s->apply(n->get_toFather(), false);

                --busy;
            }
        };

        std::vector<std::thread> workers;
        workers.reserve(K);
        for (int w = 0; w < K; ++w)
            workers.emplace_back(workerLoop, w);
        for (auto &t : workers)
            t.join();
    }

    // cleanup of the ramp-up tree- - - - - - - - - - - - - - - - - - - - - - - -
    std::function<void(ExploringNode *)> deleteTree =
        [&](ExploringNode *node)
    {
        for (auto child : node->get_children())
            deleteTree(child);
        node->get_children().clear();
        delete node;
    };
    deleteTree(rootNode);

    if ((result.load() == Solver::kOK) && (nodeBdg <= 0))
        return (Solver::kStopIter);
    return (result.load());

} // end( BranchAndXSolver::ParallelDFSSolve )

/*--------------------------------------------------------------------------*/

void BranchAndXSolver::discardRetainedTree(void)
{
    if (!f_treeRoot)
        return;

    std::function<void(ExploringNode *)> deleteTree =
        [&](ExploringNode *node)
    {
        for (auto child : node->get_children())
            deleteTree(child);
        node->get_children().clear();
        delete node;
    };
    deleteTree(f_treeRoot);

    f_treeRoot = nullptr;
    subOptimalNodes.clear();
    infeasibleNodes.clear();
    integerNodes.clear();

} // end( BranchAndXSolver::discardRetainedTree )

/*--------------------------------------------------------------------------*/

void BranchAndXSolver::process_outstanding_Modification(void)
{
    Lst_sp_Mod v_mod_tmp; // temporary list of modifications

    // try to acquire lock, spin on failure
    while (f_mod_lock.test_and_set(std::memory_order_acquire))
        ;

    for (auto mod : v_mod)
        v_mod_tmp.push_back(mod); // copy v_mod in v_mod_tmp

    v_mod.clear();

    f_mod_lock.clear(std::memory_order_release); // release lock

    if (v_mod_tmp.empty())
        return;

    // classify the changes by asking a RelaxationSolver: what a Modification
    // does to the fencing certificates of the tree is problem-specific
    // knowledge [see RelaxationSolver::classify()], so this Solver never
    // looks into the Modification itself; the classes compose bitwise, and
    // anything the RelaxationSolver cannot vouch for invalidates everything
    if (f_RelaxationSolvers.empty())
    {
        changes = 4;
        return;
    }
    auto rs = f_RelaxationSolvers.front().first;
    int acc = changes == 4 ? RelaxationSolver::eModEverything : changes;
    for (const auto &mod : v_mod_tmp)
    {
        const int c = rs->classify(mod);
        if (c == RelaxationSolver::eModEverything)
        {
            acc = c;
            break;
        }
        acc |= c;
    }
    changes = char(acc);

} // end( BranchAndXSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/
/*------------------ End File BranchAndXSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/

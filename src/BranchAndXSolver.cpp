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
#include <cmath>
#include <deque>
#include <memory>
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

SMSpp_insert_in_factory_cpp_0(BranchAndXSolver);

/*--------------------------------------------------------------------------*/
/*-------------------------- INTERNAL FUNCTIONS ----------------------------*/
/*--------------------------------------------------------------------------*/

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
/* template <typename SolverPtr>
static bool computeAllParallel(const std::vector<SolverPtr> &slvrs,
                               std::vector<int> &zs, int maxThreads);
 */
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
    std::vector<RelaxationSolver *> *f_RelaxationSolvers,
    Node *currentNode,
    const bool minimizing,
    double &bestBound,
    Solution *&bestSol,
    bool &toPrune,
    RelaxationSolver *&branchSolver,
    BranchAndXSolver *tree //,
    // int nThreads = 1, std::mutex *incumbentMutex = nullptr
)
{
    for (auto s : *f_RelaxationSolvers)
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

        // see if the primal bound improves thanks to a true solution
        if (s->has_true_var_solution())
        {
            double primal_bound = minimizing ? s->get_true_ub() : s->get_true_lb();
            /*             if (minimizing ? primal_bound < bestBound : primal_bound > bestBound)
                        {
                            // in the parallel exploration the incumbent is shared between the
                            // workers: re-check the improvement under the mutex
                            std::unique_lock<std::mutex> guard;
                            if (incumbentMutex)
                                guard = std::unique_lock<std::mutex>(*incumbentMutex);
             */
            if (minimizing ? primal_bound < bestBound : primal_bound > bestBound)
            {
                bestBound = primal_bound;
                delete bestSol;
                bestSol = s->get_true_solution();
                tree->globalInfoWrite("incumbent", "bestBound", bestBound);
            }
        }

        // update the dual bound of the node and select the branching solver
        auto dualBound = minimizing ? s->get_lb() : s->get_ub();
        if (tree->cannot_improve(dualBound, bestBound, minimizing))
        {
            toPrune = true;
            return (Solver::kOK);
        }
        /*         if (minimizing ? dualBound > currentNode->get_dual_bound()
                               : dualBound < currentNode->get_dual_bound())
                {
                    std::unique_lock<std::mutex> guard;
                    if (incumbentMutex)
                        guard = std::unique_lock<std::mutex>(*incumbentMutex);
         */
        if (minimizing ? dualBound > currentNode->get_dual_bound()
                       : dualBound < currentNode->get_dual_bound())
        {
            currentNode->set_dual_bound(dualBound);
            branchSolver = s;
            //}
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
    std::vector<ChangeSolver *> *f_HeuristicSolvers,
    Node *currentNode, const bool minimizing,
    double &bestBound, Solution *&bestSol, bool &toPrune,
    BranchAndXSolver *tree
    // int nThreads = 1, std::mutex *incumbentMutex = nullptr
)
{
    for (auto s : *f_HeuristicSolvers)
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
            /*             std::unique_lock<std::mutex> guard;
                        if (incumbentMutex)
                            guard = std::unique_lock<std::mutex>(*incumbentMutex);
                        if (minimizing ? primal_bound < bestBound : primal_bound > bestBound)
                        {
             */
            bestBound = primal_bound;
            delete bestSol;
            bestSol = s->get_Solution();
            tree->globalInfoWrite("incumbent", "bestBound", bestBound);
            //}
        }
    }
    return (Solver::kOK);
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

    private:
        std::priority_queue<ExploringNode *, std::vector<ExploringNode *>, Cmp>
            best;

        std::vector<ExploringNode *> children; // pushed since the last pop()

    }; // end( class( DiveOpenList ) )

} // anonymous namespace

/*--------------------------------------------------------------------------*/
/*----------------- METHODS OF BranchAndXSolver ------------------------*/
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
    // TODO capire come fare quando voglio "infiniti" nodi (o finché non esplode la macchina)
    nodeBudget = get_int_par(intMaxNodes);
    timeBudget = get_dbl_par(dblMaxTime);
    /* 	relTol = get_dbl_par(dblRelAcc);
        absTol = get_dbl_par(dblAbsAcc); */

    // incumbent-dependent local fixing folded into the branching Changes is
    // unsafe when the tree is retained across re-solves with different
    // incumbents: forbid it in that case, allow it otherwise [see
    // GlobalInformation::local_fixing_allowed()]

    /* 	f_globalInfo.set_local_fixing_allowed(
            !((reoptimize > 0) && (solveType == BestFS)));
     */
    if (changes == 0) // nothing changed since the last solve
        f_state = old_state;
    else
    {

        f_globalInfo.add_to_Universe<double>("DoubleProperties");
        bestBound = (f_Block->get_objective_sense() == Objective::eMax)
                        ? -Inf<double>()
                        : Inf<double>();
        auto doubleValues = f_globalInfo.get_from_Universe<double>("DoubleProperties");
        if (doubleValues)
            doubleValues->write("incumbent", bestBound);
        f_state = treeSolve(globalMutex);
    }
    changes = 0;
    unlock();
    return (f_state);

} // end( BranchAndXSolver::compute )

/*--------------------------------------------------------------------------*/

int BranchAndXSolver::explore(OpenList &open, std::mutex &globalMutex,
                              std::list<ChangeSolver *> *solvers,
                              bool minimizing, int &counter,
                              ExploringNode *rootNode,
                              ExploringNode *currentNode, // bool retain,
                              std::chrono::high_resolution_clock::time_point
                                  start)
{
    RelaxationSolver *branchSolver = nullptr;
    int res = Solver::kOK;
    ExploringNode *oldNode = nullptr;

    while ((!open.empty()) && (nodeBudget > 1) &&
           (timeBudget > std::chrono::duration<double>(
                             std::chrono::high_resolution_clock::now() - start)
                             .count()))
    {
        if (nodeBudget != INT_MAX)
            nodeBudget--;
        oldNode = currentNode;
        currentNode = open.pop();
        if (currentNode->get_parent()) // the root: the solvers are already there
            ExploringNode::moveBetweenNodes(oldNode, currentNode, solvers);
        // Lazy evaluation
        if (this->boundingProtocol == Lazy)
        {
            bool toPrune = false;

            res = computeRelaxations(&f_RelaxationSolvers, currentNode, minimizing,
                                     bestBound, bestSolution, toPrune, branchSolver,
                                     this //, maxThreadForSolvers, &globalMutex
            );
            if (res != Solver::kOK)
            {
                return (res);
            }
            if (toPrune)
            {
                currentNode = ExploringNode::prune(currentNode, solvers);
                continue;
            }
            res = computeHeuristic(&f_HeuristicSolvers, currentNode, minimizing,
                                   bestBound, bestSolution, toPrune,
                                   this
                                   // maxThreadForSolvers, &globalMutex);
            );
            if (res != Solver::kOK)
                return res;
            // keep the node only if its dual bound can improve the incumbent
            if (!cannot_improve(currentNode->get_dual_bound(), bestBound, minimizing))
                currentNode->obtainBranchList(branchSolver);
            else
            { // fenced or pruned child
                currentNode = ExploringNode::prune(currentNode, solvers);
                continue;
            }
        }

        // branching and evaluation of the new children
        if (!cannot_improve(currentNode->get_dual_bound(), bestBound, minimizing))
        {
            auto &branches = currentNode->getBranches();
            // std::vector<ExploringNode *> kept; // survivors, pushed after the loop
            for (auto br : branches)
            {
                ExploringNode *new_node = new ExploringNode(br, currentNode,
                                                            currentNode->get_level() + 1, ++counter);
                if (this->boundingProtocol == Eager)
                {
                    bool toPrune = false;
                    computeRelaxations(&f_RelaxationSolvers, new_node, minimizing,
                                       bestBound, bestSolution, toPrune, branchSolver,
                                       this //,
                                            // maxThreadForSolvers, &globalMutex
                    );
                    if (res != Solver::kOK)
                    {
                        return (res);
                    }
                    if (!toPrune)
                        computeHeuristic(&f_HeuristicSolvers, new_node, minimizing,
                                         bestBound, bestSolution, toPrune,
                                         this //,
                                              // maxThreadForSolvers, &globalMutex
                        );
                    // keep the node only if its dual bound can improve the incumbent
                    if ((!toPrune) &&
                        (!cannot_improve(new_node->get_dual_bound(), bestBound,
                                         minimizing)))
                    {
                        new_node->obtainBranchList(branchSolver);
                        currentNode->get_children().push_back(new_node);
                        // the kept child now owns the branching Change as its f_change: null
                        // the entry so that ~Node does not double-delete it on teardown
                        // TODO controllare che l'idea sia giusta, che non mi pare che lo sia così tanto
                        *(std::find(branches.begin(), branches.end(), new_node->get_f_change())) = nullptr;
                        moveSolverToFather(new_node, solvers);
                        // kept.push_back(new_node); // pushed to the open set after the loop
                        open.push(new_node); // pushed to the open set after the loop
                    }
                    else
                    { // fenced or pruned child
                        moveSolverToFather(new_node, solvers);
                        // TODO capire quanto possa essere costoso e come evitarlo
                        *(std::find(branches.begin(), branches.end(), new_node->get_f_change())) = nullptr;
                        delete new_node;
                    }
                }
                else
                {
                    new_node->set_dual_bound(currentNode->get_dual_bound());
                    open.push(new_node);
                }
            }

            branches.erase(std::remove(branches.begin(), branches.end(), nullptr), branches.end());
            if (currentNode->get_children().empty() && currentNode->get_toFather())
            {

                currentNode = ExploringNode::prune(currentNode, solvers);
            }
        }
        else if (currentNode->get_toFather())
        {

            currentNode = ExploringNode::prune(currentNode, solvers);
        }
    }

    // move the :ChangeSolver back to the root
    while (currentNode->get_toFather())
    {
        for (const auto s : *solvers)
            s->apply(currentNode->get_toFather(), false);
        currentNode = currentNode->get_parent();
    }

    /* 	if (retain)
        {
            // retain the tree for future reoptimizations: the nodes still in the open
            // set (early stops) are open work, hence part of the frontier to re-seed
            while (!open.empty())
                subOptimalNodes.push_back(open.pop());
            f_treeRoot = rootNode;
        }
        else
        { */
    // the nodes still in the open set are also children in the tree: the
    // recursive deletion of the tree covers them (deleting them from the open
    // set too would be a double delete)

    // TODO understand if it's right or we can do something better
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
    //}
    delete solvers;

    // consider also root node that don't decrease nodeBudget
    if (nodeBudget <= 1)
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

    if (!f_HeuristicSolvers.empty())
        minimizing = f_HeuristicSolvers.front()->get_Block()->get_objective_sense() == Objective::eMin;
    else if (!f_RelaxationSolvers.empty())
        minimizing = f_RelaxationSolvers.front()->get_Block()->get_objective_sense() == Objective::eMin;
    else
        throw(std::logic_error("BranchAndXSolver::initializeVariables: "
                               "both the HeuristicSolvers and the RelaxationSolvers are empty"));

    solvers->insert(solvers->end(), f_RelaxationSolvers.begin(),
                    f_RelaxationSolvers.end());
    solvers->insert(solvers->end(), f_HeuristicSolvers.begin(),
                    f_HeuristicSolvers.end());

    auto start = std::chrono::high_resolution_clock::now();
    int counter = 0;
    // the open set is a plain FIFO queue
    OpenList *open = nullptr;
    switch (solveType)
    {
    case (DFS):
    {
        open = new StackOpenList();
        break;
    }
    case (BFS):
    {
        open = new QueueOpenList();
        break;
    }
    case (BestFS):
    {
        auto cmp = [minimizing](ExploringNode *a, ExploringNode *b)
        {
            return (minimizing ? a->get_dual_bound() > b->get_dual_bound()
                               : a->get_dual_bound() < b->get_dual_bound());
        };
        open = new PriorityOpenList(cmp);
        break;
    }
    case (BestFSDive):
    {
        auto cmp = [minimizing](ExploringNode *a, ExploringNode *b)
        {
            return (minimizing ? a->get_dual_bound() > b->get_dual_bound()
                               : a->get_dual_bound() < b->get_dual_bound());
        };
        open = new DiveOpenList(cmp);
        break;
    }
    default:
        throw(std::invalid_argument("BranchAndXSolver::compute: invalid "
                                    "intSolveMethod"));
    }
    ExploringNode *rootNode = new ExploringNode(nullptr, nullptr, 0);
    RelaxationSolver *branchSolver = nullptr;
    // initialize variables
    rootNode->initializeBound(minimizing);
    if (boundingProtocol == Eager)
    {
        bool toPrune = false;
        int res = computeRelaxations(&f_RelaxationSolvers, rootNode, minimizing,
                                     bestBound, bestSolution, toPrune, branchSolver, this); // nThreads, &globalMutex);
        if (toPrune || (res != ThinComputeInterface::kOK))
            return (res);
        res = computeHeuristic(&f_HeuristicSolvers, rootNode, minimizing,
                               bestBound, bestSolution, toPrune, this); // nThreads, &globalMutex);
        if (toPrune || (res != ThinComputeInterface::kOK))
            return (res);
        rootNode->obtainBranchList(branchSolver);
    }
    open->push(rootNode);
    int returnValue = explore(*open, globalMutex, solvers, minimizing, counter,
                              rootNode, rootNode, start);
    delete open;
    return (returnValue);

} // end( BranchAndXSolver::treeSolve )

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

    /*
    if (v_mod_tmp.empty())
        return;
    */

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
    v_mod_tmp.clear(); // clear the temporary list of Modification

} // end( BranchAndXSolver::process_outstanding_Modification )

/*--------------------------------------------------------------------------*/
/*------------------ End File BranchAndXSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/

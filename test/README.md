# test

A tester for `BranchAndXSolver` that needs nothing but the core SMS++ library.

It builds in memory an `AbstractBlock` with a few `ColVariable`, a
`BoxConstraint` on each of them and a separable quadratic `DQuadFunction` as
`Objective`, to be minimized (convex) or maximized (concave) over the integer
points of the box. The integer optimum is known in closed form (in each
coordinate, the integer of the box nearest to the unconstrained optimum), and
the closed form is checked against the enumeration of all the integer points.

## Running the tests

With testing enabled in the CMake build, run these commands from the module
or SMS++ umbrella root, using the path to its build directory:

```sh
cmake --build build --target BranchAndXSolver_unit_test
ctest --test-dir build -R '^BranchAndXSolver_unit_test$' --output-on-failure
```

CTest runs the executable with this source directory as its working
directory, so the `InnerBSCfg*.txt` files are found. When running the executable
directly, use this directory as the working directory as well.

The exit code is 0 when every check passes, printing `All tests passed!!`, and
1 otherwise. The `makefile` also builds the executable with the
`BranchAndXSolver` module and the core SMS++ library.

## Test solvers

The core has no `RelaxationSolver`, so the tester defines the inner `Solver`
that `BranchAndXSolver` drives, all based on `BoxSolver`: `BoxRelaxation`
solves the continuous relaxation (the `ColVariable` are continuous in the
`Block`, their integrality is the knowledge of this `Solver`) and branches on
its most fractional coordinate by tightening the `BoxConstraint`, while
`BoxRounding` is a heuristic `ChangeSolver` rounding the solution of the
relaxation inside the box. `BoxEstimate` is the `BoxRelaxation` that is also a
`BranchEstimator` (`BranchingRule.h` of the core), estimating the best value
under each child of a branching as the value of the rounding in the box of
the child, and `BoxWrongEstimate` swaps those estimates and at times gives
none. `BoxUnresolved` models an advanced partial inner computation at every
third `compute()` and reports `kOK` with its usable bound and branching
information. They are given to `BranchAndXSolver` by the `BlockSolverConfig`
in `InnerBSCfg.txt` (the relaxation alone), `InnerBSCfg-rounding.txt` (the
relaxation and the heuristic), `InnerBSCfg-estimate.txt` and
`InnerBSCfg-wrong-estimate.txt` (the relaxations with the estimates) and
`InnerBSCfg-unresolved.txt` (the advanced partial-result wrapper), which the
tester reads from this directory.

Additional test solvers provide controlled return codes, missing bounds or
solutions, inconsistent certificates, and failures at the root or a child.
They let the regressions check the behavior of `BranchAndXSolver` independently
of a particular external optimization library.

## Coverage

The cases are:

- every exploration strategy (`intSolveMethod`) under both bounding protocols
  (`intBoundingProtocol`), minimizing, and every strategy maximizing, each run
  having to find the optimum, with both bounds equal to it and the optimal
  integer solution written back in the `Block`;
- the best-estimate search (`BestEstimate`) with the estimates of the
  relaxation, right or misleading, under both bounding protocols and in both
  senses: the estimates have to be asked for, the optimum has to be the
  same, and the order of the branchings has to depend on the estimates;
- advanced partial inner computations, with every strategy under both
  bounding protocols, in both senses and with a retained tree: usable
  results must let the exploration reach the optimum, and an unchanged
  repeat solve must reuse the completed result;
- the rounding heuristic, whose incumbent is the solution written back;
- an empty box and a box with no integer point (whose relaxation is not
  empty), both infeasible, with both bounds at the infinity of the sense;
- a root whose relaxation is already integer, solved within a budget of one
  node;
- a node budget (`intMaxNodes`) and a time budget (`dblMaxTime`) stopping the
  enumeration, with the status telling which one and no bound claimed on the
  side of the nodes not explored;
- the parameters: names, defaults, a `ComputeConfig` setting them by name, the
  factory reset, invalid values, a `BlockSolverConfig` replaced on an attached
  `Solver`;
- the `Solver` detached from a `Block`, attached to another one and back;
- a change of the `Objective` between two `compute()`, solving from scratch
  and reoptimizing the retained tree (`intReoptimize`), and a `compute()` with
  nothing changed;
- a change of variable integrality, retaining the tree when valid and
  rebuilding it when necessary, with the result compared against a solve
  from scratch.

### Partial results and branching

`test_partial_status_range` checks the interval `[kOK, kError)` directly:

| Returned status | Expected behavior |
| --- | --- |
| `kOK` | Continue the search with the relaxation result. |
| `kStopTime`, `kStopIter`, `kCutOff`, and `kError - 1` | Continue branching when the node has a finite bound and the relaxation has a solution. |
| `kLowPrecision` | Preserve valid bounds and incumbent, without authorizing branching solely from that status. |
| `kError`, `kBlockLocked` | Propagate the error without branching or publishing a new incumbent from that result. |

The cases run in both optimization senses, under eager and lazy bounding,
with and without tree retention. Other regressions check that:

- an early stop without a finite bound retains the available incumbent and
  leaves the node unresolved, rather than claiming a finite certificate;
- a child without a relaxation solution cannot reuse the branching information
  of its parent, and the solve can resume when a solution becomes available;
- bound-only results preserve their bounds and the first unresolved status,
  without a priority among statuses; a cutoff without usable branching
  information or a pruning certificate becomes `kLowPrecision`;
- mixed usable and bound-only relaxations allow the search to reach the
  optimum in either solver order, with serial or parallel evaluation and
  eager or lazy bounding; the cases include raw recoverable statuses,
  `kError - 1`, and results normalized to `kOK` by the wrapper;
- raw `kLowPrecision` remains unresolved in the mixed-result cases, without
  losing its finite values or overriding the first recorded status;
- sufficiently advanced inner results that the concrete relaxation reports
  as `kOK` remain usable, including on an unchanged repeat solve.

`test_cutoff_without_incumbent` checks that a cutoff with neither an incumbent nor
usable branching information leaves the node unresolved with its bound, under
both bounding protocols and with serial or parallel relaxation evaluation.

`test_inner_time_budget` checks that each inner solve receives the smaller of
its own time limit and the time left to the outer solve, and that its original
limit is restored afterwards. With no time left, the inner solver is not run
and its previous bounds and solutions are not read; a later solve can resume.

### Certificates, logging and cleanup

Invalid certificates (`NaN`, an infinite bound contradicting a finite true
solution, or inconsistent primal and relaxation bounds) must return `kError`
without publishing an incumbent. Separate cases exercise large finite
objective values and closure within numerical tolerance, retaining the
computed bounds on repeated solves.

The CSV checks require one row per evaluated serial node, including an eager
root, infeasibility and failed evaluations. Lazy nodes pruned without evaluation
must not add evaluation rows. Root and child failures also exercise the cleanup
of nodes and their branching and undo Changes.

`test_compute_unlock` probes the recursive Solver mutex from a separate thread
after a successful solve, an early `kError`, and an exception during worker
configuration. It also checks the `kError` state after the exception and a
successful subsequent solve.

## Additional configurations

The additional configurations, read from this directory, are:

| Configuration | Regression |
| --- | --- |
| [InnerBSCfg-bad-certificate.txt](InnerBSCfg-bad-certificate.txt) | Invalid certificates must not publish an incumbent. |
| [InnerBSCfg-log-failure.txt](InnerBSCfg-log-failure.txt) | Root and child failures must be logged without leaking the failed child. |
| [InnerBSCfg-mixed-partial.txt](InnerBSCfg-mixed-partial.txt) | Usable and bound-only results are combined in either solver order, serially and in parallel. |
| [InnerBSCfg-near-bound.txt](InnerBSCfg-near-bound.txt) | Closing within tolerance must retain the computed bound. |
| [InnerBSCfg-parallel-partial.txt](InnerBSCfg-parallel-partial.txt) | Bound-only results retain their bounds and the first unresolved status. |
| [InnerBSCfg-partial.txt](InnerBSCfg-partial.txt) | Return-code boundaries, missing bounds or solutions, sufficiently advanced results, and resumption. |
| [InnerBSCfg-time-probe.txt](InnerBSCfg-time-probe.txt) | Propagated time limits, restoration of the original limit, and no stale result extraction after an expired budget. |

## Parallelism covered

The mixed-result tests exercise parallel evaluation of multiple relaxations
at a node (`intThreadForDifferentSolvers`). The bound providers' `compute()`
methods are read-only, alongside a single solver that changes the `Block`.
During serial result extraction, each provider writes its own feasible
zero-cost integer solution before capturing it, so a fractional point left
by another solver cannot be published as that incumbent.

The fixture does not supply independent Block copies. `test_parallel_no_copy`
checks that parallel depth-first exploration (`intMaxThread` > 1) is refused
with an exception in that case. Full parallel worker exploration on copyable
Blocks is not exercised here. The lock-release regression also checks a
setup exception before exploration.


## Authors

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.

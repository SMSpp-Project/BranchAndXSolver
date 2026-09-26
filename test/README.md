# test

A tester for `BranchAndXSolver` that needs nothing but the core SMS++ library.

It builds in memory an `AbstractBlock` with a few `ColVariable`, a
`BoxConstraint` on each of them and a separable quadratic `DQuadFunction` as
`Objective`, to be minimized (convex) or maximized (concave) over the integer
points of the box. The integer optimum is known in closed form (in each
coordinate, the integer of the box nearest to the unconstrained optimum), and
the closed form is checked against the enumeration of all the integer points.

The core has no `RelaxationSolver`, so the tester defines the inner `Solver`
that `BranchAndXSolver` drives, both `BoxSolver` underneath: `BoxRelaxation`
solves the continuous relaxation (the `ColVariable` are continuous in the
`Block`, their integrality is the knowledge of this `Solver`) and branches on
its most fractional coordinate by tightening the `BoxConstraint`, while
`BoxRounding` is a heuristic `ChangeSolver` rounding the solution of the
relaxation inside the box. They are given to `BranchAndXSolver` by the
`BlockSolverConfig` in `InnerBSCfg.txt` (the relaxation alone) and
`InnerBSCfg-rounding.txt` (the relaxation and the heuristic), which the tester
reads from its working directory, this one.

The cases are:

- every exploration strategy (`intSolveMethod`) under both bounding protocols
  (`intBoundingProtocol`), minimizing, and every strategy maximizing, each run
  having to find the optimum, with both bounds equal to it and the optimal
  integer solution written back in the `Block`;
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
  nothing changed.

The relaxation applies the branching to the `Block` rather than to its own
state, so the parallel depth-first exploration (`intMaxThread` > 1), whose
workers move their own inner `Solver` independently, is not exercised.

The exit code is 0 when every check passes, printing `All tests passed!!`, and
1 otherwise. The `makefile` builds the executable including the
`BranchAndXSolver` module and the core SMS++ library.


## Authors

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.

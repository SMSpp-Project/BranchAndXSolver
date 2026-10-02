# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- a unit test of the module that needs nothing but the core: an integer box
  problem with a separable quadratic objective, whose optimum is known in
  closed form, solved with a continuous relaxation and a rounding heuristic
  built on BoxSolver, under every exploration strategy and bounding protocol,
  minimizing and maximizing, with infeasible boxes, an integer root, node and
  time budgets, the parameters, a detach and reattach and a change of the
  objective; `ctest -L BranchAndXSolver` runs it, and the CI of the module
  runs only it

### Changed

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

### Fixed

- a node whose evaluation does not end with `kOK` (e.g., a relaxation stopped
  by its time limit) no longer leaves the Changes of its path applied: the
  solvers go back to the root before the tree is discarded, so that the best
  solution is then written into a `Block` without the fixings of that node

- an exploration stopped by the node budget reported kOK, i.e., optimality,
  and kInfeasible when no incumbent had been found yet: the status is now
  kStopIter or kStopTime whenever open nodes that can improve the incumbent
  are left, and kOK otherwise, whichever budget is left

- dblMaxTime, dblRelAcc and dblAbsAcc were never stored, the base Solver
  classes not doing it, so the time budget and the tolerances of the pruning
  were always the default ones

- a compute() after a stop, nothing having changed, returned the old status
  instead of going on

- a Solver detached from a Block could not be attached again (the
  GlobalInformation was declared twice), and one attached to another Block
  kept the inner Solver of the old one: the inner Solver, the retained tree
  and the best solution are now dropped when the Block changes

- the name of the BlockSolverConfig file given as a temporary was ignored,
  the string parameter being taken by the wrong overload; the factory reset
  (its empty default) threw, and a second BlockSolverConfig added its inner
  Solver to those of the first instead of replacing them

- the incumbent found by a heuristic ChangeSolver was saved as an empty
  Solution, which could not be written back

- the bounds of a proven infeasible problem were not both the infinity of
  the sense, and invalid intSolveMethod / intBoundingProtocol were accepted,
  to throw inside compute() with the Solver locked

- the branching Changes of a node evaluated again, as in a reoptimization,
  were leaked

- on macOS a program linking the module lost the classes the module
  registers in the factories when the linker dropped the library, as it
  does under `-dead_strip_dylibs`, which conda sets: the target now asks the
  linker for the symbol that forces the module in (`-u`), which ld64,
  unlike the ELF linker, counts as a use of the library

## [0.1.0] - 2026-09-12

### Added

- the working Branch-and-Bound core merged from the original
  prototype repository: DFS / BFS / BestFS exploration driven by the ChangeSolver /
  RelaxationSolver concepts, inner Solver provided via BlockSolverConfig,
  private to the BranchAndXSolver (Modification forwarded to them)

- tolerance-based pruning on the inherited dblRelAcc / dblAbsAcc; node and
  time budgets on the inherited intMaxIter / dblMaxTime

- SMS++-style CMake build system and module makefiles

### Changed

- the ChangeSolver / RelaxationSolver concepts and GroupChange moved to the
  SMS++ core (BnXSolver branch), where they belong

- the R3-Block-per-node prototype of the original BranchAndXSolver is
  superseded

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

[Unreleased]: https://gitlab.com/smspp/branchandxsolver/-/compare/0.1.0...develop
[0.1.0]: https://gitlab.com/smspp/branchandxsolver/-/tags/0.1.0


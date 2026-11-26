# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- the value `GlobalInformation::str_BestBound` that the relaxations read is
  the global relaxation bound of the exploration under way, the best among
  the bounds of the open nodes, of the node being evaluated and of those set
  aside, and the incumbent, which every open list keeps in constant time
  (but the best-estimate one, in logarithmic time) and the parallel
  depth-first search takes from its pool and from the subtrees claimed by
  the workers, and the search of PBeFS gives each time it improves, that of
  the root standing until it gives one

- the parameters `intPBeFSSpawnDepth` and `intPBeFSBudget`, the spawn depth
  of the schemes `depthbounded` and `ordered` of PBeFS and the backtracks
  after which the scheme `budget` spawns tasks, whose defaults (those of
  PBeFS) leave most of the workers idle on more than 4 of them

- the parameter `intPBeFSSplitShared`, by which the workers of the scheme
  `selfsplit` of PBeFS share their incumbent

- `dblUpCutOff` (`dblLwCutOff` when maximizing) is the incumbent that the
  exploration starts from, so that the nodes are pruned against it, and a
  solve with no solution better than it ends in `kCutOff`

- the race of PBeFS before its scheme, by `vintPBeFSRace` (the strategy of
  each racer) and `intPBeFSRaceNodes`, and
  `vstrWorkerBlockSolverConfigurationFile`, the BlockSolverConfig of each
  worker of the parallel explorations but the first, so that the racers
  (and the workers) may have differently configured inner Solver

- the cutoff probing, by `intCutOffProbing`, `intCutOffProbingRuns`,
  `dblCutOffProbingGap` and `dblCutOffProbingFactor`: the solve is a
  sequence of explorations, each against a target incumbent close to the
  bound of the root, which the relaxations are given as their cutoff, and
  moved away from it after each exploration that proves no solution better
  than it, until one finds a solution, whose optimality the bound proven by
  the last refuted target makes quicker to prove

### Changed

- the makefiles give the shared library of the module by default, a DLL on
  Windows, as CMake does, and its objects with `SMSppSTATIC` [see
  `lib/makefile-shared` of SMS++]

- the library is kept in whoever links it by a source that calls its
  `SMSpp_force_load_` function, the tests and the tools take the static
  libraries whole by the `WHOLE_ARCHIVE` feature of CMake and copy the DLLs
  next to them only where there are any, and the installed run-time path is
  the one of `smspp_install_rpath()`, with no test on the system or on the
  compiler [see `SMSppModule.cmake` of SMS++]; CMake 3.24 or later is needed

- inner relaxation statuses in `[kOK, kError)` can drive branching, including
  time, iteration and cutoff stops and custom recoverable statuses; an early
  stop also requires a finite node bound and an available relaxation solution

- `kLowPrecision`, which lies outside that interval, contributes valid bounds
  and incumbents without authorizing branching; unpruned nodes without usable
  branching information remain unresolved, and their first recorded status
  is reported without assigning priorities (an unresolved cutoff becomes
  `kLowPrecision`)

### Fixed

- an inner Solver that throws during the exploration leaves none of the
  Changes of the nodes on the Block: the Solver are taken back to the root
  and the tree is deleted before the exception goes on, where the Block was
  left at the node, so that the Solver computed after it saw a restricted
  problem and gave a bound beyond the optimum; in the parallel explorations
  the exception of a worker stops the others and is thrown by `compute()`
  once they are back, where it ended the program

- reject `NaN` relaxation bounds, non-finite true solution values and
  contradictory certificates before updating the incumbent; construct a
  replacement `Solution` before deleting the previous one, using the same
  validation in serial and parallel evaluation

- preserve computed bounds when pruning within tolerance and when repeating
  an unchanged completed solve, and retain the bounds and feasible solutions
  supplied by interrupted relaxations

- a cutoff prunes a node only when its bound proves that it cannot improve
  the incumbent, rather than treating every `kCutOff` as a pruning certificate;
  an unpruned node without usable branching information keeps its bound and
  is set aside, including when there is no incumbent

- release the Solver with RAII on normal returns, early errors and exceptions;
  failed explorations leave its state at `kError` while propagating the
  exception

- restore the parent state before a recursive worker returns an error or a
  deeper resource-limit stop, and delete the branching Changes left unexplored

- write one CSV row per evaluated serial node, including an eager root and
  failed evaluations, without counting lazy nodes pruned before evaluation

- delete an eagerly evaluated child that fails before it is linked into the
  tree, releasing its undo Change without double-deleting its branching Change

- the workers of the parallel exploration no longer share the Block: their
  inner Solver were all attached to it, and a relaxation that applies the
  branching to the Block (say, `LagrangianDualRelaxationSolver` fixing the
  variables of its sub-Blocks) made each worker see the fixings of the
  others; the first worker now drives the inner Solver of the Block, each of
  the others an identical copy of the Block (`get_R3_Block()`) with its own
  inner Solver, a Block that gives no copy is refused with an exception, and
  a worker that stops leaves its Solver at the root of its subtree

- `dblMaxTime` bounds the whole solve: it was only checked between the nodes,
  and a long `compute()` of an inner Solver ran past it; each one now has the
  `dblMaxTime` of its own or the time left, whichever is smaller, a node that
  it leaves unsolved is set aside with the bound the Solver gives, and the
  solve ends with `kStopTime`

## [0.4.0] - 2026-11-16

### Changed

- `BranchAndXSolver::concerned_by()` is what any of its inner Solver
  reads (all the kinds before they are created), and `add_Modification()`
  passes a Modification only to the inner Solver that read its kind, as the
  Block does with the Solver registered to it

### Fixed

- the installed module `smspp.branchandxsolver` could not be imported, since
  its package gave the interface of the module the target of the library by
  its name in the build tree, which an installed package does not have

- `BranchAndXSolver.h` compiles without warnings under `-Wall`: the
  constructor initializes the members in the order they are declared

- `reseedFrontier()` marks as evaluated the nodes of the frontier it
  evaluates, which `explore()` evaluated a second time when extracting them

- the lazy bounding protocol fences, without evaluating it, a node extracted
  with a bound that can no longer improve the incumbent, as its
  documentation says, rather than evaluating it first

- `explore()` compiles without the `-Wenum-compare` warning of its return

## [0.3.0] - 2026-11-03

### Added

- the makefiles build the C++20 module `smspp.branchandxsolver` when
  `SMSppMOD` is not empty (e.g., `make SMSppMOD=1`), with clang++ 16 or
  later [see `lib/makefile-mod` of SMS++]

- the best-estimate search (`intSolveMethod` `BestEstimate`): the open node
  of best estimate of the value under it is extracted first, the estimate
  being given by the RelaxationSolver that branched, if it is also a
  `BranchEstimator` (`BranchingRule.h` of SMS++), and the bound of the
  parent otherwise, and never better than the bound of the node, which
  still prunes it

- the parameter `intTreeExplorer`, by which the tree is explored by the
  branch-and-bound of PBeFS instead of the explorations of the
  BranchAndXSolver, on `intMaxThread` workers and by the scheme of PBeFS
  that `strPBeFSScheme` names (`"shared"` by default), which is given to
  PBeFS as it is, with the same strategies (`BestEstimate` included),
  bounding protocols, tolerances and budgets; each worker moves its own
  inner Solver along the Changes of the path of a node and evaluates it
  against the incumbent that PBeFS gives it, so that the outcome is that of
  the scheme, also when the scheme makes it the same at every run; PBeFS is
  found with `find_package` or downloaded at a fixed commit, and in the
  makefiles through `extlib/makefile-libPBeFS`

- the C++20 module `smspp.branchandxsolver` (`src/BranchAndXSolver.cppm`),
  which exports every name that the headers declare at namespace scope and
  re-exports `smspp`, so that a source file can write `import
  smspp.branchandxsolver;` in place of the `#include` of the headers; it is
  built with the CMake option `SMSpp_BUILD_MODULE`, as the target
  `SMS++::BranchAndXSolver_module`

- the documentation of the class gives the formulation of the search: the
  tree and the algebra of the undo Changes, the move along the common
  ancestor, the pruning test with its tolerances, the bound and the
  incumbent of the Solver of a node, the exploration strategies, the two
  kinds of parallelism and the tightening, with the references

### Changed

- the README leaves the channels that distribute SMS++ to the
  installation guide of the wiki (Installing-SMS++), keeping the names
  of the packages of the module

- a node whose relaxation returns `kStopTime`, `kStopIter` or
  `kLowPrecision` (e.g., a Lagrangian Dual not solved to the required
  accuracy whose solution gives nothing to branch upon) no longer stops the
  whole exploration: it is set aside, neither branched upon nor fenced, the
  other nodes are explored, and its relaxation bound enters the bound of the
  solve as that of an open node; `compute()` returns `kLowPrecision` if some
  node so set aside can still improve the incumbent, while a relaxation
  returning `kError` still stops the solve

- the pipeline builds the module also with its makefiles, after CMake, and
  runs on a change of the makefiles too [see `ci/makefile` of the umbrella]

- `SMSpp_BUILD_MODULE` no longer stops the configuration when the C++20
  module cannot be built, i.e., with a CMake older than 3.28, a generator
  other than Ninja or Visual Studio, a compiler other than Clang 16 (with
  `clang-scan-deps`) or MSVC 19.34, or a module it imports that is not
  there: the option is turned off with a warning that says why, and the
  build goes on with the headers [see `cmake/CheckModule.cmake` of SMS++]

- the exploration allocates less, with the same nodes in the same order: a
  move between two nodes collects the path in a buffer of the thread rather
  than in a list built each time, the children of a node are a vector, the
  surviving children of the nodes share one buffer, the nodes of a
  depth-first worker live in its stack frames, the best-first open lists
  compare with a function object rather than through a std::function, and
  an idle worker of the parallel depth-first search waits on a condition
  variable rather than polling the pool every 100 microseconds

- the outstanding Modification are read through `Modification::changes()`
  rather than classified by the RelaxationSolver, each one by itself rather
  than by the bits of all of them or-ed together: the retained tree survives
  those that relax or restrict the problem without touching what its
  branching Changes act upon (added Variable or Constraint, Variable
  becoming integer, changed sides or Objective), and it is discarded when
  the data of some Variable change, when Variable or Constraint are
  removed, when some Variable becomes continuous or when a physical
  Modification does not say what it changes; the infeasible part of its
  frontier stays fenced whenever the feasible region cannot grow (not only
  when the Objective alone changes)

### Fixed

- on Windows, the executables that copy their DLLs next to themselves
  are linked one at a time (a job pool of Ninja), since two copies of the
  same DLL into one directory collided (e.g., with the clang++ of MSYS2)

- with the compilers of MinGW (MSYS2) the target asks whoever links it
  for the symbol that forces the module in, as with MSVC and on macOS:
  the linker imports a DLL only for the symbols asked for, so that a
  program that used none had the module left out and its classes
  missing from the factory

- the makefiles of the tests and of the tools build also with g++: they
  take the option that stops at the first error and the one of the
  debugging information, which clang++ and g++ spell differently
  (`-ferror-limit=1` and `-glldb` against `-fmax-errors=1` and `-ggdb`),
  from `SMSppSWERR` and `SMSppSWDBG` of `lib/makefile-OS` of SMS++

- the root is fenced with `dblRelAcc` and `dblAbsAcc`, as any other node,
  rather than with no tolerance at all: a root whose relaxation bound was
  the value of the true solution it gave up to the accuracy of the
  relaxation (e.g., a Lagrangian dual solved to 1e-5 whose convexified
  solution is integral, on the ORLib instances of the single-source CFL)
  was branched, and the RelaxationSolver had nothing left to branch upon

## [0.2.0] - 2026-10-09

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

- the pruning tests and the checks of an improved incumbent read the shared
  incumbent (the atomic of the GlobalInformation) rather than the plain
  double kept by the Solver, so that the workers of the parallel
  exploration never read a value being written

- the incumbent is given to the relaxations of each node as their cutoff
  (`dblUpCutOff` when minimizing, `dblLwCutOff` when maximizing), so that a
  relaxation proven beyond it stops there and the node is fenced as
  infeasible; a root fenced this way, with an incumbent already there, proves
  that incumbent optimal instead of the problem infeasible; a relaxation that
  returns `kCutOff` fences the node by bound, with the bound it reports

- the bound a node gets from the RelaxationSolver is called the relaxation
  bound, which says where it comes from and nothing about how it is
  computed: `Node::get_dual_bound()` / `set_dual_bound()` are now
  `get_bound()` / `set_bound()`; the column of the per-node log keeps its
  name `dualBound`

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

### Fixed

- in the parallel depth-first exploration a worker that hands a child over
  to the pool asks it for its branches before moving its Solver back to the
  father, as everywhere else: the branches came out of the relaxation of the
  father, and a Solver that keeps a state across the Changes (say, the
  incremental greedy of BinaryKnapsackBlock) gave a child a wrong state, and
  a value beyond the optimum

- an exploration stopped by a budget (`intMaxNodes`, `dblMaxTime`) returns as
  the bound on the side of the open nodes (`get_lb()` when minimizing,
  `get_ub()` when maximizing) the best among their relaxation bounds and the
  incumbent, rather than an infinite one, also when it is a node whose
  relaxation is not solved (e.g., stopped by its own budget) that stops it,
  the bound of that node counting among them; the parallel depth-first search
  still claims none

- a child starts from the relaxation bound of its parent, which is valid for
  it, and its evaluation can only improve it: with the eager bounding protocol
  and with the parallel depth-first search it started from an infinite one, so
  that a child whose relaxation stopped early carried a bound weaker than that
  of its parent

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

[Unreleased]: https://gitlab.com/smspp/branchandxsolver/-/compare/0.4.0...develop
[0.4.0]: https://gitlab.com/smspp/branchandxsolver/-/compare/0.3.0...0.4.0
[0.3.0]: https://gitlab.com/smspp/branchandxsolver/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/branchandxsolver/-/compare/0.1.0...0.2.0
[0.1.0]: https://gitlab.com/smspp/branchandxsolver/-/tags/0.1.0

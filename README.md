# BranchAndXSolver

`BranchAndXSolver` implements the Solver interface for a Relaxation-Agnostic
Branch-and-X (RABaX) Solver within the SMS++ framework: a generic enumerative
solver where "X" stands for Bound / Cut / Price (cuts and pricing still to
come). It is built on top of two abstract Solver concepts that live in the SMS++
core (`ChangeSolver.h`, currently on the `BnXSolver` core branch):

- `ChangeSolver`, a :Solver that can apply() a `Change` [see Change.h in the
  SMS++ core] to the Block it is attached to - and the undo Change that
  apply() returns - so that the same Solver object can be efficiently moved
  between the nodes of an enumeration tree;

- `RelaxationSolver`, a :ChangeSolver solving a *relaxation* of the problem:
  besides the relaxation value (a valid dual bound) it can produce *true*
  solutions of the original problem (valid primal bounds) and, foremost, it
  can branch(), i.e., produce the Changes generating the children of the
  current node.

The module provides:

- `BranchAndXSolver`, the enumerative :Solver itself: depth-first,
  breadth-first or best-first exploration (`intSolveMethod`), node / time
  limits (the inherited `intMaxIter` / `dblMaxTime`), tolerance-based pruning
  (the inherited `dblRelAcc` / `dblAbsAcc`), inner Solver provided either directly or through a
  BlockSolverConfig (`strNameOfBlockSolverConfigurationFile`);

- `GroupChange`, a Change grouping several Changes applied as one (not used
  yet);

- `ParallelSolver`, a sketch of the master / clones machinery for the future
  parallel exploration (not used yet).

**WARNING: WORK IN PROGRESS.** The Branch-and-Bound core is functional and
validated (see the batches of the SMS++ tests project), but cuts, pricing,
the reoptimization machinery and the parallel exploration are still to come.

## Getting started

These instructions will let you build the `BranchAndXSolver` module on
your system.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.


### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).
When built from the SMS++ umbrella project, enable it with
`-DBUILD_BranchAndXSolver=ON`.

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(BranchAndXSolver)
target_link_libraries(<my_target> SMS++::BranchAndXSolver)
```

### Build and install with makefiles

The `makefile` exports the usual SMS++ module macros (`$(BAXSLVOBJ)`,
`$(BAXSLVH)`, `$(BAXSLVINC)`) given `$(BAXSLVSDR)`, the SMS++ core macros and the
core SMS++ ones (`$(SMS++OBJ)`, `$(SMS++H)`, `$(SMS++INC)`).

## Getting help

If you need support, you want to submit bugs or propose a new feature, you
can [open a new issue](https://gitlab.com/smspp/BranchAndXSolver/-/issues/new).

## Contributing

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Filippo Magi**  
  Dipartimento di Informatica  
  Università di Pisa

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa

## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html).

## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.

# BranchAndXSolver

`BranchAndXSolver` implements the Solver interface for a Relaxation-Agnostic
Branch-and-X (RABaX) Solver within the SMS++ framework: a generic enumerative
solver where "X" stands for Bound / Cut / Price. It is built on top of two
abstract Solver concepts that live in the SMS++ core (`Change.h` /
`ChangeSolver.h`):

- `ChangeSolver`, a :Solver that can apply() a `Change` to the Block it is
  attached to, and the undo Change that apply() returns, so that the same
  Solver object can be efficiently moved between the nodes of an enumeration
  tree;

- `RelaxationSolver`, a :ChangeSolver solving a *relaxation* of the problem:
  besides the relaxation value (a valid dual bound) it can produce *true*
  solutions of the original problem (valid primal bounds), branch() to
  generate the children of the current node, and separate() valid cuts.

The module provides `BranchAndXSolver`, the enumerative :Solver itself:
depth-first, breadth-first or best-first exploration, cut separation,
reoptimization of the enumeration tree, parallel exploration, and an inner
Solver provided either directly or through a BlockSolverConfig. See the class
documentation for the full list of algorithmic parameters.

`GroupChange`, the Change grouping the per-node branching Changes, lives in
the SMS++ core together with the `ChangeSolver` / `RelaxationSolver` concepts.

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

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `BranchAndXSolver` has to include a "main makefile" of
the module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or [makefile-s](makefile-s)
including all necessary libraries but not the "core SMS++" one (for the common
case in which this is used together with other modules that already include
them). The makefiles in turn recursively include all the required other
makefiles, hence one should only need to edit the "main makefile" for
compilation type (C++ compiler and its options) and it all should be good to
go. In case some of the external libraries are not at their default location,
it should only be necessary to create the `../extlib/makefile-paths` out of the
`extlib/makefile-default-paths-*` for your OS `*` and edit the relevant bits
(commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.

## Getting help

If you need support, you want to submit bugs or propose a new feature, you
can [open a new issue](https://gitlab.com/smspp/BranchAndXSolver/-/issues/new).

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.

## Authors

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

### Contributors

## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.

## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.

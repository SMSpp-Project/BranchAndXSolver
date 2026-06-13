# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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


# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- initial module skeleton generated from ModuleTemplate

- `CSPBlock`, the :Block of the constraint satisfaction problems of
  constraint programming (CP) and of their weighted version: its physical
  representation (universes and current domains of the variables, a
  catalogue of constraints with weights and activity, costs of the values,
  groups of the variables), the check of an assignment and its objective
  value, the reader and the writer of a subset of FlatZinc, netCDF,
  `CSPBlockMod` with `CSPBlockRngdMod` and `CSPBlockSbstMod`, and
  `CSPBlockChange` for the domains, with its undo

- the unit tests and a battery on the n queens, a map colouring, sudoku and
  the scheduling of tasks on a resource, whose solutions are enumerated and
  counted

- the abstract representation of `CSPBlock`: an integer ColVariable per
  variable, a box per variable on its current domain, an FRowConstraint per
  linear constraint, and the global constraints as `CSPConstraint`, which
  check an assignment through the physical representation and give a linear
  relaxation where a cheap one is valid

- `FlatZincCSPSolver`, the Solver that solves a CSPBlock by writing its
  FlatZinc and running an external FlatZinc solver on it as a process, its
  executable being the parameter `strFlatZincSolver` and its arguments
  `vstrFlatZincArgs`, and reading the solution back; a Modification of the
  CSPBlock makes the next `compute()` write and solve the model again; the
  battery checks it against the enumeration of the solutions, and
  `tools/install_flatzinc_backends.sh` installs the MiniZinc bundle

### Changed

- the makefiles give the shared library of the module by default on Linux
  and on macOS, as CMake does, and its objects with `SMSppSTATIC` [see
  `lib/makefile-shared` of SMS++]

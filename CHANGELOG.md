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

- `FlatZincCSPSolver`, the Solver that solves a CSPBlock by writing its
  FlatZinc and running an external FlatZinc solver on it as a process
  (OR-Tools CP-SAT, Gecode, Chuffed, ...), its executable being the parameter
  `strFlatZincSolver` and its arguments `vstrFlatZincArgs`, and reading the
  solution back; a default configuration and a test on the 4 queens, skipped
  when no backend is installed

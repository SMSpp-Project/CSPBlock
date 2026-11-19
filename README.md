# CSPBlock

[![Pipeline](https://gitlab.com/smspp/cspblock/badges/develop/pipeline.svg)](https://gitlab.com/smspp/cspblock/-/pipelines?ref=develop)
[![Release](https://gitlab.com/smspp/cspblock/-/badges/release.svg)](https://gitlab.com/smspp/cspblock/-/releases)
[![Documentation](https://img.shields.io/badge/docs-Doxygen-blue)](https://smspp.gitlab.io/smspp-project/)
[![C++ Standard](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

`CSPBlock` is a SMS++ :Block for the constraint satisfaction problems of
constraint programming (CP), i.e., finding values of a set of variables, each
with a finite set of integer values, that satisfy a set of constraints of a
catalogue (positive and negative tables, linear equalities and inequalities,
alldifferent, element and cumulative), and for their weighted version, where
the soft constraints may be violated at the cost of their weight and each
value of each variable may have a cost. Its physical representation is the
universe of each variable, the sorted values it may ever take, and its
current domain within it, so that a restriction of a domain (e.g., by a
branching) keeps the indices of all the values; the constraints, each with
its kind, its scope and its integer parameters, stored one after the other;
their weights and whether each is active; the costs of the values; and the
groups the variables may be dealt out to. The `CSPBlock` has no abstract
representation: its solutions are vectors of values, one per variable, which
it checks against the constraints and whose objective value it computes.

A `CSPBlock` is read from a subset of FlatZinc, the language into which
MiniZinc flattens a model for a solver (the documentation of `load()` lists
the constraints and the annotations it reads), and written back into it; it
is serialized into and deserialized out of netCDF like any other :Block. Its
domains, universes, costs, weights and constraints can be changed or
extended, with the `CSPBlockMod` that tell the Solver what changed, and the
`CSPBlockChange` restricts the domains as a branching does (assigning,
removing or splitting a value) and brings them back.

A `CSPBlock` is solved by `FlatZincCSPSolver`, which writes its FlatZinc and
runs an external FlatZinc solver on it as a separate process (OR-Tools CP-SAT
by default, or Gecode, Chuffed and others), chosen by the `strFlatZincSolver`
and `vstrFlatZincArgs` parameters, and reads the solution back.


## Getting started

These instructions will let you build the `CSPBlock` module on
your system.

The channels that distribute SMS++ are listed in its [installation
guide](https://gitlab.com/smspp/smspp-project/-/wikis/Installing-SMS++); what
follows is about building the module yourself.

### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.

- To solve a `CSPBlock` with `FlatZincCSPSolver`, a FlatZinc backend on the
  `PATH`: `tools/install_flatzinc_backends.sh` installs the MiniZinc bundle
  (the `minizinc` driver and the backends OR-Tools CP-SAT, Gecode and Chuffed)
  under a prefix and with no root, which fits a cluster account such as
  fenchel; on macOS, `brew install minizinc`. The `strFlatZincSolver` and
  `vstrFlatZincArgs` parameters select the driver and the backend.

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

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(CSPBlock)
target_link_libraries(<my_target> SMS++::CSPBlock)
```

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `CSPBlock` has to include a "main makefile" of
the module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). These in turn recursively include all the
required other makefiles, hence one should only need to edit the "main
makefile" for compilation type (C++ compiler and its options) and it all
should be good to go. In case some of the external libraries are not at their
default location, it should only be necessary to create the
`../extlib/makefile-paths` out of the `extlib/makefile-default-paths-*` for
your OS `*` and edit the relevant bits (commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-external-libraries)
for further details.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you
can [open a new issue](https://gitlab.com/smspp/cspblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

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

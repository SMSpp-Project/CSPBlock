#!/bin/bash
# --------------------------------------------------------------------------- #
#    Install the FlatZinc backends of FlatZincCSPSolver                       #
#                                                                             #
#    FlatZincCSPSolver runs an external FlatZinc solver on the FlatZinc of a  #
#    CSPBlock [see include/FlatZincCSPSolver.h]. This script installs the     #
#    official MiniZinc bundle, which carries the `minizinc` driver and the    #
#    backends OR-Tools CP-SAT, Gecode and Chuffed, under a prefix and with no #
#    root; it is idempotent, doing nothing if the bundle of the asked version #
#    is there.                                                                #
#                                                                             #
#        ./install_flatzinc_backends.sh [--prefix=<dir>] [--version=<ver>]    #
#                                                                             #
#    It then prints the line that puts `minizinc` on the PATH and how the     #
#    backends map to the strFlatZincSolver and vstrFlatZincArgs parameters.   #
#    On macOS use `brew install minizinc`, which bundles the same backends.   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
set -euo pipefail

version=2.8.5
prefix="$HOME/.local/minizinc"

for arg in "$@"; do
  case "$arg" in
    --prefix=*) prefix="${arg#*=}" ;;
    --version=*) version="${arg#*=}" ;;
    *) echo "unknown argument: $arg" >&2; exit 1 ;;
  esac
done

bundle="MiniZincIDE-${version}-bundle-linux-x86_64"
base="https://github.com/MiniZinc/MiniZincIDE/releases/download"
url="${base}/${version}/${bundle}.tgz"

if "${prefix}/bin/minizinc" --version 2>/dev/null | grep -q "${version}"; then
  echo "MiniZinc ${version} already installed in ${prefix}."
else
  echo "Installing the MiniZinc ${version} bundle in ${prefix}..."
  tmp="$(mktemp -d)"
  trap 'rm -rf "$tmp"' EXIT
  curl -fsSL "$url" -o "${tmp}/bundle.tgz"
  # the prefix is emptied only if it holds an older bundle
  if [ -n "$(ls -A "$prefix" 2>/dev/null)" ]; then
    if [ ! -x "${prefix}/bin/minizinc" ]; then
      echo "${prefix} is not empty and holds no MiniZinc bundle" >&2
      exit 1
    fi
    rm -rf "$prefix"
  fi
  mkdir -p "$prefix"
  tar -xzf "${tmp}/bundle.tgz" -C "$prefix" --strip-components=1
fi

echo
echo "Put the driver on the PATH with:"
echo "    export PATH=\"${prefix}/bin:\$PATH\""
echo
echo "The backends, for the FlatZincCSPSolver parameters, are:"
echo "    strFlatZincSolver minizinc"
echo "    vstrFlatZincArgs  2 --solver com.google.ortools.sat"
echo "    (or org.gecode.gecode, or org.chuffed.chuffed)"
# --------------------------------------------------------------------------- #

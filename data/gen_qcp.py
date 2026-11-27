#!/usr/bin/env python3
"""
Writes in FlatZinc some quasigroup completion problems of the MiniZinc
benchmarks: the first 5 instances of order 10 (67 holes) and of order 15
(120 holes) of their directory QCP, each a partial latin square to be
completed, a satisfaction problem.

An instance is a variable v_<k> per cell k, row by row, its value in
0 .. n - 1 or the one value of a cell that is given, and a
fzn_all_different_int per row and per column, as in the model of the
benchmarks; the solve item is satisfy.

usage: gen_qcp.py [ <output directory> ]

The output directory is fzn/qcp next to this script if none is given; it
gets the instances, optima.txt (the instance and 0, the value of any of
its solutions, the instances being satisfiable) and the README of the
family. The models are downloaded from the repository of the benchmarks,
at a fixed commit.
"""

import os
import re
import sys
import urllib.request

COMMIT = "26bcd0a78433025f7b6896a6fa8caf128795760b"
SITE = ("https://raw.githubusercontent.com/MiniZinc/minizinc-benchmarks/"
        + COMMIT + "/QCP/")

NAMES = (["qcp-10-67-{}".format(i) for i in range(5)] +
         ["qcp-15-120-{}".format(i) for i in range(5)])

README = """\
data/fzn/qcp: quasigroup completion in FlatZinc
===============================================

10 quasigroup completion problems of the MiniZinc benchmarks
(https://github.com/MiniZinc/minizinc-benchmarks, directory QCP, MIT
license): the first 5 of order 10, with 67 holes, and of order 15, with 120
holes. Each is a partial latin square to be completed, i.e., one variable
per cell, a given cell having one value, and one fzn_all_different_int per
row and per column; the solve item is satisfy.

All the instances are satisfiable, and optima.txt gives for each the value
0 of any of its solutions, which a Solver that proves the instance
unsatisfiable, or that finds no solution, is checked against.

The instances are written by gen_qcp.py in data/. The problem is described
in

  C.P. Gomes, B. Selman, Problems with structure: the quasigroup completion
  problem, Proceedings of AAAI-97, 221-226, 1997.
"""


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode("latin-1")


def write_fzn(name, text):
    out = ["% {} of the MiniZinc benchmarks: complete the latin square"
           .format(name)]
    for lo, hi, v in re.findall(r"^var\s+(-?\d+)\s*\.\.\s*(-?\d+)\s*:\s*"
                                r"(\w+)\s*;", text, re.M):
        out.append("var {}..{}: {};".format(lo, hi, v))
    for args in re.findall(r"constraint\s+all_different\s*\(\s*\[([^\]]*)\]",
                           text):
        out.append("constraint fzn_all_different_int([{}]);".format(
            ",".join(a.strip() for a in args.split(","))))
    out.append("solve satisfy;")
    return "\n".join(out) + "\n"


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    odir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "fzn",
                                                               "qcp")
    os.makedirs(odir, exist_ok=True)

    rows = []
    for name in NAMES:
        with open(os.path.join(odir, name + ".fzn"), "w") as f:
            f.write(write_fzn(name, fetch(SITE + name + "_ext.mzn")))
        rows.append("{}.fzn 0".format(name))

    with open(os.path.join(odir, "optima.txt"), "w") as f:
        f.write("# the instance and the value of any of its solutions\n")
        f.write("\n".join(rows) + "\n")
    with open(os.path.join(odir, "README"), "w") as f:
        f.write(README)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Writes in FlatZinc the Golomb ruler problems with n = 5, ..., 10 marks: the
marks m1 = 0 < m2 < ... < mn in 0 .. n^2, all the differences between two
of them distinct, and the length mn as short as possible.

An instance is a variable m<i> per mark, m1 fixed to 0 by its domain, with
int_lt( m<i> , m<i+1> ); a variable d<i>_<j> in 1 .. n^2 per pair i < j,
with int_lin_eq( [ 1 , -1 , -1 ] , [ m<j> , m<i> , d<i>_<j> ] , 0 ); a
fzn_all_different_int over the differences; and int_lt( d1_2 , d<n-1>_<n> ),
which breaks the symmetry of the mirror ruler. The solve item minimizes mn.

usage: gen_golomb.py [ <output directory> ]

The output directory is fzn/golomb next to this script if none is given; it
gets the instances, optima.txt (the instance and the length of the optimal
ruler, one per line) and the README of the family.
"""

import os
import sys

# the length of the optimal Golomb ruler of n marks
OPT = {5: 11, 6: 17, 7: 25, 8: 34, 9: 44, 10: 55}

README = """\
data/fzn/golomb: Golomb rulers in FlatZinc
==========================================

The Golomb ruler problems with n = 5, ..., 10 marks (problem 6 of CSPLib,
https://www.csplib.org/Problems/prob006/): n marks in 0 .. n^2, the first at
0, whose pairwise differences are all distinct, the last one as small as
possible. The differences are variables linked to the marks by int_lin_eq
and held by one fzn_all_different_int, and the first difference is smaller
than the last one, which breaks the symmetry of the mirror ruler.
optima.txt gives the length of the optimal ruler of each n, as CSPLib
gives it.

The instances are written by gen_golomb.py in data/.
"""


def write_fzn(n):
    L = n * n
    pairs = [(i, j) for i in range(1, n + 1) for j in range(i + 1, n + 1)]
    out = ["% the Golomb ruler of {} marks: minimize its length".format(n)]
    out.append("var 0..0: m1;")
    for i in range(2, n + 1):
        out.append("var 1..{}: m{};".format(L, i))
    for i, j in pairs:
        out.append("var 1..{}: d{}_{};".format(L, i, j))
    for i in range(1, n):
        out.append("constraint int_lt(m{},m{});".format(i, i + 1))
    for i, j in pairs:
        out.append("constraint int_lin_eq([1,-1,-1],[m{},m{},d{}_{}],0);"
                   .format(j, i, i, j))
    out.append("constraint fzn_all_different_int([{}]);".format(
        ",".join("d{}_{}".format(i, j) for i, j in pairs)))
    out.append("constraint int_lt(d1_2,d{}_{});".format(n - 1, n))
    out.append("solve minimize m{};".format(n))
    return "\n".join(out) + "\n"


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    odir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "fzn",
                                                               "golomb")
    os.makedirs(odir, exist_ok=True)

    rows = []
    for n, z in OPT.items():
        name = "golomb{}".format(n)
        with open(os.path.join(odir, name + ".fzn"), "w") as f:
            f.write(write_fzn(n))
        rows.append("{}.fzn {}".format(name, z))

    with open(os.path.join(odir, "optima.txt"), "w") as f:
        f.write("# the instance and the length of its optimal ruler\n")
        f.write("\n".join(rows) + "\n")
    with open(os.path.join(odir, "README"), "w") as f:
        f.write(README)


if __name__ == "__main__":
    main()

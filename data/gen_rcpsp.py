#!/usr/bin/env python3
"""
Writes in FlatZinc a selection of the single-mode RCPSP instances of the set
j30 of PSPLIB, 30 jobs and 4 renewable resources each, with the optimal
makespans of the set (j30opt.sm). The selection is the first instance of
each parameter set p with p mod 5 equal to 1 or 3, i.e., 20 of the 48 sets,
spread over all the levels of network complexity, resource factor and
resource strength of the set.

An instance is a variable s<j> per non-dummy job j, its start in
0 .. H - d_j, H being the horizon of the instance (the sum of the
durations) and d_j the duration of the job, and the variable makespan in
0 .. H, which the solve item minimizes; a precedence i -> j is
int_lin_le( [ 1 , -1 ] , [ s<i> , s<j> ] , - d_i ), the makespan following
the predecessors of the dummy sink alike, and each resource is a
fzn_cumulative over the jobs that use it, with the durations, the usages and
the capacity as parameters.

usage: gen_rcpsp.py [ <output directory> ]

The output directory is fzn/rcpsp next to this script if none is given; it
gets the instances, optima.txt (the instance and its optimal makespan, one
per line) and the README of the family. The instances and the optima are
downloaded from PSPLIB.
"""

import io
import os
import sys
import tarfile
import urllib.request

SITE = "https://www.om-db.wi.tum.de/psplib/"
DATA = SITE + "download_dataset.php?set=j30&mode=sm&format=tgz"
OPT = SITE + "download_solution.php?set=j30&mode=sm&type=opt"

SETS = [p for p in range(1, 49) if p % 5 in (1, 3)]

README = """\
data/fzn/rcpsp: single-mode RCPSP instances of PSPLIB in FlatZinc
==================================================================

20 instances of the set j30 of PSPLIB (https://www.om-db.wi.tum.de/psplib/),
30 jobs and 4 renewable resources each: the first instance of each parameter
set p with p mod 5 equal to 1 or 3, which spreads them over all the levels
of network complexity, resource factor and resource strength of the set.

Each instance minimizes the makespan under the precedences, written as
int_lin_le, and one fzn_cumulative per resource, whose durations, usages and
capacity are parameters. optima.txt gives the optimal makespan of each, as
j30opt.sm of PSPLIB gives it.

The instances are written by gen_rcpsp.py in data/. PSPLIB asks to cite

  R. Kolisch, A. Sprecher, PSPLIB - A project scheduling problem library,
  European Journal of Operational Research 96(1), 205-216, 1997.
"""


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read()


def parse_sm(text):
    """the durations, the usages, the successors, the capacities and the
    horizon of a .sm file of PSPLIB"""
    lines = text.splitlines()
    horizon = None
    succ, dur, use, cap = {}, {}, {}, []
    sec = None
    for ln in lines:
        s = ln.strip()
        if s.startswith("horizon"):
            horizon = int(s.split(":")[1])
        elif s.startswith("PRECEDENCE RELATIONS"):
            sec = "prec"
        elif s.startswith("REQUESTS/DURATIONS"):
            sec = "req"
        elif s.startswith("RESOURCEAVAILABILITIES"):
            sec = "cap"
        elif s.startswith("*"):
            sec = None
        elif not s:
            continue
        elif sec == "prec" and s[0].isdigit():
            f = [int(x) for x in s.split()]
            succ[f[0]] = f[3:3 + f[2]]
        elif sec == "req" and s[0].isdigit():
            f = [int(x) for x in s.split()]
            dur[f[0]] = f[2]
            use[f[0]] = f[3:]
        elif sec == "cap" and s[0].isdigit():
            cap = [int(x) for x in s.split()]
    return dur, use, succ, cap, horizon


def write_fzn(name, dur, use, succ, cap, H):
    jobs = sorted(dur)
    src, snk = jobs[0], jobs[-1]
    real = [j for j in jobs if j not in (src, snk)]
    out = io.StringIO()
    out.write("% {} of the set j30 of PSPLIB, single mode: minimize the "
              "makespan\n".format(name))
    for j in real:
        out.write("var 0..{}: s{};\n".format(H - dur[j], j))
    out.write("var 0..{}: makespan;\n".format(H))
    for i in real:
        for j in succ[i]:
            other = "makespan" if j == snk else "s{}".format(j)
            out.write("constraint int_lin_le([1,-1],[s{},{}],{});\n"
                      .format(i, other, -dur[i]))
    for k, c in enumerate(cap):
        t = [j for j in real if use[j][k] > 0]
        if not t:
            continue
        out.write("constraint fzn_cumulative([{}],[{}],[{}],{});\n".format(
            ",".join("s{}".format(j) for j in t),
            ",".join(str(dur[j]) for j in t),
            ",".join(str(use[j][k]) for j in t), c))
    out.write("solve minimize makespan;\n")
    return out.getvalue()


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    odir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "fzn",
                                                               "rcpsp")
    os.makedirs(odir, exist_ok=True)

    opt = {}
    for ln in fetch(OPT).decode("latin-1").splitlines():
        f = ln.split()
        if len(f) >= 3 and f[0].isdigit() and f[1].isdigit():
            opt[(int(f[0]), int(f[1]))] = int(f[2])

    tar = tarfile.open(fileobj=io.BytesIO(fetch(DATA)), mode="r:gz")
    files = {os.path.basename(m.name): m for m in tar.getmembers()
             if m.isfile()}

    rows = []
    for p in SETS:
        name = "j30{}_1".format(p)
        text = tar.extractfile(files[name + ".sm"]).read().decode("latin-1")
        dur, use, succ, cap, H = parse_sm(text)
        with open(os.path.join(odir, name + ".fzn"), "w") as f:
            f.write(write_fzn(name, dur, use, succ, cap, H))
        rows.append("{}.fzn {}".format(name, opt[(p, 1)]))

    with open(os.path.join(odir, "optima.txt"), "w") as f:
        f.write("# the instance and its optimal makespan [j30opt.sm of "
                "PSPLIB]\n")
        f.write("\n".join(rows) + "\n")
    with open(os.path.join(odir, "README"), "w") as f:
        f.write(README)


if __name__ == "__main__":
    main()

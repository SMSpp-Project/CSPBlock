#!/usr/bin/env python3
"""
Writes in FlatZinc the minimum colouring of some graphs of the DIMACS
graph colouring benchmark whose chromatic number is known: the Mycielski
graphs myciel3, myciel4 and myciel5, the queen graphs queen5_5, queen6_6 and
queen7_7, the book graphs anna, david, homer, huck and jean, the game graph
games120 and the mile graphs miles250 and miles500.

An instance is a variable x<v> per vertex v, its colour in 1 .. U, U being
the largest degree of a vertex plus 1, and the variable
k in q .. U, q being the size of a clique found greedily, which the solve
item minimizes; each edge u v is int_ne( x<u> , x<v> ), and each vertex is
int_le( x<v> , k ). The vertices of the clique have the colours 1 .. q,
fixed by their domain, which breaks the symmetry of the colours and cuts
no optimal colouring.

usage: gen_colouring.py [ <output directory> ]

The output directory is fzn/colouring next to this script if none is given;
it gets the instances, optima.txt (the instance and its chromatic number,
one per line) and the README of the family. The graphs are downloaded from
the site of the benchmark.
"""

import os
import sys
import urllib.request

SITE = "https://mat.tepper.cmu.edu/COLOR/instances/"

# the chromatic number of each graph [see the site]
CHI = {"myciel3": 4, "myciel4": 5, "myciel5": 6,
       "queen5_5": 5, "queen6_6": 7, "queen7_7": 7,
       "anna": 11, "david": 11, "homer": 13, "huck": 11, "jean": 10,
       "games120": 9, "miles250": 8, "miles500": 20}

README = """\
data/fzn/colouring: minimum graph colouring in FlatZinc
=======================================================

14 graphs of the DIMACS graph colouring benchmark
(https://mat.tepper.cmu.edu/COLOR/instances.html) whose chromatic number is
known: the Mycielski graphs myciel3, myciel4 and myciel5, the queen graphs
queen5_5, queen6_6 and queen7_7, the book graphs anna, david, homer, huck
and jean and the graphs games120, miles250 and miles500 of the Stanford
GraphBase.

Each instance minimizes the number k of colours, with int_ne between the
colours of the ends of each edge and int_le between the colour of each
vertex and k; the colours of a clique found greedily are fixed to 1 .. q,
which breaks their symmetry. optima.txt gives the chromatic number of each
graph, as the site of the benchmark gives it.

The instances are written by gen_colouring.py in data/. The benchmark is
described in

  D.S. Johnson, M.A. Trick (eds.), Cliques, Coloring, and Satisfiability:
  Second DIMACS Implementation Challenge, DIMACS Series in Discrete
  Mathematics and Theoretical Computer Science 26, AMS, 1996.

and the book, game and mile graphs come from

  D.E. Knuth, The Stanford GraphBase: A Platform for Combinatorial
  Computing, ACM Press, 1993.
"""


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode("latin-1")


def parse_col(text):
    """the number of vertices and the edges, each once, of a .col file"""
    n, edges = 0, set()
    for ln in text.splitlines():
        f = ln.split()
        if not f:
            continue
        if f[0] == "p":
            n = int(f[2])
        elif f[0] == "e":
            u, v = int(f[1]), int(f[2])
            if u != v:
                edges.add((min(u, v), max(u, v)))
    return n, sorted(edges)


def clique(n, adj):
    """a clique, grown greedily from each vertex by the largest degree"""
    best = []
    for s in range(1, n + 1):
        c, cand = [s], set(adj[s])
        while cand:
            v = max(cand, key=lambda v: (len(adj[v] & cand), -v))
            c.append(v)
            cand &= adj[v]
        if len(c) > len(best):
            best = c
    return best


def write_fzn(name, n, edges):
    adj = {v: set() for v in range(1, n + 1)}
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    U = max(len(a) for a in adj.values()) + 1
    q = clique(n, adj)
    fixed = {v: c + 1 for c, v in enumerate(q)}
    out = [("% {} of the DIMACS graph colouring benchmark: minimize the "
            "number of colours").format(name)]
    for v in range(1, n + 1):
        d = fixed.get(v)
        out.append("var {}..{}: x{};".format(d or 1, d or U, v))
    out.append("var {}..{}: k;".format(len(q), U))
    for u, v in edges:
        out.append("constraint int_ne(x{},x{});".format(u, v))
    for v in range(1, n + 1):
        out.append("constraint int_le(x{},k);".format(v))
    out.append("solve minimize k;")
    return "\n".join(out) + "\n"


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    odir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "fzn",
                                                               "colouring")
    os.makedirs(odir, exist_ok=True)

    rows = []
    for name, chi in CHI.items():
        n, edges = parse_col(fetch(SITE + name + ".col"))
        with open(os.path.join(odir, name + ".fzn"), "w") as f:
            f.write(write_fzn(name, n, edges))
        rows.append("{}.fzn {}".format(name, chi))

    with open(os.path.join(odir, "optima.txt"), "w") as f:
        f.write("# the instance and the chromatic number of its graph\n")
        f.write("\n".join(rows) + "\n")
    with open(os.path.join(odir, "README"), "w") as f:
        f.write(README)


if __name__ == "__main__":
    main()

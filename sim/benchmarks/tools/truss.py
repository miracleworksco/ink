#!/usr/bin/env python3
"""Reference for benchmark 12: pin-jointed Pratt truss by the direct stiffness
method (numpy-free, small dense solve). Prints member forces, stresses and the
midspan deflection under the benchmark's load.

usage: python3 benchmarks/tools/truss.py
"""
import math

E, r, P = 200000.0, 3.0, 600.0          # MPa, member radius mm, N at the centre bottom node
A = math.pi * r * r
panel, h, n = 60.0, 60.0, 6
nodes = {f"b{i}": (i * panel, 0.0) for i in range(n + 1)}
nodes.update({f"t{i}": (i * panel, h) for i in range(1, n)})
members = [(f"b{i}", f"b{i+1}") for i in range(n)]                    # bottom chord
members += [(f"t{i}", f"t{i+1}") for i in range(1, n - 1)]            # top chord
members += [("b0", "t1"), (f"b{n}", f"t{n-1}")]                         # end posts
members += [(f"b{i}", f"t{i}") for i in range(1, n)]                  # verticals
members += [("t1", "b2"), ("t2", "b3"), ("t4", "b3"), ("t5", "b4")]   # diagonals, sloping to the centre

names = list(nodes)
idx = {k: i for i, k in enumerate(names)}
N = 2 * len(names)
K = [[0.0] * N for _ in range(N)]
geo = []
for a, b in members:
    (xa, ya), (xb, yb) = nodes[a], nodes[b]
    L = math.hypot(xb - xa, yb - ya)
    c, s = (xb - xa) / L, (yb - ya) / L
    k = E * A / L
    geo.append((a, b, L, c, s))
    dofs = [2 * idx[a], 2 * idx[a] + 1, 2 * idx[b], 2 * idx[b] + 1]
    v = [-c, -s, c, s]
    for i in range(4):
        for j in range(4):
            K[dofs[i]][dofs[j]] += k * v[i] * v[j]
F = [0.0] * N
F[2 * idx["b3"] + 1] = -P
fixed = {2 * idx["b0"], 2 * idx["b0"] + 1, 2 * idx[f"b{n}"] + 1}       # pin + roller
free = [i for i in range(N) if i not in fixed]
# Gaussian elimination on the free system
M = [[K[i][j] for j in free] + [F[i]] for i in free]
m = len(free)
for c in range(m):
    p = max(range(c, m), key=lambda r_: abs(M[r_][c]))
    M[c], M[p] = M[p], M[c]
    for r_ in range(m):
        if r_ != c and M[r_][c] != 0:
            f = M[r_][c] / M[c][c]
            M[r_] = [x - f * y for x, y in zip(M[r_], M[c])]
u = [0.0] * N
for i, dof in enumerate(free):
    u[dof] = M[i][m] / M[i][i]
print(f"A = {A:.3f} mm^2")
for a, b, L, c, s in geo:
    du = [u[2 * idx[b]] - u[2 * idx[a]], u[2 * idx[b] + 1] - u[2 * idx[a] + 1]]
    force = E * A / L * (c * du[0] + s * du[1])
    print(f"{a:>3}-{b:<3} L {L:6.1f}  force {force:8.1f} N  stress {force / A:7.2f} MPa")
print(f"midspan deflection (b3): {u[2 * idx['b3'] + 1]:.4f} mm")

# A bar's solid answer, by CalculiX (ccx, or $CCX), C3D20R on a structured mesh: the truth of the holds
# table's beam rows (truths.json, "CalculiX, C3D20R ..."), in place of a beam formula.
# usage: solid.py CASE.case [NL,NB,NT ...]   (default meshes 40,2,2 80,4,4 160,8,8 320,12,12; RUNS: work folder)
# It reads the case: its box (L along x, B, T), moved and turned about z; E, nu, density; each fix at an end
# (its region the slab at x = 0 or x = L) and what it leaves free; the load at x = L (a force along -x, as a
# uniform pressure on that face); and its claim: the first buckling factor, or with (modes n) the first
# frequency (Hz). A fix with no free motion holds every node of its end face. A pivot or slide holds its face's
# mean motion, as ink does: of the face's rigid motions, those left after its free ones are taken out
# W-orthogonally (W the face's area), each as one equation, the integral of psi . u over the face equal to 0
# (ccx 2.23's *BUCKLE ignores *DISTRIBUTING, so not that).
import sys, os, re, math, json, subprocess, time

CCX = os.environ.get("CCX", "ccx")


def read_case(path):
    s = open(path).read()
    L, B, T = map(float, re.search(r"\(box ([-0-9.e]+) ([-0-9.e]+) ([-0-9.e]+)\)", s).groups())
    ang = float(re.search(r"\(rotate-z \(move \(box [^)]*\) [^)]*\) ([-0-9.e]+)\)", s).group(1))
    tx, ty, tz = map(float, re.search(r"\(rotate-z \(move \(box [^)]*\) [^)]*\) [-0-9.e]+\) ([-0-9.e]+) ([-0-9.e]+) ([-0-9.e]+)\)", s).groups())
    c, sn = math.cos(math.radians(ang)), math.sin(math.radians(ang))
    body = lambda p: ((p[0] - tx) * c + (p[1] - ty) * sn, -(p[0] - tx) * sn + (p[1] - ty) * c, p[2] - tz)
    nu = float(re.search(r"\(nu ([-0-9.e]+)\)", s).group(1))
    E = float(re.search(r"\(E ([-0-9.e]+)MPa\)", s).group(1))
    rho = re.search(r"\(density ([-0-9.e]+)kg/m3\)", s)
    modes = re.search(r"\(modes (\d+)\)", s)
    fixes = {}
    for line in s.splitlines():
        m = re.match(r"\s*\(fix (\S+) (.*)", line)
        if not m:
            continue
        end = "far" if re.search(r"\(- [0-9.]+9+ \(\+", line) else "near"
        piv = re.search(r"\(pivot \(at ([-0-9.e]+) ([-0-9.e]+) ([-0-9.e]+)\) \(axis ([-0-9.e]+) ([-0-9.e]+) ([-0-9.e]+)\)\)", line)
        slide = "(slide" in line
        free = []
        if piv:
            at = body(tuple(map(float, piv.groups()[:3])))
            ax = tuple(map(float, piv.groups()[3:]))
            ax = (ax[0] * c + ax[1] * sn, -ax[0] * sn + ax[1] * c, ax[2])
            free.append(("turn", at, ax))
        if slide:
            free.append(("slide", None, (1.0, 0.0, 0.0)))
        fixes[end] = free
    load = re.search(r"\(force ([-0-9.e]+)N", s)
    return dict(L=L, B=B, T=T, E=E, nu=nu, rho=float(rho.group(1)) * 1e-12 if rho else None,
                modes=int(modes.group(1)) if modes else 0, fixes=fixes, load=-float(load.group(1)) if load else 0.0)


def solve(p, nl, nb, nt, wd):
    L, B, T = p["L"], p["B"], p["T"]
    nid, X = {}, []
    def node(i, j, k):
        if (i, j, k) not in nid:
            nid[(i, j, k)] = len(X) + 1
            X.append((L * i / (2 * nl), -B / 2 + B * j / (2 * nb), -T / 2 + T * k / (2 * nt)))
        return nid[(i, j, k)]
    E = []
    for a in range(nl):
        for b in range(nb):
            for cc in range(nt):
                i, j, k = 2 * a, 2 * b, 2 * cc
                cor = [(i, j, k), (i + 2, j, k), (i + 2, j + 2, k), (i, j + 2, k), (i, j, k + 2), (i + 2, j, k + 2), (i + 2, j + 2, k + 2), (i, j + 2, k + 2)]
                mid = [(i + 1, j, k), (i + 2, j + 1, k), (i + 1, j + 2, k), (i, j + 1, k), (i + 1, j, k + 2), (i + 2, j + 1, k + 2), (i + 1, j + 2, k + 2), (i, j + 1, k + 2),
                       (i, j, k + 1), (i + 2, j, k + 1), (i + 2, j + 2, k + 1), (i, j + 2, k + 1)]
                E.append((a, [node(*q) for q in cor + mid]))
    lines = ["*NODE, NSET=NALL"] + [f"{i+1},{x:.15g},{y:.15g},{z:.15g}" for i, (x, y, z) in enumerate(X)]
    lines.append("*ELEMENT, TYPE=C3D20R, ELSET=EALL")
    lines += [f"{e+1}," + ",".join(map(str, ns[:15])) + ",\n" + ",".join(map(str, ns[15:])) for e, (_, ns) in enumerate(E)]
    lines.append(f"*MATERIAL, NAME=M\n*ELASTIC\n{p['E']},{p['nu']}" + (f"\n*DENSITY\n{p['rho']:.6e}" if p["rho"] else ""))
    lines.append("*SOLID SECTION, ELSET=EALL, MATERIAL=M")
    g, gw = [-math.sqrt(0.6), 0.0, math.sqrt(0.6)], [5 / 9, 8 / 9, 5 / 9]
    def n8(r, s):
        cc = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        out = [0.25 * (1 + r * a) * (1 + s * b) * (r * a + s * b - 1) for a, b in cc]
        return out + [0.5 * (1 - r * r) * (1 - s), 0.5 * (1 + r) * (1 - s * s), 0.5 * (1 - r * r) * (1 + s), 0.5 * (1 - r) * (1 - s * s)]
    bc, eqs = [], []
    for end, free in p["fixes"].items():
        i0, x0 = (0, 0.0) if end == "near" else (2 * nl, L)
        if not free:
            bc += [f"{n},1,3,0" for (i, j, k), n in nid.items() if i == i0]
            continue
        # the face's weights: node -> (int N, int N y, int N z)
        w = {}
        for b in range(nb):
            for cc in range(nt):
                j, k = 2 * b, 2 * cc
                quad = [(j, k), (j + 2, k), (j + 2, k + 2), (j, k + 2), (j + 1, k), (j + 2, k + 1), (j + 1, k + 2), (j, k + 1)]
                hy, hz = B / nb / 2, T / nt / 2
                yc, zc = -B / 2 + B * (j + 1) / (2 * nb), -T / 2 + T * (k + 1) / (2 * nt)
                for r, wr in zip(g, gw):
                    for t, wt in zip(g, gw):
                        N = n8(r, t); y, z = yc + r * hy, zc + t * hz; dA = wr * wt * hy * hz
                        for q, (jj, kk) in enumerate(quad):
                            v = w.setdefault(nid[(i0, jj, kk)], [0.0, 0.0, 0.0])
                            v[0] += N[q] * dA; v[1] += N[q] * y * dA; v[2] += N[q] * z * dA
        # motions on the face (x = x0), each component as (a, b, c): a + b y + c z
        def motion(om, at, v):
            # om x (p - at) + v at p = (x0, y, z)
            r0 = (x0 - at[0], -at[1], -at[2])
            comp = []
            for d in range(3):
                e1, e2 = (d + 1) % 3, (d + 2) % 3
                a = om[e1] * r0[e2] - om[e2] * r0[e1] + v[d]
                # om x (0, y, z): y and z parts
                by = om[e1] * (1 if e2 == 1 else 0) - om[e2] * (1 if e1 == 1 else 0)
                cz = om[e1] * (1 if e2 == 2 else 0) - om[e2] * (1 if e1 == 2 else 0)
                comp.append((a, by, cz))
            return comp
        z3 = (0.0, 0.0, 0.0)
        rigid = [motion(z3, z3, (1, 0, 0)), motion(z3, z3, (0, 1, 0)), motion(z3, z3, (0, 0, 1)),
                 motion((1, 0, 0), (x0, 0, 0), z3), motion((0, 1, 0), (x0, 0, 0), z3), motion((0, 0, 1), (x0, 0, 0), z3)]
        frees = [motion(ax, at, z3) if kind == "turn" else motion(z3, z3, ax) for kind, at, ax in free]
        A, Iy, Iz = B * T, B ** 3 * T / 12, B * T ** 3 / 12
        def inner(m1, m2):  # int m1 . m2 dA over the face (centroid at y = z = 0)
            return sum(a1 * a2 * A + b1 * b2 * Iy + c1 * c2 * Iz for (a1, b1, c1), (a2, b2, c2) in zip(m1, m2))
        def comb(ms, cs):
            return [tuple(sum(cf * m[d][q] for m, cf in zip(ms, cs)) for q in range(3)) for d in range(3)]
        basis = []
        def orth(m):
            for o in basis:
                k = inner(m, o)
                m = [tuple(m[d][q] - k * o[d][q] for q in range(3)) for d in range(3)]
            n = inner(m, m)
            return m, n
        for f in frees:
            m, n = orth(f)
            if n > 1e-12 * inner(f, f):
                basis.append([tuple(v / math.sqrt(n) for v in comp) for comp in m])
        nfree = len(basis)
        for r_ in rigid:
            m, n = orth(r_)
            if n > 1e-9 * inner(r_, r_):
                basis.append([tuple(v / math.sqrt(n) for v in comp) for comp in m])
        for m in basis[nfree:]:
            eq = [(n, d + 1, m[d][0] * v[0] + m[d][1] * v[1] + m[d][2] * v[2]) for n, v in w.items() for d in range(3)]
            eqs.append([t for t in eq if abs(t[2]) > 1e-14 * max(abs(x[2]) for x in eq)])
    used = set()
    for eq in eqs:
        # (its dependent dof: one no other equation's terms touch, else the first free one; a dof two equations
        # shared as their dependents' left ccx's buckling with a spurious factor near 1)
        others = {(n, d) for e2 in eqs if e2 is not eq for n, d, _ in e2}
        first = next((t for t in eq if (t[0], t[1]) not in others and (t[0], t[1]) not in used), None) or next(t for t in eq if (t[0], t[1]) not in used)
        used.add((first[0], first[1]))
        eq = [first] + [t for t in eq if t is not first]
        lines.append("*EQUATION"); lines.append(str(len(eq)))
        terms = [f"{n},{d},{cf:.15g}" for n, d, cf in eq]
        lines += [",".join(terms[i:i + 4]) for i in range(0, len(terms), 4)]
    if bc:
        lines.append("*BOUNDARY"); lines += bc
    if p["modes"]:
        lines.append("*STEP\n*FREQUENCY\n3")
    else:
        hi = [e + 1 for e, (a, _) in enumerate(E) if a == nl - 1]
        lines.append("*STEP\n*BUCKLE\n2\n*DLOAD")
        lines += [f"{e},P4,{p['load'] / (B * T):.15g}" for e in hi]
    lines.append("*NODE FILE\nU\n*END STEP")
    os.makedirs(wd, exist_ok=True)
    open(os.path.join(wd, "job.inp"), "w").write("\n".join(lines) + "\n")
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("THREADS", "8"), CCX_NPROC_EQUATION_SOLVER=os.environ.get("THREADS", "8"))
    t0 = time.time()
    r = subprocess.run([CCX, "-i", "job"], cwd=wd, env=env, capture_output=True, text=True)
    dt = time.time() - t0
    txt = open(os.path.join(wd, "job.dat")).read() if os.path.exists(os.path.join(wd, "job.dat")) else ""
    vals = []
    if p["modes"] and "E I G E N V A L U E   O U T P U T" in txt:
        for ln in txt.split("E I G E N V A L U E   O U T P U T")[1].splitlines():
            s_ = ln.split()
            if len(s_) == 5 and s_[0].isdigit():
                vals.append(float(s_[3]))
    if not p["modes"] and "B U C K L I N G" in txt:
        for ln in txt.split("B U C K L I N G")[1].splitlines():
            s_ = ln.split()
            if len(s_) == 2 and s_[0].isdigit():
                vals.append(float(s_[1]))
    return dict(mesh=[nl, nb, nt], nodes=len(X), value=vals[0] if vals else None, s=round(dt, 1), err=None if vals else (r.stdout + r.stderr)[-600:])


if __name__ == "__main__":
    case = sys.argv[1]
    p = read_case(case)
    meshes = sys.argv[2:] or ["40,2,2", "80,4,4", "160,8,8", "320,12,12"]
    for m in meshes:
        nl, nb, nt = map(int, m.split(","))
        r = solve(p, nl, nb, nt, os.path.join(os.environ.get("RUNS", "runs"), f"{os.path.basename(case)}_{nl}_{nb}_{nt}"))
        print(json.dumps(dict(case=os.path.basename(case), **r)), flush=True)

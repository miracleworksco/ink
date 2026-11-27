#!/usr/bin/env python3
"""The MBB half-beam of Andreassen et al. 2011 (their 88-line code's example:
60 x 20 Q4 elements, volume fraction 0.5, penalty 3, filter radius 2.4, a unit
load at the top of the symmetry line, rolled at the far bottom corner), the
reference of ink-ext/tests/topology.rs's topology_matches_the_88_line_codes_mbb_beam.

Prints each run's compliance over the solid domain's:
- the 88-line code ported to numpy (ft 1 sensitivity, ft 2 density, ft 3
  Heaviside filters);
- the same code with ink's scheme: density filter, tanh projection at 1/2
  sharpened 1 to 16 (doubling when settled or after a tenth of the budget),
  MMA (Svanberg) with the multiplier meeting the projected volume, stop when
  no physical density moves by 0.01 at 16.

Run: uv run --with numpy --with scipy python benchmarks/tools/top88_mbb.py
(Ryzen AI Max+ 395, 2026-11-01: solid 125.878; ft1 1.722, ft2 1.857, ft3 1.595;
ink's scheme at budgets 60, 100, 300: 1.542, 1.512, 1.509).
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def top88(nelx, nely, volfrac, penal, rmin, ft, maxit=2000):
    E0, Emin, nu = 1.0, 1e-9, 0.3
    A11 = np.array([[12, 3, -6, -3], [3, 12, 3, 0], [-6, 3, 12, -3], [-3, 0, -3, 12]])
    A12 = np.array([[-6, -3, 0, 3], [-3, -6, -3, -6], [0, -3, -6, 3], [3, -6, 3, -6]])
    B11 = np.array([[-4, 3, -2, 9], [3, -4, -9, 4], [-2, -9, -4, -3], [9, 4, -3, -4]])
    B12 = np.array([[2, -3, 4, -9], [-3, 2, 9, -2], [4, 9, 2, 3], [-9, -2, 3, 2]])
    KE = 1 / (1 - nu**2) / 24 * (np.block([[A11, A12], [A12.T, A11]]) + nu * np.block([[B11, B12], [B12.T, B11]]))
    nodenrs = np.arange((1 + nelx) * (1 + nely)).reshape(1 + nelx, 1 + nely).T
    edofVec = (2 * nodenrs[:-1, :-1] + 2).T.reshape(-1)  # 0-based of matlab's 2*n+1
    edofMat = edofVec[:, None] + np.array([0, 1, 2 * nely + 2, 2 * nely + 3, 2 * nely, 2 * nely + 1, -2, -1])
    iK = np.kron(edofMat, np.ones((8, 1))).flatten().astype(int)
    jK = np.kron(edofMat, np.ones((1, 8))).flatten().astype(int)
    ndof = 2 * (nelx + 1) * (nely + 1)
    F = np.zeros(ndof); F[1] = -1
    fixed = np.union1d(np.arange(0, 2 * (nely + 1), 2), [ndof - 1])
    free = np.setdiff1d(np.arange(ndof), fixed)
    # filter
    rows, cols, vals = [], [], []
    for i1 in range(nelx):
        for j1 in range(nely):
            e1 = i1 * nely + j1
            for i2 in range(max(i1 - int(np.ceil(rmin)) + 1, 0), min(i1 + int(np.ceil(rmin)), nelx)):
                for j2 in range(max(j1 - int(np.ceil(rmin)) + 1, 0), min(j1 + int(np.ceil(rmin)), nely)):
                    w = max(0, rmin - np.hypot(i1 - i2, j1 - j2))
                    rows.append(e1); cols.append(i2 * nely + j2); vals.append(w)
    H = coo_matrix((vals, (rows, cols)), shape=(nelx * nely, nelx * nely)).tocsr()
    Hs = np.asarray(H.sum(1)).ravel()
    x = volfrac * np.ones(nelx * nely)
    beta = 1.0
    if ft in (1, 2):
        xPhys = x.copy()
    else:
        xTilde = x.copy(); xPhys = 1 - np.exp(-beta * xTilde) + xTilde * np.exp(-beta)
    loop, change, loopbeta = 0, 1.0, 0
    while change > 0.01 and loop < maxit:
        loop += 1; loopbeta += 1
        sK = (KE.flatten()[None, :] * (Emin + xPhys[:, None] ** penal * (E0 - Emin))).flatten()
        K = coo_matrix((sK, (iK, jK)), shape=(ndof, ndof)).tocsc()
        U = np.zeros(ndof)
        U[free] = spsolve(K[free][:, free], F[free])
        ce = np.einsum('ij,jk,ik->i', U[edofMat], KE, U[edofMat])
        c = ((Emin + xPhys ** penal * (E0 - Emin)) * ce).sum()
        dc = -penal * (E0 - Emin) * xPhys ** (penal - 1) * ce
        dv = np.ones(nelx * nely)
        if ft == 1:
            dc = H @ (x * dc) / Hs / np.maximum(1e-3, x)
        elif ft == 2:
            dc = H @ (dc / Hs); dv = H @ (dv / Hs)
        else:
            dx = beta * np.exp(-beta * xTilde) + np.exp(-beta)
            dc = H @ (dc * dx / Hs); dv = H @ (dv * dx / Hs)
        l1, l2, move = 0.0, 1e9, 0.2
        while (l2 - l1) / (l1 + l2) > 1e-3:
            lmid = 0.5 * (l2 + l1)
            xnew = np.maximum(0, np.maximum(x - move, np.minimum(1, np.minimum(x + move, x * np.sqrt(-dc / dv / lmid)))))
            if ft == 1:
                xPhys = xnew
            elif ft == 2:
                xPhys = (H @ xnew) / Hs
            else:
                xTilde = (H @ xnew) / Hs
                xPhys = 1 - np.exp(-beta * xTilde) + xTilde * np.exp(-beta)
            if xPhys.sum() > volfrac * nelx * nely:
                l1 = lmid
            else:
                l2 = lmid
        change = np.abs(xnew - x).max()
        x = xnew
        if ft == 3 and beta < 512 and (loopbeta >= 50 or change <= 0.01):
            beta *= 2; loopbeta = 0; change = 1.0
    grey = np.mean(4 * xPhys * (1 - xPhys))
    return c, loop, grey



def setup(nelx, nely, rmin):
    nu = 0.3
    A11 = np.array([[12, 3, -6, -3], [3, 12, 3, 0], [-6, 3, 12, -3], [-3, 0, -3, 12]])
    A12 = np.array([[-6, -3, 0, 3], [-3, -6, -3, -6], [0, -3, -6, 3], [3, -6, 3, -6]])
    B11 = np.array([[-4, 3, -2, 9], [3, -4, -9, 4], [-2, -9, -4, -3], [9, 4, -3, -4]])
    B12 = np.array([[2, -3, 4, -9], [-3, 2, 9, -2], [4, 9, 2, 3], [-9, -2, 3, 2]])
    KE = 1 / (1 - nu**2) / 24 * (np.block([[A11, A12], [A12.T, A11]]) + nu * np.block([[B11, B12], [B12.T, B11]]))
    nodenrs = np.arange((1 + nelx) * (1 + nely)).reshape(1 + nelx, 1 + nely).T
    edofVec = (2 * nodenrs[:-1, :-1] + 2).T.reshape(-1)
    edofMat = edofVec[:, None] + np.array([0, 1, 2 * nely + 2, 2 * nely + 3, 2 * nely, 2 * nely + 1, -2, -1])
    iK = np.kron(edofMat, np.ones((8, 1))).flatten().astype(int)
    jK = np.kron(edofMat, np.ones((1, 8))).flatten().astype(int)
    ndof = 2 * (nelx + 1) * (nely + 1)
    F = np.zeros(ndof); F[1] = -1
    fixed = np.union1d(np.arange(0, 2 * (nely + 1), 2), [ndof - 1])
    free = np.setdiff1d(np.arange(ndof), fixed)
    rows, cols, vals = [], [], []
    for i1 in range(nelx):
        for j1 in range(nely):
            e1 = i1 * nely + j1
            for i2 in range(max(i1 - int(np.ceil(rmin)) + 1, 0), min(i1 + int(np.ceil(rmin)), nelx)):
                for j2 in range(max(j1 - int(np.ceil(rmin)) + 1, 0), min(j1 + int(np.ceil(rmin)), nely)):
                    w = max(0, rmin - np.hypot(i1 - i2, j1 - j2))
                    rows.append(e1); cols.append(i2 * nely + j2); vals.append(w)
    H = coo_matrix((vals, (rows, cols)), shape=(nelx * nely, nelx * nely)).tocsr()
    Hs = np.asarray(H.sum(1)).ravel()
    return KE, edofMat, iK, jK, ndof, F, free, H, Hs

def solve(KE, edofMat, iK, jK, ndof, F, free, E):
    sK = (KE.flatten()[None, :] * E[:, None]).flatten()
    K = coo_matrix((sK, (iK, jK)), shape=(ndof, ndof)).tocsc()
    U = np.zeros(ndof)
    U[free] = spsolve(K[free][:, free], F[free])
    return U

class Mma:
    def __init__(s): s.k = 0
    def update(s, x, df, dg, step, over):
        n = len(x)
        if s.k < 2:
            s.low, s.upp = x - 0.5, x + 0.5
        else:
            zz = (x - s.x1) * (s.x1 - s.x2)
            gam = np.where(zz < 0, 0.7, np.where(zz > 0, 1.2, 1.0))
            s.low = np.clip(x - gam * (s.x1 - s.low), x - 10, x - 0.01)
            s.upp = np.clip(x + gam * (s.upp - s.x1), x + 0.01, x + 10)
        low, upp = s.low, s.upp
        alfa = np.maximum(np.maximum(low + 0.1 * (x - low), x - step), 0)
        beta = np.minimum(np.minimum(upp - 0.1 * (upp - x), x + step), 1)
        def pq(d):
            p, q = np.maximum(d, 0), np.maximum(-d, 0); r = 0.001 * (p + q) + 1e-5
            return (p + r) * (upp - x) ** 2, (q + r) * (x - low) ** 2
        p0, q0 = pq(df); p1, q1 = pq(dg)
        at = lambda lam: np.clip((np.sqrt(p0 + lam * p1) * low + np.sqrt(q0 + lam * q1) * upp) / (np.sqrt(p0 + lam * p1) + np.sqrt(q0 + lam * q1)), alfa, beta)
        gt = lambda y: 1.0 if over(y) else -1.0
        y = at(0.0)
        if gt(y) > 0:
            lo, hi = 0.0, 1.0
            while gt(at(hi)) > 0 and hi < 1e30: lo, hi = hi, hi * 10
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if gt(at(mid)) > 0: lo = mid
                else: hi = mid
                if hi - lo <= 1e-9 * hi: break
            y = at(hi)
        s.x2 = getattr(s, 'x1', x.copy()); s.x1 = x.copy(); s.k += 1
        return y

def ink_scheme(nelx, nely, volfrac, rmin, maxit, penal=3.0, emin=1e-9):
    KE, edofMat, iK, jK, ndof, F, free, H, Hs = setup(nelx, nely, rmin)
    n = nelx * nely
    x = volfrac * np.ones(n)
    proj = lambda b, t: ((np.tanh(b * 0.5) + np.tanh(b * (t - 0.5))) / (2 * np.tanh(b * 0.5)), b * (1 - np.tanh(b * (t - 0.5)) ** 2) / (2 * np.tanh(b * 0.5)))
    beta, every, since, prev, mma, c0 = 1.0, max(maxit // 10, 1), 0, None, Mma(), None
    for it in range(1, maxit + 1):
        xt = H @ x / Hs
        xb, dxb = proj(beta, xt)
        E = emin + xb ** penal * (1 - emin)
        U = solve(KE, edofMat, iK, jK, ndof, F, free, E)
        ce = np.einsum('ij,jk,ik->i', U[edofMat], KE, U[edofMat])
        c = (E * ce).sum()
        c0 = c0 or c
        dc = H @ ((-penal * (1 - emin) * xb ** (penal - 1) * ce * dxb) / Hs)
        dv = H @ (dxb / Hs)
        change = 1.0 if prev is None else np.abs(xb - prev).max()
        prev = xb
        x = mma.update(x, n * dc / c0, dv, 0.2, lambda y: proj(beta, H @ y / Hs)[0].mean() > volfrac)
        since += 1
        settled = change < 0.01
        if beta < 16 and (settled or since >= every):
            beta *= 2; since = 0
        elif settled and beta >= 16:
            break
    # final state at the final densities
    xb, _ = proj(beta, H @ x / Hs)
    E = emin + xb ** penal * (1 - emin)
    U = solve(KE, edofMat, iK, jK, ndof, F, free, E)
    return F @ U, it, np.mean(4 * xb * (1 - xb)), xb.mean()


if __name__ == "__main__":
    nelx, nely = 60, 20
    KE, edofMat, iK, jK, ndof, F, free, H, Hs = setup(nelx, nely, 2.4)
    csolid = F @ solve(KE, edofMat, iK, jK, ndof, F, free, np.ones(nelx * nely))
    print(f"solid compliance {csolid:.3f}")
    for ft in (1, 2, 3):
        c, it, grey = top88(nelx, nely, 0.5, 3.0, 2.4, ft)
        print(f"88-line ft{ft}: c {c:.2f} ({c / csolid:.3f} of solid), {it} iterations, greyness {grey:.3f}")
    for maxit in (60, 100, 300):
        c, it, grey, vol = ink_scheme(nelx, nely, 0.5, 2.4, maxit)
        print(f"ink's scheme, budget {maxit}: c {c:.3f} ({c / csolid:.4f} of solid), {it} iterations, greyness {grey:.4f}, volume {vol:.4f}")

#!/usr/bin/env python3
"""Benchmark 03: Euler-Bernoulli estimate for the blade (same section math as the .ink).

Integrates the loft field's cross-sections numerically (area, centroid, 2x2
inertia tensor per z), then the coupled bending deflection under the 50 N +y
tip load with the root clamped at z = 2. Shear and torsion are ignored."""
import math
SPAN, CR, TR, CT, TT, TW, E, F = 150, 60, 0.18, 30, 0.10, 30, 69000.0, 50.0
def naca(x, y, c, th):
    u = x / c + 0.25; uc = min(max(u, 0), 1)
    yt = 5*th*c*(0.2969*math.sqrt(uc) - 0.1260*uc - 0.3516*uc**2 + 0.2843*uc**3 - 0.1036*uc**4)
    return max(abs(y) - yt, -c*u, c*(u - 1))
def field(x, y, z):   # normalised morph, as blade-field in the .ink
    s = min(max(z / SPAN, 0), 1); a = math.radians(TW * s); c = CR + (CT - CR) * s
    xr, yr = math.cos(a)*x + math.sin(a)*y, math.cos(a)*y - math.sin(a)*x
    g0, g1 = naca(xr / c, yr / c, 1, TR), naca(xr / c, yr / c, 1, TT)
    return c * (g0 + (g1 - g0) * s)
def section(z, h=0.1):
    A = Sx = Sy = 0.0; pts = []
    nx, ny = int(65 / h), int(30 / h)
    for i in range(nx):
        x = -18 + (i + 0.5) * h
        for j in range(ny):
            y = -15 + (j + 0.5) * h
            if field(x, y, z) < 0: pts.append((x, y))
    A = len(pts) * h * h
    cx = sum(p[0] for p in pts) * h * h / A; cy = sum(p[1] for p in pts) * h * h / A
    Ixx = sum((p[1]-cy)**2 for p in pts) * h * h   # about x axis
    Iyy = sum((p[0]-cx)**2 for p in pts) * h * h
    Ixy = sum((p[0]-cx)*(p[1]-cy) for p in pts) * h * h
    return A, cx, cy, Ixx, Iyy, Ixy
zs = [2 + k * 4 for k in range(38)] + [150]
secs = {z: section(z) for z in zs}
# Moment of the +y tip force about the section: M = (-(L-z) F, 0) (about x).
# Unsymmetric bending: curvature k = (1/E) I^-1 M with I = [[Iyy, -Ixy], [-Ixy, Ixx]]
# in the (bending about y, bending about x) convention; deflections by double integration.
L = 150 - 1  # load centroid ~1 mm below the tip
def curv(z):
    A, cx, cy, Ixx, Iyy, Ixy = secs[z]
    Mx = -(L - z) * F
    det = Ixx * Iyy - Ixy ** 2
    # v'' (y deflection) and u'' (x deflection) for moment about x only:
    v2 = -Mx * Iyy / (E * det)
    u2 = -Mx * Ixy / (E * det)
    return u2, v2
# deflection at tip = integral (L - z) * curvature dz from the clamp (z = 2)
u = v = 0.0
for k in range(len(zs) - 1):
    z0, z1 = zs[k], zs[k + 1]; zm = 0.5 * (z0 + z1)
    c0, c1 = curv(z0), curv(z1)
    u += (L - zm) * 0.5 * (c0[0] + c1[0]) * (z1 - z0)
    v += (L - zm) * 0.5 * (c0[1] + c1[1]) * (z1 - z0)
for z in (2, 74, 150):
    A, cx, cy, Ixx, Iyy, Ixy = secs[z]
    print(f"z={z}: A={A:.1f} centroid=({cx:.2f},{cy:.2f}) Ixx={Ixx:.0f} Iyy={Iyy:.0f} Ixy={Ixy:.0f}")
print(f"beam tip deflection at the centroid: |u| = {math.hypot(u, v):.3f} mm, uy = {v:.3f} mm")
print(f"tip centroid at z=149: ({secs[150][1]:.2f}, {secs[150][2]:.2f})")
# Peak bending stress per section (unsymmetric bending, moment about x only).
def peak_stress(z, h=0.1):
    A, cx, cy, Ixx, Iyy, Ixy = section(z)
    Mx = (L - z) * F; det = Ixx * Iyy - Ixy ** 2
    best = 0.0
    for i in range(int(65 / h)):
        x = -18 + (i + 0.5) * h
        for j in range(int(30 / h)):
            y = -15 + (j + 0.5) * h
            if field(x, y, z) < 0:
                best = max(best, abs(Mx * (Iyy * (y - cy) - Ixy * (x - cx)) / det))
    return best
for z in (2, 50, 80, 102, 120, 140):
    print(f"z={z}: beam peak bending stress {peak_stress(z):.1f} MPa")

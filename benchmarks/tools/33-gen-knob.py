#!/usr/bin/env python3
"""Benchmark 33: write a large binary STL of a knurled, finned knob, and
print its exact volume.

The part (mm, axis z, z 0..60), a stack of prisms sharing one angular grid:
  - a through bore, a regular n-gon inscribed in r = 4;
  - z 0..20: a straight knurl, 120 triangular teeth with straight flanks,
    tips on r = 15.0 and roots on r = 14.4 (a 240-gon);
  - z 20..60: a hub, an n-gon inscribed in r = 8;
  - 19 fins, n-gons inscribed in r = 22, 1 mm thick at a 2 mm pitch
    (z 22..23, 24..25, ..., 58..59), so the gaps are 1 mm too.

Every cross-section is a polygon, so the mesh's volume is known exactly:
the sum of polygon areas times heights. The true-circle volume is printed
too (the n-gons differ from circles by 2e-6 relative at n = 2160).

usage: python3 benchmarks/tools/33-gen-knob.py OUT.stl [n]     (n: angular samples, a multiple of 240)
"""
import math, struct, sys

out = sys.argv[1]
n = int(sys.argv[2]) if len(sys.argv) > 2 else 2160
assert n % 240 == 0
TEETH, R_TIP, R_ROOT = 120, 15.0, 14.4
R_BORE, R_HUB, R_FIN = 4.0, 8.0, 22.0

def circle(r):
    return [(r * math.cos(2 * math.pi * i / n), r * math.sin(2 * math.pi * i / n)) for i in range(n)]

def knurl():
    # straight flanks between alternating tip and root vertices; the n samples
    # lie on the flanks (n / 240 per flank), so they share the angular grid's order
    corners = []
    for k in range(2 * TEETH):
        r = R_TIP if k % 2 == 0 else R_ROOT
        a = 2 * math.pi * k / (2 * TEETH)
        corners.append((r * math.cos(a), r * math.sin(a)))
    per = n // (2 * TEETH)
    pts = []
    for k in range(2 * TEETH):
        (x0, y0), (x1, y1) = corners[k], corners[(k + 1) % (2 * TEETH)]
        for j in range(per):
            s = j / per
            pts.append((x0 + s * (x1 - x0), y0 + s * (y1 - y0)))
    return pts, corners

def area(poly):
    return 0.5 * sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(poly, poly[1:] + poly[:1]))

bore, hub, fin = circle(R_BORE), circle(R_HUB), circle(R_FIN)
kn, kn_corners = knurl()
# segments (z0, z1, outer profile)
segs = [(0.0, 20.0, kn), (20.0, 22.0, hub)]
for k in range(19):
    z = 22.0 + 2 * k
    segs.append((z, z + 1, fin))
    segs.append((z + 1, z + 2, hub))

tris = []
def quad(a, b, c, d):        # counter-clockwise seen from outside: a b c d
    tris.append((a, b, c)); tris.append((a, c, d))

for z0, z1, prof in segs:
    for i in range(n):
        j = (i + 1) % n
        p, q = prof[i], prof[j]
        quad((p[0], p[1], z0), (q[0], q[1], z0), (q[0], q[1], z1), (p[0], p[1], z1))       # outer wall, normal out
        p, q = bore[i], bore[j]
        quad((p[0], p[1], z0), (p[0], p[1], z1), (q[0], q[1], z1), (q[0], q[1], z0))       # bore wall, normal in
# horizontal faces: at each level, the ring between the profile below and the profile above
levels = [(0.0, None, segs[0][2])] + [(segs[k][1], segs[k][2], segs[k + 1][2]) for k in range(len(segs) - 1)] + [(60.0, segs[-1][2], None)]
for z, below, above in levels:
    if below is None:                       # bottom face (normal -z): ring from bore to the profile
        for i in range(n):
            j = (i + 1) % n
            o0, o1, b0, b1 = above[i], above[j], bore[i], bore[j]
            quad((b0[0], b0[1], z), (b1[0], b1[1], z), (o1[0], o1[1], z), (o0[0], o0[1], z))
    elif above is None:                     # top face (normal +z)
        for i in range(n):
            j = (i + 1) % n
            o0, o1, b0, b1 = below[i], below[j], bore[i], bore[j]
            quad((b0[0], b0[1], z), (o0[0], o0[1], z), (o1[0], o1[1], z), (b1[0], b1[1], z))
    else:                                    # a step: the bigger profile's exposed ring
        ab = area(above) > area(below)
        big, small = (above, below) if ab else (below, above)
        for i in range(n):
            j = (i + 1) % n
            o0, o1, s0, s1 = big[i], big[j], small[i], small[j]
            if ab:   # ring faces down (underside of a fin / the grip's top is below)
                quad((s0[0], s0[1], z), (s1[0], s1[1], z), (o1[0], o1[1], z), (o0[0], o0[1], z))
            else:    # ring faces up
                quad((s0[0], s0[1], z), (o0[0], o0[1], z), (o1[0], o1[1], z), (s1[0], s1[1], z))

with open(out, "wb") as fh:
    fh.write(b"ink benchmark 33 knurled finned knob".ljust(80, b" "))
    fh.write(struct.pack("<I", len(tris)))
    for a, b, c in tris:
        ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
        vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        l = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
        fh.write(struct.pack("<12fH", nx / l, ny / l, nz / l, *a, *b, *c, 0))

v_mesh = sum((z1 - z0) * (area(prof) - area(bore)) for z0, z1, prof in segs)
# the true-circle part (the knurl is a true 240-gon either way)
a_kn = area(kn_corners)
v_true = 20 * a_kn + 2 * math.pi * R_HUB**2 + 19 * math.pi * (R_FIN**2 + R_HUB**2) - 60 * math.pi * R_BORE**2
print(f"triangles {len(tris)}  n {n}  file {out}")
print(f"knurl section area {a_kn:.6f} (240-gon), sampled-flank polygon {area(kn):.6f}")
print(f"volume: mesh exact {v_mesh:.4f}  true circles {v_true:.4f}  (rel {v_mesh / v_true - 1:.2e})")
# areas for reference
print(f"fins only: {19 * math.pi * (R_FIN**2 - R_HUB**2):.4f} of the true volume")

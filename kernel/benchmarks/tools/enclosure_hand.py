#!/usr/bin/env python3
"""Hand volumes for benchmarks/26-enclosure.ink, independent of the kernel.

Horizontal sections are integrated in z: each section of the drafted,
R-rounded cup is a rounded rectangle (core half-widths a, b grown by rho), so
A = 4ab + 4 rho (a + b) + pi rho^2. Bosses and lid notches (disks cut by
the rounded-rectangle cavity) are integrated on a fine 2D grid per slice.
"""
import math

L, W, H, wall, R, draft = 80, 60, 30, 2.0, 3.0, 1.0
boss_r, hole_r, bx, by = 3.5, 1.25, 35.5, 25.5
rib_t, rib_h = 1.2, 8.0
tl, lip_h, lip_t, gap, screw_r = 2.0, 3.0, 1.2, 0.2, 1.7
fillet = 0.5
ta = math.tan(math.radians(draft)); ca = math.cos(math.radians(draft))


def core(z):  # half-widths of the drafted core box at height z
    return L / 2 - R + (z - H) * ta, W / 2 - R + (z - H) * ta


def rrect_area(a, b, rho):
    return 4 * a * b + 4 * rho * (a + b) + math.pi * rho * rho


def rho_at(z, r):  # horizontal growth of a core offset by r (bottom face at z = R)
    if z >= R:
        return r / ca
    dz = R - z
    return math.sqrt(max(r * r - dz * dz, 0.0))


def integrate(f, z0, z1, n=3000):
    h = (z1 - z0) / n
    return sum(f(z0 + (i + 0.5) * h) for i in range(n)) * h


def in_rrect(px, py, a, b, rho):  # inside core (a, b) grown by rho (rho may be < 0)
    qx, qy = abs(px) - a, abs(py) - b
    d = math.hypot(max(qx, 0), max(qy, 0)) + min(max(qx, qy), 0)
    return d < rho


def disk_in_region(cx, cy, r, inside, h=0.02):
    n = int(2 * r / h) + 1
    s = 0
    for i in range(n):
        px = cx - r + (i + 0.5) * h
        for j in range(n):
            py = cy - r + (j + 0.5) * h
            if (px - cx) ** 2 + (py - cy) ** 2 < r * r and inside(px, py):
                s += 1
    return s * h * h


# ---- cup
v_out = integrate(lambda z: rrect_area(*core(max(z, R)), rho_at(z, R)), 0, H)
v_in = integrate(lambda z: rrect_area(*core(max(z, R)), rho_at(z, R - wall)), wall, H)
cup = v_out - v_in


# ---- ribs inside the cavity (floor z = wall to wall + rib_h; bottom fillet r 1 ignored)
def cav_half(z):
    a, b = core(z)
    return a + (R - wall) / ca, b + (R - wall) / ca


ribs = -6 * (1 - math.pi / 4) * (R - wall) ** 2 * rib_t   # rib ends inside the floor fillet

for k in range(400):
    z = wall + (k + 0.5) * rib_h / 400
    A, B = cav_half(z)
    ribs += (2 * 2 * B + 2 * A - 2 * rib_t) * rib_t * rib_h / 400   # 2 ribs along y, 1 along x minus crossings


# ---- bosses inside the cavity (per slice), minus pilot holes
bosses = 0.0
zs = [wall + (k + 0.5) * 0.1 for k in range(10)] + [R + (k + 0.5) * (H - R) / 27 for k in range(27)]
for z in zs:
    dz = 0.1 if z < R else (H - R) / 27
    r = boss_r + (H - z) * ta / ca
    a, b = core(max(z, R))
    rho = rho_at(z, R - wall)
    bosses += 4 * disk_in_region(bx, by, r, lambda px, py: in_rrect(px, py, a, b, rho)) * dz
holes = 4 * math.pi * hole_r ** 2 * (H - wall - 1)

base = cup + ribs + bosses - holes

# ---- lid
a, b = core(H)
plate = rrect_area(a, b, R / ca) * tl - 4 * math.pi * screw_r ** 2 * tl
a3, b3 = core(H - lip_h)
rho_o = (R - wall) / ca - gap
ring_area = rrect_area(a3, b3, rho_o) - 4 * (a3 + rho_o - lip_t) * (b3 + rho_o - lip_t)
c = boss_r + gap + 0.15            # lip stops square: removed where |x| > bx - c and |y| > by - c
h = 0.005
nx, ny = int((a3 + 1 - (bx - c)) / h), int((b3 + 1 - (by - c)) / h)
cnt = 0
for i in range(nx):
    px = bx - c + (i + 0.5) * h
    for j in range(ny):
        py = by - c + (j + 0.5) * h
        if in_rrect(px, py, a3, b3, rho_o) and not in_rrect(px, py, a3, b3, rho_o - lip_t):
            cnt += 1
notch = 4 * cnt * h * h
lid = plate + (ring_area - notch) * lip_h

print(f"cup {cup:.1f}  ribs {ribs:.1f}  bosses(in cavity) {bosses:.1f}  holes -{holes:.1f}")
print(f"base (no 0.5 fillets) {base:.1f}")
print(f"lid plate {plate:.1f}  lip ring {ring_area * lip_h:.1f}  notches -{notch * lip_h:.1f}  lid {lid:.1f}")

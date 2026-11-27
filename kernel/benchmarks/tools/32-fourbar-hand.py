"""Hand calculation for benchmark 32 (four-bar crank-rocker).

Ground A = (0,0), D = (d,0). Crank AB = a turns at 360 deg per unit t.
Coupler BC = b, rocker DC = c, open configuration (C above the x axis).
Prints Grashof class, rocker angle range (closed form at the two
crank-coupler collinear positions, and by dense sampling), the transmission
angle range, and the minimum distances that matter for collisions.
"""
import math, sys

a, b, c, d = 15.0, 60.0, 50.0, 55.0

def pose(th):
    bx, by = a * math.cos(th), a * math.sin(th)
    f = math.hypot(bx - d, by)                      # |BD|
    phi = math.pi + math.atan2(-by, d - bx)         # direction D -> B, away from the atan2 cut
    gam = math.acos((c * c + f * f - b * b) / (2 * c * f))
    th4 = phi - gam                                 # open branch: C above (A,D) line
    cx, cy = d + c * math.cos(th4), c * math.sin(th4)
    th3 = math.atan2(cy - by, cx - bx)
    mu = math.acos((b * b + c * c - f * f) / (2 * b * c))
    return bx, by, cx, cy, th3, th4, mu

s, l = min(a, b, c, d), max(a, b, c, d)
p_q = a + b + c + d - s - l
print(f"links a={a} b={b} c={c} d={d}: s+l={s+l} p+q={p_q} -> {'Grashof' if s+l <= p_q else 'non-Grashof'}; shortest is the {'crank' if s == a else 'other'} -> crank-rocker")

# rocker extremes: crank and coupler collinear, AC = b + a (extended) and b - a (folded)
for name, L in (("extended", b + a), ("folded", b - a)):
    psi = math.acos((d * d + c * c - L * L) / (2 * d * c))   # angle at D between DA and DC
    print(f"  {name}: AC = {L}: rocker angle = {180 - math.degrees(psi):.4f} deg")

# transmission angle extremes: |BD| = d - a (th = 0) and d + a (th = 180)
for name, f in (("min (crank at 0)", d - a), ("max (crank at 180)", d + a)):
    print(f"  transmission angle {name}: {math.degrees(math.acos((b*b + c*c - f*f) / (2*b*c))):.4f} deg")

N = 200000
th4s, th3s, mus = [], [], []
for i in range(N):
    th = 2 * math.pi * i / N
    bx, by, cx, cy, th3, th4, mu = pose(th)
    th4s.append(math.degrees(th4)); th3s.append(math.degrees(th3)); mus.append(math.degrees(mu))
i4lo, i4hi = th4s.index(min(th4s)), th4s.index(max(th4s))
print(f"sampled: rocker {min(th4s):.4f} (t={i4lo/N:.5f}) .. {max(th4s):.4f} (t={i4hi/N:.5f}) deg, swing {max(th4s)-min(th4s):.4f}")
print(f"sampled: coupler {min(th3s):.4f} .. {max(th3s):.4f} deg")
print(f"sampled: transmission {min(mus):.4f} .. {max(mus):.4f} deg")
fwd = (i4hi - i4lo) % N
print(f"time ratio (rocker forward/back stroke): {fwd/N:.5f} / {1-fwd/N:.5f}")
for t in (0, 0.25, 0.5, 0.75):
    bx, by, cx, cy, th3, th4, mu = pose(2 * math.pi * t)
    print(f"t={t}: B=({bx:.4f},{by:.4f}) C=({cx:.4f},{cy:.4f}) th3={math.degrees(th3):.4f} th4={math.degrees(th4):.4f} mu={math.degrees(mu):.4f}")

# segment-segment min distance in the plane between crank AB and rocker DC (centre lines)
def segdist(p, q, r, s_):
    def pt_seg(p, a_, b_):
        ax, ay = b_[0] - a_[0], b_[1] - a_[1]
        u = max(0, min(1, ((p[0]-a_[0])*ax + (p[1]-a_[1])*ay) / (ax*ax + ay*ay)))
        return math.hypot(p[0] - a_[0] - u*ax, p[1] - a_[1] - u*ay)
    return min(pt_seg(p, r, s_), pt_seg(q, r, s_), pt_seg(r, p, q), pt_seg(s_, p, q))
md = min(segdist((0, 0), pose(2*math.pi*i/2000)[:2], (d, 0), pose(2*math.pi*i/2000)[2:4]) for i in range(2000))
print(f"min centre-line distance crank AB to rocker DC over the turn: {md:.3f}")

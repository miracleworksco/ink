; Benchmark 09: bearing housing with press-fit claims.
;
; Task: an aluminium flange housing for a 6204 deep-groove ball bearing
; (d 20, D 47, B 14) on a steel shaft. The outer ring is pressed into the
; housing bore (N7-like: bore 46.980, 0.010 mm radial interference) and
; located by a shoulder; the inner ring is pressed on a k5 shaft (20.008,
; 0.004 mm radial interference). The fits are intended overlaps, declared
; as contact zones; interfere must find no other overlap.
;
; Reference values:
;   - outer fit: radial interference 0.010 mm (min 0.005, max 0.020 claimed);
;   - inner fit: radial interference 0.004 mm (min 0.002, max 0.006 claimed);
;   - shoulder bore 2 x 16 = 32 > inner ring OD 2 x 13.9: >= 2 mm clearance
;     to the inner ring (rotating part) is required;
;   - housing volume by hand (worked at the end of the file): 59 535.5 mm^3
;     -> 160.75 g at 2.70 g/cm^3.
;
; Pass criteria (ink assert --res 0.05, ink interfere, ink mass):
;   - all fit, shoulder, size and coaxiality claims proven; the mass claims
;     hold as estimates;
;   - interfere at --res 0.1: the only overlapping pairs are the two
;     declared press fits. The zones' surfaces touch the parts' over more
;     than the search's budget covers at the model's finest cells;
;   - housing mass within 0.5% of 160.75 g (mass-lo, mass-hi).

(import "std")

(param bore-r 23.490)       ; housing bore radius (bearing OD 47 / 2 = 23.5)
(param shaft-r 10.004)      ; shaft radius (bearing bore 20 / 2 = 10)
(param shoulder-r 16)       ; shoulder hole radius

; ---------------------------------------------------------------- bearing 6204
; Rings with ball grooves, 8 balls of 7.938 (5/16") on a 16.75 pitch radius.
(def ball-r 3.969)
(def pitch-r 16.75)
(def groove (torus-z pitch-r (+ ball-r 0.02)))
(def outer-ring (subtract (emax (max (tag ring-od (- (len2 x y) 23.5)) (- 19.3 (len2 x y))) (- (abs z) 7)) groove))
(def inner-ring (subtract (emax (max (- (len2 x y) 13.9) (tag ring-bore (- 10 (len2 x y)))) (- (abs z) 7)) groove))
(def balls (pattern-polar (move (sphere ball-r) pitch-r 0 0) 8))
(body bearing (shape (min outer-ring inner-ring balls)) (attr material "steel") (attr color "#B0B4BA"))

; ---------------------------------------------------------------- housing
; Axis z. Flange plate 80 x 80 x 8 with 4 x M8 clearance holes (9) on a
; 60 x 60 square; boss r 32 up to z = 22. Bore from z = 4 (shoulder) to 22.
(def plate (extrude-z (rounded-rect-2d 80 80 6) 0 8))
(def boss (move (cylinder-z 32 22) 0 0 11))
(def fixing-holes (symmetric-x (symmetric-y (move (cylinder-z 4.5 50) 30 30 0))))
(def bore (move (cylinder-z bore-r 18.001) 0 0 13.0005))       ; z 4 .. 22 (and past the top)
(def shoulder-hole (cylinder-z shoulder-r 100))
(def housing-field (subtract (min plate boss) (min fixing-holes bore shoulder-hole)))
(body housing (shape housing-field) (attr material "aluminum-6061") (attr density 2.7g/cm3) (attr color "#9DB4C8"))

; ---------------------------------------------------------------- shaft
(def shaft-field (move (cylinder-z shaft-r 60) 0 0 10))
(body shaft (shape shaft-field) (attr material "steel") (attr density 7.85g/cm3) (attr color "#707070"))

; The press fits' overlap is intended inside thin bands at the ring's
; outside (N7, 0.010 radial) and bore (k5, 0.004), in the bearing's coordinates.
(fn fit-band (r) (max (- (abs (- (len2 x y) r)) 0.02) (- (abs z) 7.1)))
(assembly root
  (instance housing housing)
  (instance bearing bearing (at (move 0 0 11)))          ; seated on the shoulder at z = 4
  (instance shaft shaft)
  (contact bearing housing (fit-band 23.5))
  (contact bearing shaft (fit-band 10)))

; ---------------------------------------------------------------- fit claims
; The fits as geometric proofs over the seat: 2-10 um bands, proven at
; --res 0.05.
; Outer ring seat, over the bearing width away from the chamfer-free edges.
(def seat-z0 4.5) (def seat-z1 17.5)
; min interference 0.005: housing material fills r in [23.495, 23.5]
(assert outer-fit-min (subset (extrude-z (ring-2d 23.495 23.5) seat-z0 seat-z1) (shape-of root/housing)) (accuracy 0.001))
; max interference 0.020: no housing material inside r = 23.48
(assert outer-fit-max (disjoint (extrude-z (circle-2d 23.48) seat-z0 seat-z1) (shape-of root/housing)) (accuracy 0.001))
; inner fit: the shaft fills r in [10, 10.002] (>= 0.002) and stops before 10.006
(assert inner-fit-min (subset (extrude-z (ring-2d 10 10.002) seat-z0 seat-z1) (shape-of root/shaft)) (accuracy 0.0005))
(assert inner-fit-max (disjoint (extrude-z (ring-2d 10.006 12) seat-z0 seat-z1) (shape-of root/shaft)) (accuracy 0.001))
; the bearing sits on the shoulder: outer ring touches z = 4 but does not cut it
; The bearing seats on the shoulder: no search separates the faces at the
; model's floor within its budget, so it is decided at cells of 0.05 (to 0.1
; mm).
(assert seated-on-shoulder (disjoint (shape-of root/bearing) (max (shape-of root/housing) (- z 4))) (res 0.05))
; the shoulder clears the rotating inner ring (true gap 2.1 mm; 2.0 proves at
; res 0.05, 2.2 is disproven with a witness at 2.109). The ring is taken from
; the placed bearing by name, so it follows the assembly.
(def inner-ring-world (shape-of root/bearing inner-ring))
(assert shoulder-clears-inner (clearance (shape-of root/housing) inner-ring-world 2.0))
; bore and shaft radii within 0.0005 of nominal, by point probes at four
; azimuths at both ends of the seat: so they are coaxial within 0.001
(def probe-at (vec (vec 1 0 5) (vec 0 1 5) (vec -1 0 5) (vec 0 -1 5)
                   (vec 1 0 17) (vec 0 1 17) (vec -1 0 17) (vec 0 -1 17)))
(fn ring-probe (r i) (let ((d (nth probe-at i))) (probe (* r (nth d 0)) (* r (nth d 1)) (nth d 2))))
(for i (range 8) (assert bore-open-$i (disjoint (ring-probe 23.4895 i) (shape-of root/housing))))
(for i (range 8) (assert bore-metal-$i (subset (ring-probe 23.4905 i) (shape-of root/housing))))
(for i (range 8) (assert shaft-metal-$i (subset (ring-probe 10.0035 i) (shape-of root/shaft))))
(for i (range 8) (assert shaft-open-$i (disjoint (ring-probe 10.0045 i) (shape-of root/shaft))))

; Housing volume by hand (a solid of revolution plus the plate):
;   plate outside r 32: (80^2 - (4 - pi) 6^2) 8 - pi 32^2 8   = 25216.9
;   r < 32 over z 0..22: pi 32^2 22                          = 70774.4
;   minus bore pi 23.49^2 18, shoulder pi 16^2 4, holes 4 pi 4.5^2 8
;   = 25216.9 + 70774.4 - 31202.5 - 3217.0 - 2035.8           = 59535.5 mm^3
;   mass at 2.70 g/cm^3: 160.75 g
; within 0.5% (measured):
(assert mass-lo (integral housing (attr density) > 159.95))
(assert mass-hi (integral housing (attr density) < 161.55))

; Benchmark 32: pin-jointed four-bar linkage (crank-rocker) over a full turn.
;
; Task: crank a = 15 (ground pivot A at the origin) turning once per unit t,
; coupler b = 60, rocker c = 50 (ground pivot D at (55, 0)), ground d = 55.
; Real pins in holes at all four joints: pins r 3, holes r 3.1 (0.1 radial
; clearance), links 12 wide. Layers in z: ground plate 0..5 with pins A and D
; up to 10.5; crank and rocker 6..10 with pins B and C up to 15.5; coupler
; 11..15. Kinematics in closed form, checked over the motion (t-span 0 1)
; with t ranged; tools/32-solved.ink checks them against (solve) at sampled t.
;
; Hand calculation (tools/32-fourbar-hand.py):
;   - Grashof: s + l = 15 + 60 = 75 < p + q = 50 + 55 = 105, the crank is the
;     shortest link: crank-rocker.
;   - Rocker extremes at the crank-coupler collinear positions, from the
;     triangle A D C with AC = b + a = 75 and b - a = 45:
;     cos psi = (d^2 + c^2 - AC^2) / (2 d c) = -0.01818 and 0.63636,
;     rocker angle 180 - psi = 88.958 deg (t = 0.1161) and 129.521 deg
;     (t = 0.6639): swing 40.563 deg, time ratio 0.548 / 0.452.
;   - Transmission angle mu (coupler to rocker): cos mu = (b^2 + c^2 - BD^2)
;     / (2 b c) with BD = d -+ a: min 41.410 deg (t = 0), max 78.463 deg
;     (t = 0.5). It equals rocker angle minus coupler angle.
;   - Links never cross: crank and rocker centre lines stay 39.7 apart (> 12).
;
; Pass criteria (ink assert over the file's t-span, one claim per call):
;   - proven: grashof, rocker-range-lo/hi (88.9 / 129.6), transmission-min
;     (> 40), transmission-tight (> 41.3), transmission-from-placements
;     (> 41.3), loop-closes-x/y, crank-clears-rocker, coupler-clears-ground;
;   - holds numerically (measured at steps: subset and overlaps claims, and
;     shapes turned by rotate-z of t, not placed by an instance):
;     pin-a-seated, pin-d-seated, pin-c-engaged, rel-pin-b-seated,
;     rel-pin-c-seated, rel-joint-b-disjoint, rel-joint-c-disjoint;
;   - disproven: the three value controls, control-joint-c-clearance-disproven,
;     control-rel-joint-c-clearance-disproven;
;   - tools/32-solved.ink at t = 0, 0.11611, 0.25, 0.5, 0.66387, 0.75:
;     solved and closed-form angles agree within 1e-6 deg: holds, an
;     estimate (accuracy null: solved values are numerical).

(import "std")
(t-span 0 1)

; ---------------------------------------------------------------- linkage
(param a 15)          ; crank A-B (A: ground pivot at the origin)
(param b 60)          ; coupler B-C
(param c 50)          ; rocker D-C (D: ground pivot at (d, 0))
(param d 55)          ; ground A-D
(param pin-r 3)       ; pins
(param hole-r 3.1)    ; holes: 0.1 radial clearance
(param w 12)          ; link width

(def r2d (/ 180 pi))
(def th (* 2 pi t))                                   ; crank angle, one turn per unit t
(def bx (* a (cos th)))
(def by (* a (sin th)))
(def f (sqrt (- (+ (square a) (square d)) (* 2 a d (cos th)))))      ; |BD|, 40 .. 70
; direction D -> B, written as pi + atan2 of the reversed vector, whose x
; (d - bx >= 40) stays away from atan2's branch cut
(def phi (+ pi (atan2 (neg by) (- d bx))))
(def gam (acos (/ (- (+ (square c) (square f)) (square b)) (* 2 c f))))
(def th4 (- phi gam))                                 ; rocker angle, open branch (C above AD)
(def cx (+ d (* c (cos th4))))
(def cy (* c (sin th4)))
(def th3 (atan2 (- cy by) (- cx bx)))                 ; coupler angle
(def crank-deg (* 360 t))
(def coupler-deg (* r2d th3))
(def rocker-deg (* r2d th4))
; transmission angle (coupler to rocker), law of cosines in triangle BCD
(def mu-deg (* r2d (acos (/ (- (+ (square b) (square c)) (square f)) (* 2 b c)))))

; ---------------------------------------------------------------- parts (own frames)
; A link of length l along +x from its own origin, z0..z1, with rounded ends.
(fn link (l z0 z1) (extrude-z (move (rounded-rect-2d (+ l w) w (* 0.5 w)) (* 0.5 l) 0 0) z0 z1))
(fn pin-at (px z0 z1) (move (cylinder-z pin-r (- z1 z0)) px 0 (* 0.5 (+ z0 z1))))
(fn hole-at (px) (move (cylinder-z hole-r 100) px 0 0))

; ground: plate z 0..5 under both pivots, pins up to z 10.5 (through the crank and rocker eyes)
(def plate (extrude-z (move (rounded-rect-2d (+ d 24) 24 12) (* 0.5 d) 0 0) 0 5))
(def pin-a (pin-at 0 0 10.5))
(def pin-d (pin-at d 0 10.5))
(body ground (shape (min plate pin-a pin-d)) (attr color "#999999")
  (frame pivot-a (at)) (frame pivot-d (at (move d 0 0))))

; crank: layer z 6..10, eye at A, pin at B up through the coupler (z 6..15.5)
(def crank-hole (hole-at 0))
(def crank-pin (pin-at a 6 15.5))
(def crank-shape (min (subtract (link a 6 10) crank-hole) crank-pin))
(body crank (shape crank-shape) (attr color "#D79B00") (frame pin-b (at (move a 0 0))))

; rocker: layer z 6..10, eye at D, pin at C up through the coupler
(def rocker-hole (hole-at 0))
(def rocker-pin (pin-at c 6 15.5))
(def rocker-shape (min (subtract (link c 6 10) rocker-hole) rocker-pin))
(body rocker (shape rocker-shape) (attr color "#82B366"))

; coupler: layer z 11..15, eyes at B (own origin) and C (b, 0)
(def coupler-hole-b (hole-at 0))
(def coupler-hole-c (hole-at b))
(def coupler-shape (subtract (link b 11 15) (min coupler-hole-b coupler-hole-c)))
(body coupler (shape coupler-shape) (attr color "#6C8EBF"))

(assembly root
  (instance ground ground)
  (instance crank crank (at (rot-z crank-deg) (frame-of ground pivot-a)))
  (instance coupler coupler (at (rot-z (- coupler-deg crank-deg)) (frame-of crank pin-b)))
  (instance rocker rocker (at (rot-z rocker-deg) (frame-of ground pivot-d))))

; ---------------------------------------------------------------- claims (run with --t-span 0:1)
; linkage class (params only)
(assert grashof (value (- (+ a b) (+ c d)) < 0))
; output (rocker) angle range against the hand calculation 88.958 .. 129.521
(assert rocker-range-lo (value rocker-deg > 88.9))
(assert rocker-range-hi (value rocker-deg < 129.6))
; expect-status: control-rocker-reaches-lo-disproven fails
(assert control-rocker-reaches-lo-disproven (value rocker-deg > 89.05))
; expect-status: control-rocker-reaches-hi-disproven fails
(assert control-rocker-reaches-hi-disproven (value rocker-deg < 129.4))
; transmission angle (hand: 41.410 .. 78.463)
(assert transmission-min (value mu-deg > 40))
(assert transmission-tight (value mu-deg > 41.3))
; expect-status: control-transmission-disproven fails
(assert control-transmission-disproven (value mu-deg > 41.5))
; the same angle read off the placements (rocker minus coupler): the interval
; dependency between the two angles is the test here
(assert transmission-from-placements (value (- rocker-deg coupler-deg) > 41.3))
; loop closure: the coupler's C eye and the rocker's C pin coincide
(assert loop-closes-x (value (abs (- (+ bx (* b (cos th3))) cx)) < 0.001))
(assert loop-closes-y (value (abs (- (+ by (* b (sin th3))) cy)) < 0.001))
; no collision between the links that never touch
(assert crank-clears-rocker (disjoint (shape-of root/crank) (shape-of root/rocker)))
(assert coupler-clears-ground (disjoint (shape-of root/coupler) (shape-of root/ground)))
; pins seated: each pin's length inside the eye stays inside the hole
(assert pin-a-seated (subset (max (shape-of root/ground pin-a) (slab-z 6 10)) (shape-of root/crank crank-hole)))
(assert pin-d-seated (subset (max (shape-of root/ground pin-d) (slab-z 6 10)) (shape-of root/rocker rocker-hole)))
(assert pin-b-seated (subset (max (shape-of root/crank crank-pin) (slab-z 11 15)) (shape-of root/coupler coupler-hole-b)))
(assert pin-c-seated (subset (max (shape-of root/rocker rocker-pin) (slab-z 11 15)) (shape-of root/coupler coupler-hole-c)))
; and engaged (the pin is in the hole at every moment, not beside it)
(assert pin-c-engaged (overlaps (shape-of root/rocker rocker-pin) (max (shape-of root/coupler coupler-hole-c) (slab-z 11 15))))
; Joint running clearance: 0.1 radial; test 0.05 and a failing 0.11 control.
; Near-gap motion queries remain slow.
(assert joint-c-clearance (clearance (shape-of root/coupler) (shape-of root/rocker) 0.05))
; expect-status: control-joint-c-clearance-disproven fails
(assert control-joint-c-clearance-disproven (clearance (shape-of root/coupler) (shape-of root/rocker) 0.11))

; ---------------------------------------------------------------- the joints in the coupler's frame
; World-frame joint claims at B and C must bound two parts that both move, and
; interval bounds over a time piece do not know they move together. In the
; coupler's frame (B at its origin, C at (b, 0)) the crank and the rocker only
; turn about the pins, by the relative angles below. Exact given loop closure
; (proven above): crank q = R(th3 - th) p + (a, 0); rocker q = R(th3 - th4) (p - (b, 0)) + (c, 0).
(def crank-in-coupler (rotate-z (move crank-shape (neg a) 0 0) (- crank-deg coupler-deg)))
(def rocker-in-coupler (move (rotate-z (move rocker-shape (neg c) 0 0) (- rocker-deg coupler-deg)) b 0 0))
(def crank-pin-in-coupler (rotate-z (move crank-pin (neg a) 0 0) (- crank-deg coupler-deg)))
(def rocker-pin-in-coupler (move (rotate-z (move rocker-pin (neg c) 0 0) (- rocker-deg coupler-deg)) b 0 0))
(assert rel-pin-b-seated (subset (max crank-pin-in-coupler (slab-z 11 15)) coupler-hole-b))
(assert rel-pin-c-seated (subset (max rocker-pin-in-coupler (slab-z 11 15)) coupler-hole-c))
(assert rel-joint-b-disjoint (disjoint coupler-shape crank-in-coupler))
(assert rel-joint-c-disjoint (disjoint coupler-shape rocker-in-coupler))
(assert rel-joint-c-clearance (clearance coupler-shape rocker-in-coupler 0.09))
; expect-status: control-rel-joint-c-clearance-disproven fails
(assert control-rel-joint-c-clearance-disproven (clearance coupler-shape rocker-in-coupler 0.11))

; A CoreXY's moving parts, reduced to what its dynamics rests on: a head
; (x) on a gantry (y), and two steppers turning (x + y) / R and (x - y) / R
; (R = 20 / pi mm: a 20-tooth GT2 pulley's pitch radius), each with a
; pulley and a grub screw in its pulley's radial hole riding its shaft, its
; axis along gravity. Each belt is open, both ends clamped to the head: from
; the head along the gantry, round a wheel on the gantry's end, up the side
; to its motor, across the top, down the far side round a foot idler, back
; up round a second gantry wheel and along the gantry to the head. Its route
; moves with x and y (the gantry's wheels and the head carry it), its
; material anchored at the head's clamp. B is A mirrored in x, its top 50 mm
; lower (100 mm shorter); each belt's six spans are couplings (15 lbf
; installed, efficiency 0.95). Idlers have spin-inertia 0.
(import "std")
(import "stepper.ink" as motor)
(def R (/ 20 pi))
(def W 350)
(def F 350)
(def h 20)
(def HA 400)
(def HB 350)
; A point of a route, its x mirrored by s.
(fn at-s (s a b) (vec (* s a) b 0))
; Belt A's route (s = 1) or B's (s = -1, A's mirrored in x) in its plane,
; top H, the head at (x, y): its clamps at x +- h on the line y - R.
(fn route (s x y H)
  (let ((u (* s x)) (k (vec 0 0 s)))
    (path (vec (path-line (at-s s (+ u h) (- y R)) (at-s s W (- y R)))
               (path-arc (at-s s W y) k R (vec 0 -1 0) (/ pi 2))
               (path-line (at-s s (+ W R) y) (at-s s (+ W R) H))
               (path-arc (at-s s W H) k R (at-s s 1 0) (/ pi 2))
               (path-line (at-s s W (+ H R)) (at-s s (neg W) (+ H R)))
               (path-arc (at-s s (neg W) H) k R (vec 0 1 0) (/ pi 2))
               (path-line (at-s s (neg (+ W R)) H) (at-s s (neg (+ W R)) (neg F)))
               (path-arc (at-s s (neg W) (neg F)) k R (at-s s -1 0) pi)
               (path-line (at-s s (- R W) (neg F)) (at-s s (- R W) (- y (* 2 R))))
               (path-arc (at-s s (- (* 2 R) W) (- y (* 2 R))) (vec 0 0 (neg s)) R (at-s s -1 0) (/ pi 2))
               (path-line (at-s s (- (* 2 R) W) (- y R)) (at-s s (- u h) (- y R))))
          (vec 0 0 1))))
; Its length, the same at every x, y.
(fn route-length (H) (+ (* 4 W) (* 2 H) (* 2 F) (* -2 h) (* -4 R) (* 3 pi R)))
(body (belt (s 1) (H 400) (u 0) (v 0)) (shape (path-sweep (max (- (abs x) 0.7) (- (abs y) 3)) (route s u v H))) (attr density 1.2g/cm3)
  (attr carried (carried (max (- (abs y) 0.7) (- (abs z) 3)) (route s u v H) (route-length H) 0)))
(body gantry (shape (move (box 700 40 20) 0 0 -10)) (attr density 2.7g/cm3))
(body head (shape (move (box 50 50 30) 0 0 15)) (attr density 1g/cm3))
(body pulley (shape (subtract (cylinder-z (- R 0.4) 6) (cylinder-z 2.5 7) (move (box 4 3 3) 4.5 0 0))) (attr density 2.7g/cm3))
(body grub (shape (move (box 3 3 3) 5 0 0)) (attr density 7.85g/cm3))
(body idler (shape (cylinder-z (- R 0.4) 6)) (attr density 2.7g/cm3) (attr spin-inertia 0))
(fn turned (s) (rad->deg (/ s R)))
; A span of a belt in plane z, from a point p of its route to q.
(fn span (p q z) (vec (vec (nth p 0) (nth p 1) z) (vec (nth q 0) (nth q 1) z)))
(assembly (corexy (x 0) (y 0))
  (instance gantry gantry (at) (slide y y))
  (instance head head (frame-of gantry) (slide x x))
  ; A, in the plane z = 21: its motor top right.
  (instance motorA (motor.drive (turned (+ x y))) (at (move W HA 0)))
  (instance pulleyA pulley (at (move 0 0 -3) (frame-of motorA shaft)))
  (instance grubA grub (frame-of pulleyA))
  (instance cornerA idler (at (move (neg W) HA 21)) (turn z (turned (+ x y))))
  (instance footA idler (at (move (neg W) (neg F) 21)) (turn z (turned (+ x y))))
  (instance outA idler (at (move W 0 21) (frame-of gantry)) (turn z (turned x)))
  (instance backA idler (at (move (- (* 2 R) W) (* -2 R) 21) (frame-of gantry)) (turn z (turned (neg x))))
  (instance beltA (belt 1 HA x y) (at (move 0 0 21)))
  (relate head outA (belt-span (span (at-s 1 (+ x h) (- y R)) (at-s 1 W (- y R)) 21) 15lbf 0.95))
  (relate outA pulleyA (belt-span (span (at-s 1 (+ W R) y) (at-s 1 (+ W R) HA) 21) 15lbf 0.95))
  (relate pulleyA cornerA (belt-span (span (at-s 1 W (+ HA R)) (at-s 1 (neg W) (+ HA R)) 21) 15lbf 0.95))
  (relate cornerA footA (belt-span (span (at-s 1 (neg (+ W R)) HA) (at-s 1 (neg (+ W R)) (neg F)) 21) 15lbf 0.95))
  (relate footA backA (belt-span (span (at-s 1 (- R W) (neg F)) (at-s 1 (- R W) (- y (* 2 R))) 21) 15lbf 0.95))
  (relate backA head (belt-span (span (at-s 1 (- (* 2 R) W) (- y R)) (at-s 1 (- x h) (- y R)) 21) 15lbf 0.95))
  ; B, in the plane z = 41: A mirrored, its motor top left.
  (instance motorB (motor.drive (turned (- x y))) (at (move (neg W) HB 20)))
  (instance pulleyB pulley (at (move 0 0 -3) (frame-of motorB shaft)))
  (instance grubB grub (frame-of pulleyB))
  (instance cornerB idler (at (move W HB 41)) (turn z (turned (- x y))))
  (instance footB idler (at (move W (neg F) 41)) (turn z (turned (- x y))))
  (instance outB idler (at (move (neg W) 0 41) (frame-of gantry)) (turn z (turned x)))
  (instance backB idler (at (move (- W (* 2 R)) (* -2 R) 41) (frame-of gantry)) (turn z (turned (neg x))))
  (instance beltB (belt -1 HB x y) (at (move 0 0 41)))
  (relate head outB (belt-span (span (at-s -1 (- h x) (- y R)) (at-s -1 W (- y R)) 41) 15lbf 0.95))
  (relate outB pulleyB (belt-span (span (at-s -1 (+ W R) y) (at-s -1 (+ W R) HB) 41) 15lbf 0.95))
  (relate pulleyB cornerB (belt-span (span (at-s -1 W (+ HB R)) (at-s -1 (neg W) (+ HB R)) 41) 15lbf 0.95))
  (relate cornerB footB (belt-span (span (at-s -1 (neg (+ W R)) HB) (at-s -1 (neg (+ W R)) (neg F)) 41) 15lbf 0.95))
  (relate footB backB (belt-span (span (at-s -1 (- R W) (neg F)) (at-s -1 (- R W) (- y (* 2 R))) 41) 15lbf 0.95))
  (relate backB head (belt-span (span (at-s -1 (- (* 2 R) W) (- y R)) (at-s -1 (- (neg x) h) (- y R)) 41) 15lbf 0.95)))

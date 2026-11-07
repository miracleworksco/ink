; A four-bar (crank 15, coupler 60, rocker 50, ground 55) over one turn,
; its angles in closed form. The rocker's pin (r 3) runs in the coupler's
; far eye (r 3.1), but the eye is placed through the coupler's angle (an
; atan2 of the rocker's) and the pin through the rocker's: they meet only
; by the loop's closure, which no bound of either placement over a piece
; of the motion sees. Seen from the rocker the coupler turns about the
; pin they share, its drift off that pivot bounded over the motion
; (motion::pivot) and charged to each claim: the pin is seated in the eye,
; engaged with it and clear of the coupler by 0.05, at every moment. A
; coupler whose eye is 0.3 too far out (the loop not closed) holds the
; pin no longer: it fails.
(import "std")
(t-span 0 1)
(param a 15)
(param b 60)
(param c 50)
(param d 55)
(def th (* 2 pi t))
(def bx (* a (cos th)))
(def by (* a (sin th)))
(def f (sqrt (- (+ (square a) (square d)) (* 2 a d (cos th)))))
(def phi (+ pi (atan2 (neg by) (- d bx))))
(def th4 (- phi (acos (/ (- (+ (square c) (square f)) (square b)) (* 2 c f)))))
(def th3 (atan2 (- (* c (sin th4)) by) (- (+ d (* c (cos th4))) bx)))
(def r2d (/ 180 pi))
(fn link (l z0 z1) (extrude-z (move (rounded-rect-2d (+ l 12) 12 6) (* 0.5 l) 0 0) z0 z1))
(def rocker-pin (move (cylinder-z 3 9.5) c 0 10.75))
(body rocker (shape (min (link c 6 10) rocker-pin)))
(def eye (move (cylinder-z 3.1 100) b 0 0))
(body coupler (shape (subtract (link b 11 15) eye)))
(def far-eye (move (cylinder-z 3.1 100) (+ b 0.3) 0 0))
(body long (shape (subtract (link (+ b 0.3) 11 15) far-eye)))
(body crank (shape (link a 6 10)) (frame pin-b (at (move a 0 0))))
(assembly root
  (instance crank crank (at (rot-z (* 360 t))))
  (instance coupler coupler (at (rot-z (- (* r2d th3) (* 360 t))) (frame-of crank pin-b)))
  (instance long long (at (rot-z (- (* r2d th3) (* 360 t))) (frame-of crank pin-b)))
  (instance rocker rocker (at (rot-z (* r2d th4)) (move d 0 0))))
(assert seated (subset (max (shape-of root/rocker rocker-pin) (slab-z 11 15)) (shape-of root/coupler eye)))
(assert engaged (overlaps (shape-of root/rocker rocker-pin) (max (shape-of root/coupler eye) (slab-z 11 15))))
(assert clear (clearance (shape-of root/coupler) (shape-of root/rocker) 0.05))
(assert long-seated (subset (max (shape-of root/rocker rocker-pin) (slab-z 11 15)) (shape-of root/long far-eye)))
; expect-status: seated holds
; expect-status: engaged holds
; expect-status: clear holds
; expect-status: long-seated fails

; Benchmark 10: slider-crank mechanism with motion proofs over a full turn.
;
; Task: crank r = 20 turning once per unit time about the origin, connecting
; rod l = 70, slider running in a two-rail guide along +x. Pins are 6 mm with
; 6.4 mm holes (0.2 radial clearance). Parts stack in z: crank 0..4, rod
; 5..9, slider 10..18, rails 10..18. Kinematics are closed form (no solve):
;   B = r (cos th, sin th),  x_C = r cos th + sqrt(l^2 - r^2 sin^2 th)
; checked over the motion (t-span 0 1) with t ranged (docs/IR.md Motion).
;
; Reference values:
;   - stroke = 2 r = 40: slider pin from x = 50 (t = 0.5) to 90 (t = 0);
;   - max rod angle = asin(r / l) = 16.60 deg;
;   - slider block (20 long) stays within x 40 .. 100, inside rails 30 .. 115.
;
; Pass criteria (ink assert, ink interfere: over the file's t-span):
;   - rod angle within +-16.7 deg for all t (proven);
;   - slider stays inside the guide envelope for all t (holds numerically:
;     a within claim is measured at steps over a motion);
;   - the rod never touches the rails for all t (proven);
;   - the slider clears the rails by 0.2 above z 9.9 (holds numerically: the
;     slider cut by a half-space is not a union of placed bodies);
;   - interfere over t in [0, 1]: no overlapping pair, none unknown
;     (pins in holes have 0.2 mm clearance, so everything is disjoint);
;   - slider pin at t = 0 and t = 0.5: x = 90 and 50 (eval).

(import "std")
(t-span 0 1)

(param r 20)
(param l 70)
(def th (* 2 pi t))
(def bx (* r (cos th)))
(def by (* r (sin th)))
(def cx (+ bx (sqrt (- (square l) (square by)))))
(def rod-deg (* (/ 180 pi) (atan2 (neg by) (- cx bx))))
(def crank-deg (* 360 t))

; ---------------------------------------------------------------- parts (own frames)
(fn link-2d (len w) (rounded-rect-2d (+ len w) w (* 0.5 w)))
(def pin-hole (- 3.2 (len2 x y)))
; crank: bar from the origin to (r, 0), boss at the origin, pin at the far end
(def crank-bar (extrude-z (move (link-2d r 12) (* 0.5 r) 0 0) 0 4))
(def crank-pin (move (cylinder-z 3 5.5) r 0 6.75))                 ; z 4 .. 9.5
(body crank (shape (min crank-bar crank-pin)) (attr color "#D79B00") (frame pin (at (move r 0 0))))

; rod: from (0,0) to (l,0) with two pin holes
(def rod-bar (extrude-z (move (link-2d l 10) (* 0.5 l) 0 0) 5 9))
(body rod (shape (max rod-bar pin-hole (move pin-hole l 0 0))) (attr color "#6C8EBF"))

; slider: block 20 x 16 x 8 at z 10..18 centred on its pin, pin down to z 4.5
(def slider-block (move (box 20 16 8) 0 0 14))
(def slider-pin (move (cylinder-z 3 5.6) 0 0 7.3))                ; z 4.5 .. 10.1
(body slider (shape (min slider-block slider-pin)) (attr color "#82B366"))

; frame: two rails (|y| 8.5 .. 14, z 10 .. 18, x 30 .. 115)
(def rails (max (- (abs (- (abs y) 11.25)) 2.75) (slab-z 10 18) (- (abs (- x 72.5)) 42.5)))
; its pivot (the crank's bearing, not modelled) and its guide
(body frame (shape rails) (attr color "#999999") (frame pivot (at)) (frame guide (at)))

(assembly root
  (instance frame frame)
  (instance crank crank (at (rot-z crank-deg) (frame-of frame pivot)))
  (instance rod rod (at (rot-z (- rod-deg crank-deg)) (frame-of crank pin)))
  (instance slider slider (at (move cx 0 0) (frame-of frame guide))))

; ---------------------------------------------------------------- motion claims (all t)
(assert rod-angle (value (abs rod-deg) < 17.5))
(assert rod-angle-tight (value (abs rod-deg) < 16.7))
(assert slider-in-guide (within (shape-of root/slider) (39.9 -8.1 4.4) (100.1 8.1 18.1)))
(assert rod-clears-rails (disjoint (shape-of root/rod) (shape-of root/frame)))
(assert slider-clears-rails (clearance (shape-of root/frame) (max (shape-of root/slider) (- 9.9 z)) 0.2))

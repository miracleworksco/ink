; A crank 15 long turns once about the origin over t 0..1; its pin (r 3, at
; x = 15) carries a follower link by its eye (r 3.1: 0.1 radial clearance),
; which rocks back and forth about the pin as it goes round. In the world
; the pin and the eye both sweep a circle at 94 per unit t, so a bound over
; a piece of the motion blurs each far past the 0.1 gap; in the crank's
; frame (motion::moving_frame) the pin rests and the eye only turns about
; it. So the pin's length inside the eye is proven inside it over every
; moment, and a pin too fat for the eye (r 3.4) fails, its witness a world
; point at the joint.
(import "std")
(t-span 0 1)
(def crank-pin (move (cylinder-z 3 9.5) 15 0 10.75))
(body crank (shape (min (extrude-z (move (rounded-rect-2d 27 12 6) 7.5 0 0) 6 10) crank-pin))
  (frame pin (at (move 15 0 0))))
(def eye-hole (cylinder-z 3.1 100))
(body follower (shape (subtract (extrude-z (move (rounded-rect-2d 52 12 6) 20 0 0) 11 15) eye-hole)))
(body fat (shape (move (cylinder-z 3.4 3) 15 0 13)))
(assembly root
  (instance crank crank (at (rot-z (* 360 t))))
  (instance follower follower (at (rot-z (* 30 (sin (* 4 pi t)))) (frame-of crank pin)))
  (instance fat fat (at (rot-z (* 360 t)))))
(assert seated (subset (max (shape-of root/crank crank-pin) (slab-z 11 15)) (shape-of root/follower eye-hole)))
(assert clear (disjoint (max (shape-of root/crank crank-pin) (slab-z 11 15)) (shape-of root/follower)))
(assert too-fat (subset (max (shape-of root/fat) (slab-z 11 15)) (shape-of root/follower eye-hole)))
; expect-status: seated holds
; expect-status: clear holds
; expect-status: too-fat fails

; A crank driven from its slider through a solve, from dead centre (50 cos
; th = 50 + s, s = -20 t): its speed has no bound at t = 0, where its rate
; is not defined, and passes 1000 deg/s for t < 6.6e-4. Measured at the
; steps' middles alone it held (794.7 deg/s); checked at the span's ends
; too, it is unknown there.
(import "std")
(t-span 0 1)
(param th 0.3)
(body block (shape (box 10 10 10)))
(body bar (shape (move (box 50 4 4) 25 0 0)))
(assembly (rig (s 0))
  (instance slider block (at (move 50 0 0)) (slide x s))
  (instance crank bar (at) (turn z (rad->deg th)) (speed-limits -1000 1000)))
(assembly root (instance r (rig (* -20 t))))
(solve loop (unknowns th) (== (* 50 (cos th)) (nth (get (frame-of root/r/slider) "origin") 0)))
; expect-status: speed-limits:rig/crank unknown
; expect-status: speed-limits:root/r/crank unknown

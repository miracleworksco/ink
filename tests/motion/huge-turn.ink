; A turn by 1e20 degrees a step: 280 degrees past a whole turn, so the turn
; really jumps at t = 1, and its speed limit is unknown there (whole turns
; are counted in integers, not by 1e20 / 360 rounded to a whole number).
; expect-status: speed-limits:root/arm unknown
(import "std")
(t-span 0.9 1.1)
(body arm (shape (box 10 2 2)))
(assembly root
  (instance arm arm (at) (turn z (* 1e20 (floor t))) (speed-limits -1 1)))

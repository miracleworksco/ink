; Half of (a + b - a), a 2^53 whole turns a step and b one: the turn jumps
; half round at t = 1, so its speed limit is unknown there (whole turns are
; summed in integers; in floats 2^53 + 1 is 2^53, and the turn was lost).
; expect-status: speed-limits:root/arm unknown
(import "std")
(t-span 0.9 1.1)
(def n (floor t))
(def a (deg->rad (* 3242591731706757120 n)))
(def b (deg->rad (* 360 n)))
(def angle (rad->deg (* 0.5 (- (+ a b) a))))
(body arm (shape (box 10 2 2)))
(assembly root
  (instance arm arm (at) (turn z angle) (speed-limits -1 1)))

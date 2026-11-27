; A drive that jumps 30 degrees at t = p, p toleranced: unknown at the jump
; (unbounded speed). The same kink written with min and max holds for
; every p.
(import "std")
(t-span 0 1s)
(param p 0.5 (tol -0.05 0.05))
(body arm (shape (box 10 2 2)))
(assembly root
  (instance j arm (at (move 0 0 0)) (turn z (select (- (/ t 1s) p) (* 10deg/s t) (+ 30 (* 10deg/s t)))) (speed-limits -300deg/s 300deg/s (accuracy 0.01deg/s)))
  (instance m arm (at (move 0 40 0)) (turn z (+ (* 100 (min (/ t 1s) p)) (* 50 (max (- (/ t 1s) p) 0)))) (speed-limits -150deg/s 150deg/s (accuracy 0.01deg/s))))
; expect-status: speed-limits:root/j unknown

; A drive that jumps 30 degrees at t = p, p toleranced: unknown at the jump
; (unbounded speed). The same kink written with min and max holds for
; every p.
(import "std")
(t-span 0 1)
(param p 0.5 (tol -0.05 0.05))
(body arm (shape (box 10 2 2)))
(assembly root
  (instance j arm (at (move 0 0 0)) (turn z (select (- t p) (* 10 t) (+ 30 (* 10 t)))) (speed-limits -300 300 (accuracy 0.01)))
  (instance m arm (at (move 0 40 0)) (turn z (+ (* 100 (min t p)) (* 50 (max (- t p) 0)))) (speed-limits -150 150 (accuracy 0.01))))
; expect-status: speed-limits:root/j unknown

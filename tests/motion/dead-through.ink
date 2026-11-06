; A crank driven from its slider into dead centre and back out (50 cos th
; = 50 + s, s = -20 (1 - t)^2, at dead centre at t = 1): which way it
; turns past there (on through, or back) is not known, so the moments from
; there on are refused and the claim is unknown. Solved afresh from the
; load-time root at each step, it turned back, and the claim held.
(import "std")
(t-span 0 1.5)
(param th 0.3)
(body block (shape (box 10 10 10)))
(body bar (shape (move (box 50 4 4) 25 0 0)))
(assembly (rig (s 0))
  (instance slider block (at (move 50 0 0)) (slide x s))
  (instance crank bar (at) (turn z (rad->deg th))))
(assembly root (instance r (rig (* -20 (square (- 1 t))))))
(solve loop (unknowns th) (== (* 50 (cos th)) (nth (get (frame-of root/r/slider) "origin") 0)))
(assert above (value th > -0.01))
; expect-status: above unknown

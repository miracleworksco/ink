; A solve over a ranged param (q = 100 t + p, p in +-0.01) is certified
; at one moment, not over the motion: the block moving through the wall
; is not proven clear (at each step it is not solved again over p, so
; nothing is decided: unknown).
(import "std")
(t-span 0 1)
(param p 0 (tol -0.01 0.01))
(solve position (unknowns (q 0)) (== q (+ (* 100 t) p)))
(body block (shape (box 2 2 2)))
(body wall (shape (move (box 2 2 2) 50 0 0)))
(assembly root (instance block block (at (move q 0 0))) (instance wall wall))
(assert through-wall (disjoint (shape-of root/block) (shape-of root/wall)))
; expect-status: through-wall unknown

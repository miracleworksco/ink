; A slide whose rate is a ranged param: x = 50 s t, s in 1..2, so the
; block ends anywhere in 50..100. Its margin must use the fastest s: a
; wall at 90 is reached for some s; one at 110 never is.
(import "std")
(t-span 0 1)
(param s 1 (tol 0 1))
(body block (shape (box 2 2 2)))
(body wall (shape (move (box 2 2 2) 90 0 0)))
(body far (shape (move (box 2 2 2) 110 0 0)))
(assembly root (instance block block (at (move (* 50 s t) 0 0))) (instance wall wall) (instance far far))
(assert reaches-wall (disjoint (shape-of root/block) (shape-of root/wall)))
(assert clear-of-far (disjoint (shape-of root/block) (shape-of root/far)))
; expect-status: reaches-wall fails

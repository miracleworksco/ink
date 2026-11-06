; A nonlinear drive, x = 100 t^3: slow at first, 300 mm per unit t at the
; end, where it reaches a block at x = 99 (t > 0.99). The margin must
; follow the rate within each step; one taken from the start would miss it.
(import "std")
(t-span 0 1)
(body a (shape (box 2 2 2)))
(body near (shape (move (box 2 2 2) 99 0 0)))
(body far (shape (move (box 2 2 2) 103 0 0)))
(assembly root (instance a a (at (move (* 100 t t t) 0 0))) (instance near near) (instance far far))
(assert reaches-near (disjoint (shape-of root/a) (shape-of root/near)))
(assert clear-of-far (disjoint (shape-of root/a) (shape-of root/far)))
; expect-status: reaches-near fails

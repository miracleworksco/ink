; A cube slid by p in 0..20 (its tolerance band), a wall moving toward it
; at 20 mm/s: at t = 1, p = 20 the gap is 21.375 - 20.5 = 0.875, so the
; clearance of 1 fails there. Checked as boxes of t and p with one margin
; for both; a margin over p alone, then t's added, let the far points it
; left out at 1 + D come closer by t's (it held).
(import "std")
(t-span 0 1)
(param p 10 (tol -10 10))
(body cube (shape (box 1 1 1)))
(body wall (shape (move (box 8 8 8) 4 0 0)))
(assembly root
  (instance c cube (at (move 0 0 0)) (slide x p))
  (instance w wall (at (move 41.375 0 0)) (slide x (* -20 t))))
; expect-status: clear fails
(assert clear (clearance (shape-of root/c) (shape-of root/w) 1))

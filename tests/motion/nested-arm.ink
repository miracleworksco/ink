; The spinning arm of spin-arm.ink, carried by a parent that turns as much
; again: the arm turns once relative to the carrier, which turns once, so
; the arm's points move by both joints' rates. Its own joint alone would
; make each step's margin half what it is.
(import "std")
(t-span 0 1)
(body arm (shape (move (box 50 2 2) 25 0 0)))
(body post (shape (move (cylinder-z 0.5 4) 40 0 0)))
(body high (shape (move (cylinder-z 0.5 4) 40 0 12)))
(assembly carrier (instance arm arm (at (move 0 0 0)) (turn z (* 180 t))))
(assembly root
  (instance c carrier (at (move 0 0 0)) (turn z (* 180 t)))
  (instance post post)
  (instance high high))
(assert sweeps-post (disjoint (shape-of root/c/arm) (shape-of root/post)))
(assert clears-high (disjoint (shape-of root/c/arm) (shape-of root/high)))
; expect-status: sweeps-post fails

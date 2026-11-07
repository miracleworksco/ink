; An arm spinning up from rest as t^2, through ten turns: it does not
; repeat in t, but the pair of it and a post reads t only through its
; angle, which repeats every turn. A post out of reach is proven clear
; over one turn of the angle, and so over the ten; over t the search
; halves the span only to 1/256 of it, 28 degrees of the turn at its end.
; A post in its reach is hit at a moment it gets there.
(import "std")
(t-span 0 1)
(body arm (shape (box 30 2 2)))
(body post (shape (cylinder-z 0.5 4)))
(assembly line (instance arm arm (at (move 25 0 0))))
(assembly root
  (instance a line (at) (turn z (* 3600 (square t))))
  (instance p post (at (move 41 0 0)))
  (instance q post (at (move 0 35 0))))
(assert clear (disjoint (shape-of root/a/arm) (shape-of root/p)))
(assert hits (disjoint (shape-of root/a/arm) (shape-of root/q)))
; expect-status: hits fails
; expect-status: clear holds

; An arm 50 long and 2 wide (x 0..50) turns once about z over t 0..1; a post
; of radius 0.5 at x = 40 sits on the arm at t = 0 and 1. Sixteen uniform
; pieces of its turn, each proven over all its values, find it there, as do
; the default pieces (motion's tests).
; A post above the arm's plane is clear at every moment: proven.
(import "std")
(t-span 0 1)
(body arm (shape (move (box 50 2 2) 25 0 0)))
(body post (shape (move (cylinder-z 0.5 4) 40 0 0)))
(body high (shape (move (cylinder-z 0.5 4) 40 0 12)))
(assembly root
  (instance arm arm (at (move 0 0 0)) (turn z (* 360 t)))
  (instance post post)
  (instance high high))
(assert sweeps-post (disjoint (shape-of root/arm) (shape-of root/post)))
(assert clears-high (clearance (shape-of root/arm) (shape-of root/high) 1))
; expect-status: sweeps-post fails

; A block jumps 400 mm at each quarter of t (x = 400 floor(4 t)) and
; lands on a post at x = 1600 at t = 1 only. Never proven (the block
; travels 800 times its size, so t ranged is not tried), and no step's
; middle is the moment it lands: the span's end is measured too, and
; fails there, a witness.
(import "std")
(t-span 0 1)
(body block (shape (box 2 2 2)))
(body post (shape (move (box 2 2 2) 1600 0 0)))
(assembly root
  (instance block block (at (move 0 0 0)) (slide x (* 400 (floor (* 4 t)))))
  (instance post post))
(assert lands-on-post (disjoint (shape-of root/block) (shape-of root/post)))
; expect-status: lands-on-post fails

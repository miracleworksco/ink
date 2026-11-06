; A block 2 wide slides 100 mm along x over t 0..1, through a wall 0.2
; thick at x = 50. No step's middle puts the block in the wall (16 steps
; land it at 46.875 and 53.125): only the slide's rate in each step's
; margin finds it. A wall beside the path is clear throughout.
(import "std")
(t-span 0 1)
(body block (shape (box 2 2 2)))
(body wall (shape (move (box 0.2 20 20) 50 0 0)))
(body side (shape (move (box 0.2 20 20) 50 15 0)))
(assembly root
  (instance block block (at (move 0 0 0)) (slide x (* 100 t)))
  (instance wall wall)
  (instance side side))
(assert through-wall (disjoint (shape-of root/block) (shape-of root/wall)))
(assert beside-wall (disjoint (shape-of root/block) (shape-of root/side)))
; expect-status: through-wall fails

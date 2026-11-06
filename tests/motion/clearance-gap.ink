; A block 2 wide slides 10 mm along a wall 0.08 from its face. A proven
; clearance guarantees its distance whatever the tolerance: 0.09 fails (a
; pair 0.08 apart at a step's middle), 0.05 holds at every moment. To the
; model's 0.05 the tolerance used to come off the clearance, and 0.09 was
; proven (no point within 0.04 at each middle, with the step's margin).
(import "std")
(t-span 0 1)
(body block (shape (box 2 2 2)))
(body wall (shape (move (box 40 2 4) 5 2.08 0)))
(assembly root
  (instance block block (at (move 0 0 0)) (slide x (* 10 t)))
  (instance wall wall))
(assert keeps-0.09 (clearance (shape-of root/block) (shape-of root/wall) 0.09))
(assert keeps-0.05 (clearance (shape-of root/block) (shape-of root/wall) 0.05))
; expect-status: keeps-0.09 fails

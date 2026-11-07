; A block slides 100 mm along x by a solve that reads t (q = 100 t),
; through a wall at x = 50: with t ranged, q is enclosed over each piece
; of t by the steps it is followed by, so the claims are proven, a
; failure at a moment q's enclosure there confirms. The clearance claim
; comes first: bounds cached at one moment must be forgotten at the next
; (they read t through the solve), else the block stays where it was.
(import "std")
(t-span 0 1)
(solve position (unknowns (q 0)) (== q (* 100 t)))
(body block (shape (box 2 2 2)))
(body wall (shape (move (box 2 2 2) 50 0 0)))
(body side (shape (move (box 2 2 2) 50 15 0)))
(assembly root (instance block block (at (move q 0 0))) (instance wall wall) (instance side side))
(assert near-wall (clearance (shape-of root/block) (shape-of root/wall) 1))
(assert through-wall (disjoint (shape-of root/block) (shape-of root/wall)))
(assert beside-wall (disjoint (shape-of root/block) (shape-of root/side)))
; expect-status: near-wall fails
; expect-status: through-wall fails
; expect-status: beside-wall holds

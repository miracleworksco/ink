; A block 1.5 x 4 x 4 about x = 1e16 + 1, which is no float: no float
; point is inside it, and the float middle of its box is a face. It is
; not empty, and how deep it can be is its real half-width: so its
; emptiness is unknown, never shown to 0.1.
; expect-status: not-empty unknown
(def bead
  (max (- (abs (- (- x 1e16) 1)) 0.75)
       (- (abs y) 2)
       (- (abs z) 2)))
(assert not-empty (empty bead) (accuracy 0.1))

; A solid undecidable throughout a 4 mm cube (sin^2 + cos^2 - 1: its
; bounds straddle 0 at any cell size) with a ball 0.15 mm in radius deep
; inside it in one corner. The search's first stuck cell ends its proof,
; not its hunt: it goes on for a point of the solid deeper than the
; accuracy and finds one in the ball (it ended unknown at that cell).
; expect-status: e fails
(import "std")
(def zero (- (+ (square (sin x)) (square (cos x))) 1))
(body fog (shape (min (max (box 4 4 4) zero) (move (sphere 0.15) 1.8 1.8 1.8))))
(assert e (empty (shape-of fog)))

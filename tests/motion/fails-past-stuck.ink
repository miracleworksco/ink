; A value undecidable at its accuracy at every moment (sin^2 + cos^2 - 1,
; less the accuracy: its bounds straddle the bound however fine the
; pieces) and failing only within 0.0002 of t = 0.731234. The search's
; first stuck piece ends its proof, not its hunt: it goes on for a
; violation and finds this one (it ended unknown at that piece, and the
; moments measured then missed the violation).
; expect-status: motion fails
(import "std")
(t-span 0 1)
(def zero (- (+ (square (sin t)) (square (cos t))) 1))
(def hit (- (* 1000 (abs (- t 0.731234))) 0.2))
(assert motion (value (min (- zero 0.05) hit) >= 0) (accuracy 0.05))

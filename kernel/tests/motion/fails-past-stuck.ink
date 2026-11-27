; A value undecidable at its accuracy at every moment (sin^2 + cos^2 - 1,
; less the accuracy: its bounds straddle the bound however fine the
; pieces) and failing only within 0.0002 of t = 0.731234. The search's
; first stuck piece ends its proof, not its hunt: it goes on for a
; violation and finds this one (it ended unknown at that piece, and the
; moments measured then missed the violation).
; expect-status: motion fails
(import "std")
(t-span 0 1s)
(def zero (- (+ (square (sin (* 1rad/s t))) (square (cos (* 1rad/s t)))) 1))
(def hit (- (* 1000 (abs (- (/ t 1s) 0.731234))) 0.2))
(require motion (value (min (- zero 0.05) hit) >= 0) (accuracy 0.05))

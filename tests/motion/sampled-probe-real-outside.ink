; A probe truly outside the moving ball at the first sampled moment: a
; guaranteed failure.
; expect-status: inside fails
(import "std")
(t-span 0.01 0.02)
(def ball (move (sphere 1) (+ t 0.999999) 0 0))
(assert inside (subset (probe 0 0 0) ball))

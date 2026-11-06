; q^2 = 0.1 - t has no solution once t > 0.1: those steps are undecided
; (q keeps a value that is no answer there), not held.
(import "std")
(t-span 0 1)
(solve impossible (unknowns (q 0)) (== (square q) (- 0.1 t)))
(assert positive (value (+ q 1) > 0))
; expect-status: positive unknown

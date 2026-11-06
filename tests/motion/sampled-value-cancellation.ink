; A value whose float at a sampled moment rounds away its real one:
; (1e16 + t) - 1e16 is t in real numbers, above 0.001 over the whole span,
; and 0 in floats. Read at a moment by its bounds rounded outward, never a
; guaranteed failure (it was one, value 0).
; expect-status: positive unknown
(t-span 0.01 0.02)
(assert positive (value (- (+ 1e16 t) 1e16) > 0.001) (accuracy 0.000001))

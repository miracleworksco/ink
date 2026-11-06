; value-cancellation's large intermediate named by a def: the same.
; expect-status: positive unknown
(t-span 0.01 0.02)
(def derived (+ 1e16 t))
(assert positive (value (- derived 1e16) > 0.001) (accuracy 0.000001))

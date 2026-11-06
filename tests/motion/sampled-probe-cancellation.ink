; The probe moved by (1e16 + t) - 1e16 alone (t in real numbers, as the
; ball): the cause of far-probe-near-tangent without the far placement.
; expect-status: inside unknown
(import "std")
(t-span 0.01 0.02)
(def shift (- (+ 1e16 t) 1e16))
(def ball (move (sphere 1) (+ t 0.999999) 0 0))
(def p (place (probe 0 0 0) (at (move shift 0 0))))
(assert inside (subset p ball))

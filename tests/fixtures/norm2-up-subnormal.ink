; norm2-up is a length rounded up at every scale: of subnormals too, where
; widening by a factor does not move a float (it gave (norm2-up eta eta) =
; eta, not sqrt(2) eta, and a cone of that size slope sqrt 2 under a promise
; of 1). Subnormals are outside ink's range, so these
; claims are unknown; norm2_up_rounds_up_at_every_scale checks the values.
; expect-status: smallest unknown
; expect-status: subnormal unknown
; expect-status: near-normal unknown
; expect-status: ordinary unknown
; expect-status: ordinary-close unknown
; expect-status: infinite unknown
; expect-status: overflowed unknown
; expect-status: past-largest unknown
; expect-status: zero unknown
(import "std")
(def eta 5e-324)
(assert smallest (value (/ (norm2-up eta eta) eta) >= 1.4142))
(assert subnormal (value (/ (norm2-up 1e-320 1e-320) 1e-320) >= 1.4142))
(assert near-normal (value (/ (norm2-up 1e-310 1e-310) 1e-310) >= 1.4142))
(assert ordinary (value (norm2-up 3 4) >= 5))
(assert ordinary-close (value (norm2-up 3 4) < 5.0000000001))
; And past the largest float: an infinite or overflowing length is inf
; (inf / inf once made it the smallest float).
(assert infinite (value (norm2-up inf 1) >= 1e308))
(assert overflowed (value (norm2-up (- 1e308 -1e308) 1) >= 1e308))
(assert past-largest (value (norm2-up 1.5e308 1.5e308) >= 1e308))
(assert zero (value (norm2-up 0 0) > 0))

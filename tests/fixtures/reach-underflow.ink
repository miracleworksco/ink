; A cube of half-side 1e-200: its reach from the z axis is sqrt(2) 1e-200,
; whose squares underflow in floats. The reach bounds it from above. 1e-200
; is outside ink's range: unknown.
; expect-status: covers unknown
(import "std")
(def tiny
  (max (- (abs x) 1e-200)
       (- (abs y) 1e-200)
       (- (abs z) 1e-200)))
(def reach (profile-reach tiny))
(assert covers (value reach >= 1.4142135623730951e-200))

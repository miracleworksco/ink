; far-probe-direct with the probe's offset written (1e16 + t) - 1e16, equal
; to the ball's t in real numbers: inside at every moment, frozen 0.01 out
; in floats.
; expect-status: inside unknown
(import "std")
(t-span 0.01 0.02)
(def offset (+ 1e12 t))
(def rounded-offset (+ 1e12 (- (+ 1e16 t) 1e16)))
(def ball (move (sphere 1) offset 0 0))
(def p (place (probe -0.999999 0 0) (at (move rounded-offset 0 0))))
(assert inside (subset p ball))

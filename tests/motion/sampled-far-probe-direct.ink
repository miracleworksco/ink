; A probe 1e-6 inside a ball, both carried 1e12 + t out: inside at every
; moment. A probe moved by t is a moving point: its coordinates composed
; into the field, bounded outward, never the field at their floats (frozen
; where they round, it failed, guaranteed).
; expect-status: inside unknown
(import "std")
(t-span 0.01 0.02)
(def offset (+ 1e12 t))
(def ball (move (sphere 1) offset 0 0))
(def p (place (probe -0.999999 0 0) (at (move offset 0 0))))
(assert inside (subset p ball))

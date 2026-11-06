; A load-time number computed through an overflow is still the number it
; is: (min (exp 1000) 1) is 1, as a param's default and a res.
; exp 1000 overflows, a number floats lost: outside ink's range, unknown.
; expect-status: finite unknown
; expect-status: volume unknown
(import "std")
(param p (min (exp 1000) 1))
(assert finite (value p > .5))
(assert volume (volume (sphere 1) > 1) (res (min (exp 1000) .1)))

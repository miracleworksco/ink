; A real collision found at a sampled moment (t = 0.015625, two spheres
; near tangent 1e8 out): the same claim asked of its own nodes at that
; moment fails too, so it stays a guaranteed failure with its witness in
; the world. Asked to 1e-6 mm 1e8 mm out, it is outside ink's range: unknown.
; expect-status: apart unknown
(import "std")
(t-span 0 1)
(def shift (+ 1e8 (* 0.001 t) 0.0001))
(def a (move (sphere 1) shift 0 0))
(def b (move (sphere 1) (+ shift 1.99999) (* 100 (- t 0.015625)) 0))
(assert apart (disjoint a b) (accuracy 0.000001))

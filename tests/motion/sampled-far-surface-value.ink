; A point a micrometre inside a ball carried 1e12 + t out: its field is
; -0.0001 t in real numbers, below -5e-7 at every moment; the far
; translation's floats lose the offset. Not a guaranteed failure (it was).
; expect-status: material unknown
(import "std")
(t-span 0.01 0.02)
(def offset (+ 1e12 t))
(def ball (move (sphere 1) offset 0 0))
(def near-surface (warp ball (- offset (- 1 (* 0.0001 t))) 0 0))
(assert material (value near-surface < -0.0000005) (accuracy 0.00000001))

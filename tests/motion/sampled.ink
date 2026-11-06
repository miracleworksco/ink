; A baked (sampled) sphere moving along x, far from another: its distance
; is the samples', so the claim is measured at each step, not proven.
(import "std")
(t-span 0 1)
(def f (grid-bake (sphere 1) (res 0.5)))
(body a (shape f))
(body b (shape (move (sphere 1) 0 0 10)))
(assembly root (instance a a (at (move t 0 0))) (instance b b))
(assert separated (disjoint (shape-of root/a) (shape-of root/b)))
; expect-status: separated holds estimate

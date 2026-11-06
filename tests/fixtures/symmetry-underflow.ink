; A strip whose two forms have the same factors, (p a) b and p b with a,
; b and c applied in another order: a product of a param and constants is
; its float, and at p = 1e-200 (p a) b underflows to 0 while p b does not.
; The strip is then |v| < 0.01, not round, and turning it about z sweeps
; it through the post: no symmetry may excuse the turn. 1e-200 is outside
; ink's range: unknown.
; expect-status: no-collision unknown
(import "std")
(param p 1e-200)
(def a 1e-200)
(def b 1e300)
(def c 1e100)
(def u (* (* x (* (* p a) b)) c))
(def v (* (* (* y (* p b)) a) c))
(def strip (- (sqrt (+ (square u) (square v))) 0.01))
(body rotor (shape (max strip (cylinder-z 2 1))))
(body post (shape (sphere 0.006)))
(t-span 0 1)
(assembly root
  (instance r rotor (at) (turn z (* 360 t)))
  (instance p post (at (move 0 1.5 0))))
(assert no-collision
  (disjoint (shape-of root/r) (shape-of root/p))
  (accuracy 0.001))

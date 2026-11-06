; promise-param-overflow.ink with t warped to 1e200: the rewritten k is
; built of literals, 0.001 + inf / inf, and gap held, guaranteed.
; 1e200 is outside ink's range: unknown, never the false holds.
; expect-status: gap unknown
(import "std")
(def q (+ 1 (square t)))
(def a (warp (lipschitz (+ .001 (/ q q)) (sphere 1)) x y z 1e200))
(def b (move (sphere .01) 1.05 0 0))
(assert gap (clearance a b .1))

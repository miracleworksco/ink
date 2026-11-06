; A promise true for every real p: k is 0.001 + q / q, 1.001. At
; --param p=1e200 (or the range 1e200:1e200) q overflows, q / q is inf /
; inf, 0, and k was 0.001: gap held, guaranteed, between points 0.042
; apart. A proof takes no promise whose
; constants are not real at the query's params (commands::claims).
; expect-status: gap fails
(import "std")
(param p 1)
(def q (+ 1 (square p)))
(def a (lipschitz (+ .001 (/ q q)) (sphere 1)))
(def b (move (sphere .01) 1.05 0 0))
(assert gap (clearance a b .1))

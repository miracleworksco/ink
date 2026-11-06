; The bound q = 10 - 100t is solved at each step, not frozen at t = 0:
; t < q fails once t passes about 0.099.
(import "std")
(t-span 0 1)
(solve position (unknowns (q 10)) (== q (- 10 (* 100 t))))
(assert below (value t < q))
; expect-status: below fails estimate

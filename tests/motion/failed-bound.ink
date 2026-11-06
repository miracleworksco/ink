; A bound whose solve fails once t > 0.1 (q^2 = 0.1 - t): those steps are
; undecided, not held.
(t-span 0 1)
(solve impossible (unknowns (q 1)) (== (square q) (- 0.1 t)))
(assert below (value (* 0.1 t) < q))
; expect-status: below unknown

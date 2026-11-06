; overflow-promise.ink with every number written well inside 1e+-100: the
; ends are computed (1e308), and their difference still overflows.
; expect-warning: (value (- z1 z0) < inf)
; expect-status: gap unknown
(import "std")
(def e (* 1e100 1e100 1e100 1e8))
(def a (max (- (cone-z 1 1 (neg e) e) 0.1) (box 4 4 4)))
(def b (move (sphere 0.01) 1.05 0 0))
(assert gap (clearance a b 0.1))

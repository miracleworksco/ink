; An L plate whose re-entrant edge lies on the cells' planes (x = y = 0),
; bent by a load on its arm away from the clamp. The edge is sharp: its
; elastic stress is unbounded, so `finite-peak` and `nonzero-safety` are
; false, however low a grid reads them.
(import "std")
(def raw (subtract (box 20 20 4) (move (box 20 20 8) 10 10 0)))
(def part (rotate-z raw 0))
(fn reg (f) (rotate-z f 0))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa))
 (fix wall (reg (+ x 9)))
 (load end (force 0 1N 0) (region (reg (- 9 x)))) (res 1))
(require finite-peak (max-stress s < 1MPa))
(require nonzero-safety (min-safety-factor s > 1))

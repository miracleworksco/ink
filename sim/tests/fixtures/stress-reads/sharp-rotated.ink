; The corner of sharp-L turned 17 degrees, off the cells' planes. Sharp, so
; its stress is unbounded: `finite-peak` and `elastic-safe` are false,
; `corner-on-own-grid` true, however a grid reads them.
(import "std")
(def raw (subtract (box 10 10 4) (move (box 10 10 8) 5 5 0)))
(def part (rotate-z raw 17))
(fn reg (f) (rotate-z f 17))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa))
 (fix wall (reg (+ x 4.5)))
 (load end (force 0 0.01N 0) (region (reg (- 4.5 x)))) (res 1))
(require finite-peak (max-stress s < 1MPa))
(require elastic-safe (min-safety-factor s > 1))
(require corner-on-own-grid (max-stress s > 1MPa))

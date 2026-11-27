(import "std")
; A plate with a fin thinner than its cells: split cells, and a need that is not monotone in the cell size.
(def part (min (box 60 60 4) (move (box 40 0.3 10) 0 0.07 6.9)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3)) (fix root (- (+ x 30) 1)) (load tip (force 0 0 -10N) (region (- 29 x))))

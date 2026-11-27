(import "std")
(def part (lattice-fill (box 60 60 60) (gyroid 12) 1.8))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa)) (fix root (- (+ x 30) 1)) (load tip (force 0 0 -100N) (region (- 29 x))))

(import "std")
(def part (box 100 100 3))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa)) (fix root (- (+ x 50) 1)) (load tip (force 0 0 -100N) (region (- 49 x))))

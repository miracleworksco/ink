(import "std")
(def part (subtract (cylinder-z 40 30) (cylinder-z 25 40)))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa)) (fix root (- x -38)) (load tip (force 0 0 -100N) (region (- 38 x))))

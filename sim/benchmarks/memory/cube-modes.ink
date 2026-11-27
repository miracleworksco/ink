(import "std")
(def part (box 100 100 100))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (density 7.8g/cm3)) (fix root (- (+ x 50) 1)) (modes 4))

(import "std")
(def part (box 100 100 100))
(study s (shape part) (heat (k 50W/m/K)) (temperature root (- (+ x 50) 1) 20degC) (source h (power 10W) (region (- 49 x))))

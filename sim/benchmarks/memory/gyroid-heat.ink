(import "std")
(def part (lattice-fill (box 60 60 60) (gyroid 12) 1.8))
(study s (shape part) (heat (k 50W/m/K)) (temperature root (- (+ x 30) 1) 20degC) (source h (power 10W) (region (- 29 x))))

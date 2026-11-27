(import "std")
(def part (box 100 100 100))
(study s (shape part) (heat (k 50W/m/K) (density 7800kg/m3) (capacity 500J/kg/K)) (temperature root (- (+ x 50) 1) 20degC) (source h (power 10W) (region (- 49 x))) (transient (duration 60s) (steps 10) (initial 20degC)))

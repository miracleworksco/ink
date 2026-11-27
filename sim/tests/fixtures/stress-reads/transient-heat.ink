; An insulated rod held at 20 C at one end, heated 1 W at the other: at
; 1000 s its end is 44.19 C (one-dimensional; steady well before 1000 s).
; `true-lower` (> 44) is true.
(import "std")
(def part (move (box 20 4 4) 10 0 0))
(study s (shape part) (heat (k 50W/m/K) (density 7800kg/m3) (capacity 500J/kg/K))
 (temperature wall (- x 0.3) 20degC)
 (source tip (power 1W) (region (- 19.3 x)))
 (transient (duration 1000s) (steps 100) (initial 20degC)) (res 1))
(require true-lower (max-temperature s > 44degC))

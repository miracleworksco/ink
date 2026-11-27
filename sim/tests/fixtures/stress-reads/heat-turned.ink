; The heat bar of heat-interface-under-a-cell turned 2 degrees about z: a
; 20 x 10 x 2 bar along u, k 0.2 to u = 10.08 and 200 beyond, held at 20 C
; over its first 0.01 mm, 1 W in its last mm: exactly 1-D along u, 2539.855
; C at its end. Cells of 0.1, 0.2 and 0.4 mm read 2530.91, 2525 and 2510
; (grid error 6.55): 9 K low past it, the face's stair of cells at 2
; degrees no plane all three grids share. True: turned-over-2538-6.
; False: turned-under-2538-6.
(import "std")
(def u (+ 10 (* 0.999390827019 (- x 10)) (* 0.034899496703 y)))
(def v (- (* 0.999390827019 y) (* 0.034899496703 (- x 10))))
(def part (max (neg u) (- u 20) (- (abs v) 5) (- (abs z) 1)))
(study hot (shape part) (heat (k (* 1W/m/K (select (- u 10.08) 0.2 200))))
 (temperature cold (- u 0.01) 20degC)
 (source chip (power 1W) (region (- 19 u))) (res 0.1))
(require turned-under-2538-6 (max-temperature hot < 2538.585073degC))
(require turned-over-2538-6 (max-temperature hot > 2538.585073degC))

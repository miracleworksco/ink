; band.ink's bar turned 1 degree about z, its soft band (E 2000 in 200000)
; from u 8.12 to 16.12, pulled 0.01 mm along u (nu 0): uniform 2.451 MPa
; along u, a safety factor of 102, exactly. Cells split to 0.125 mm along
; the face's stair read 88.3 MPa there (grid error 11.4). True:
; turned-sf-over-95-9. False: turned-over-2-6.
(import "std")
(def u (+ 12 (* 0.999847695156 (- x 12)) (* 0.017452406437 (- y 5)) (* 0.000000000000 (- z 5))))
(def v (+ 0 (* -0.017452406437 (- x 12)) (* 0.999847695156 (- y 5)) (* 0.000000000000 (- z 5))))
(def w (+ 0 (* -0.000000000000 (- x 12)) (* -0.000000000000 (- y 5)) (* 1.000000000000 (- z 5))))
(def part (max (neg u) (- u 24) (- (abs v) 5) (- (abs w) 5)))
(study s (shape part) (elastic (E (* 1MPa (select (max (- 8.119408967 u) (- u 16.11940897)) 2000 200000))) (nu 0) (yield 250MPa))
 (fix a (- u 0.001) (move 0 0 0))
 (fix b (- 23.999 u) (move 0.009998476952 0.0001745240644 0)) (res 0.5))
(require turned-over-2-6 (max-stress s > 2.598039216MPa (excluding support-margin)))
(require turned-sf-over-95-9 (min-safety-factor s > 95.88 (excluding support-margin)))

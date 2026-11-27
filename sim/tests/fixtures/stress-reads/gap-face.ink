; Two plates 14 x 3 x 0.8 (nu 0) across an air gap of 0.125 mm, z 0..0.8 and
; 0.925..1.725, each bent by its ends' moves about its own mid-plane: sxx =
; E k (z - zc) in each, exactly (the gap's faces carry nothing), 8 MPa at
; every face, a safety factor of 31.25. Not glued: on 0.5 mm cells split to
; 0.0625 the gap's sides move apart as the closed form has them. The cells
; beside the gap read 11.29 MPa where the upper plate's lower face meets its
; side (sxx -10.3, szz 2.0 where the truth is -8 and 0), and the three grids,
; split to 0.0625, 0.125 and 0.25 mm, converge away from the truth (safety
; factor 22.15, 24.72, 28.71). True: gap-sf-over-30-6. False:
; gap-sf-under-30-6.
(import "std")
(def slab (max (neg x) (- x 14) (neg y) (- y 3)))
(def part (max slab (min (max (- z 0.8) (- 0 z)) (max (- z 1.725) (- 0.925 z)))))
(study s (shape part) (elastic (E 200000MPa) (nu 0) (yield 250MPa))
 (fix a1 (max (- x 0.001) (- z 0.8)) (move (* 0.0001 x (- z 0.4)) 0 (* -5e-05 x x)))
 (fix a2 (max (- x 0.001) (- 0.925 z)) (move (* 0.0001 x (- z 1.325)) 0 (* -5e-05 x x)))
 (fix b1 (max (- 13.999 x) (- z 0.8)) (move (* 0.0001 x (- z 0.4)) 0 (* -5e-05 x x)))
 (fix b2 (max (- 13.999 x) (- 0.925 z)) (move (* 0.0001 x (- z 1.325)) 0 (* -5e-05 x x))) (res 0.5))
(require gap-sf-under-30-6 (min-safety-factor s < 30.625 (excluding support-margin)))
(require gap-sf-over-30-6 (min-safety-factor s > 30.625 (excluding support-margin)))

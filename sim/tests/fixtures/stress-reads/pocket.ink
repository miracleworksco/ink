; A 60 x 40 x 6 plate with a 20 x 20 pocket 3 deep, every pocket edge r1.5
; (floor and vertical corners), its rim sharp (convex). Cantilevered, bent in
; z, on cells of 1 mm. CalculiX 23.81 -> 24.28 MPa, rising to about 24.7
; (yield 250: a safety factor of 10.1 to 10.5). True: lt-28, lt-33, gt-22,
; sf-gt-9, some, nonzero-safety. False: lt-22, sf-gt-11.
(import "std")
(def cavity (move (rounded-box 20 20 6 1.5) 0 0 3))
(def part (subtract (box 60 40 6) cavity))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa))
 (fix wall (+ x 29))
 (load end (force 0 0 50N) (region (- 29 x))) (res 1))
(require lt-22 (max-stress s < 22MPa))
(require lt-28 (max-stress s < 28MPa))
(require lt-33 (max-stress s < 33MPa))
(require gt-22 (max-stress s > 22MPa))
(require sf-gt-9 (min-safety-factor s > 9))
(require sf-gt-11 (min-safety-factor s > 11))
(require some (max-stress s > 0.01MPa))
(require nonzero-safety (min-safety-factor s > 1))

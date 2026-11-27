; A 200 x 40 x 4 plate with a d10 hole, 100 MPa gross, modelled as an eighth
; on symmetry planes: CalculiX 327.12 MPa at the hole. True: lt-335, lt-345,
; lt-360, gt-310. False: lt-310, lt-320, gt-335. On cells of 1 mm.
(import "std")
(def full (subtract (box 200 40 4) (cylinder-z 5 10)))
(body plate (shape (max full (neg x) (neg y) (neg z))) (attr material "steel"))
(study pull (shape (shape-of plate)) (elastic (material steel))
  (fix sym-x (- x 0.01) (dofs x))
  (fix sym-y (- y 0.01) (dofs y))
  (fix sym-z (- z 0.01) (dofs z))
  (load pull (force 4000N 0 0) (region (- 99 x)))
  (res 1))
(require lt-310 (max-stress pull < 310MPa))
(require lt-320 (max-stress pull < 320MPa))
(require lt-335 (max-stress pull < 335MPa))
(require lt-345 (max-stress pull < 345MPa))
(require lt-360 (max-stress pull < 360MPa))
(require gt-310 (max-stress pull > 310MPa))
(require gt-335 (max-stress pull > 335MPa))

; A pressurised quarter ring on symmetry planes (benchmark 13): r 49..51, 40
; long, aluminium, 2 MPa on the bore alone. Lame (plane stress) at the bore:
; von Mises 51.05 MPa; CalculiX (C3D10, 3 meshes) 51.06. Its cut cells' own
; reads scatter 40..54 MPa on the bore at every cell size, so a peak read from
; them sits near the top of that scatter on every grid, and the three grids
; agree. True: vessel-under-51-1. False: vessel-peak-over-52, vessel-over-51-1.
(import "std")
(def ri 49) (def ro 51) (def L 40)
(def ring (max (extrude-z (ring-2d ri ro) 0 L) (neg x) (neg y)))
(body shell (shape ring) (attr material "aluminum-6061"))
(study pressure-true (shape (shape-of shell))
  (elastic (material aluminum-6061))
  (fix sym-x (- x 0.01) (dofs x))
  (fix sym-y (- y 0.01) (dofs y))
  (fix sym-z (- z 0.01) (dofs z))
  (load inside (pressure 2MPa) (region (max (- (len2 x y) (+ ri 0.2)) (- 0.01 x) (- 0.01 y) (- 0.01 z) (- z 39.99))))
  (res 0.4))
(require vessel-peak-over-52 (max-stress pressure-true > 52MPa))

(require vessel-over-51-1 (max-stress pressure-true > 51.1MPa))
(require vessel-under-51-1 (max-stress pressure-true < 51.1MPa))

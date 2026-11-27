; Benchmark 13: thin-walled pressure vessel shell vs p r / t.
;
; Task: the cylindrical shell of an aluminium vessel, r_i = 49, t = 2 (r/t =
; 25), internal pressure 2 MPa. Modelled as a 40 mm ring slice, quarter
; symmetry (x >= 0, y >= 0) with symmetry supports (x held on the x = 0
; face, y on the y = 0 face, z on the bottom face; top free: plane stress,
; an open-ended cylinder).
;
; Study `pressure-true` loads the bore with (pressure p). Study `pressure`
; approximates it by 6 sector loads of 15 deg each on a 0.5 mm band at the
; bore, each along its sector's bisector with magnitude p r_i (pi/12) L =
; 1026.25 N.
;
; Hand calc (Lame, plane stress, E = 68.9 GPa, nu = 0.33):
;   hoop stress p r_i / t = 49.0 MPa (Lame: 50.0 at the bore, 48.0 outside)
;   radial displacement at r = 50:
;     p r_i^2 / (E (r_o^2 - r_i^2)) ((1 - nu) r + (1 + nu) r_o^2 / r) = 0.03579 mm
;
; Pass criteria (ink solve --points at mid-wall, 45 deg, mid-height):
;   - radial displacement within 5% of 0.03579 mm;
;   - von Mises at mid-wall within 10% of 49.5 MPa (hoop 49, radial -1).
;
; Both studies must pass; `pressure-true` also at 0 deg (on the x axis).

(import "std")

(param p 2)
(def ri 49) (def ro 51) (def L 40)
(def ring (max (extrude-z (ring-2d ri ro) 0 L) (neg x) (neg y)))
(body shell (shape ring) (attr material "aluminum-6061") (attr color "#A8B8C8"))

; material within 0.5 mm of the bore, between angles a0 and a1 (deg)
(fn sector-band (a0 a1)
  (max (- (len2 x y) (+ ri 0.5))
       (- (* (sin (deg->rad a0)) x) (* (cos (deg->rad a0)) y))
       (- (* (cos (deg->rad a1)) y) (* (sin (deg->rad a1)) x))))
(def f-sector (* p ri (/ pi 12) L))
(fn fx (a) (* f-sector (cos (deg->rad a))))
(fn fy (a) (* f-sector (sin (deg->rad a))))

(study pressure (shape (shape-of shell))
  (elastic (material aluminum-6061))
  (mirror sym-x x 0)
  (mirror sym-y y 0)
  (fix axial-datum (- z 0.25) (dofs z) (support))
  (load p0 (force (* 1N (fx 7.5)) (* 1N (fy 7.5)) 0) (region (sector-band 0 15)))
  (load p1 (force (* 1N (fx 22.5)) (* 1N (fy 22.5)) 0) (region (sector-band 15 30)))
  (load p2 (force (* 1N (fx 37.5)) (* 1N (fy 37.5)) 0) (region (sector-band 30 45)))
  (load p3 (force (* 1N (fx 52.5)) (* 1N (fy 52.5)) 0) (region (sector-band 45 60)))
  (load p4 (force (* 1N (fx 67.5)) (* 1N (fy 67.5)) 0) (region (sector-band 60 75)))
  (load p5 (force (* 1N (fx 82.5)) (* 1N (fy 82.5)) 0) (region (sector-band 75 90)))
  (res 0.4))

; The same with the pressure load master added during this benchmark: a
; traction -p n on the surface inside the region (bore only, clear of the
; symmetry and end faces).
(study pressure-true (shape (shape-of shell))
  (elastic (material aluminum-6061))
  (mirror sym-x x 0)
  (mirror sym-y y 0)
  (fix axial-datum (- z 0.25) (dofs z) (support))
  ; (the bore alone: a micrometre in from its rims. A nanometre was under
  ; what a surface point's place rounds to, and slivers of the end faces
  ; were pressed: no grid but the claim's own read the load whole.)
  (load inside (pressure (* 1MPa p)) (region (max (- (len2 x y) (+ ri 0.2)) (- 0.001 x) (- 0.001 y) (- 0.001 z) (- z (- L 0.001)))))
  (res 0.4))

(require holds (min-safety-factor pressure > 3))
(require holds-true (min-safety-factor pressure-true > 3))

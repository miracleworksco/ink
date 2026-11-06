; Checks a desk lamp will not tip over: its centre of mass stays over the
; base on a 10 degree slope facing any way, with the arm swung any way
; about the post, and a long-armed copy tips. The claims are integrals:
; estimates, sampled over the swing. The lamp is a steel disc base, an
; aluminium post and arm, and a printed head hanging from the arm's end.
;
; Frame: the base's bottom face is z = 0, its axis +z; the arm reaches
; along +x. Densities in g/mm^3, so masses are in grams.
;
; Hand calculation (V in mm^3, m in g, centres in mm):
;   base  steel disc r 75, 15 thick:  V = pi 75^2 15  = 265 072  m = 2080.8  at (0, 0, 7.5)
;   post  aluminium rod r 8, 300 long: V = pi 8^2 300 = 60 319   m = 162.9   at (0, 0, 165)
;   arm   aluminium bar 300 x 16 x 10: V = 48 000             m = 129.6   at (158, 0, 310)
;   head  PLA cylinder r 40, 60 long:  V = pi 40^2 60 = 301 593  m = 374.0   at (280, 0, 275)
;   total m = 2747.2; centre x = (129.6 158 + 374.0 280) / 2747.2 = 45.57,
;   y = 0, z = (2080.8 7.5 + 162.9 165 + 129.6 310 + 374.0 275) / 2747.2 = 67.52.
; Tipping: the lamp stands while the vertical through its centre of mass
; meets the floor inside the base's footprint, the r 75 disc of its bottom
; face. Tilted by b toward a direction u, that vertical meets the base's
; plane at (x, y) + z tan(b) u (the centre is above the base, z > 0). The
; regular octagon inscribed in the disc (apothem 75 cos 22.5 = 69.29) lies
; inside it. So it is enough that for each of its eight normals n,
; n.(x, y) + z tan(b) < 69.29, whatever u is.
; At b = 10: along +x, 45.57 + 67.52 0.1763 = 57.48 < 69.29; the other
; normals give less. Flat, it tips over at atan((75 - 45.57) / 67.52) =
; 23.6 degrees toward +x.
;
; The control: the arm 550 long (m = 237.6, centre x 283) and the head at
; x = 530: total 2855.2, centre x = (237.6 283 + 374.0 530) / 2855.2 =
; 92.97, beyond the base's edge at 75: it tips over even on a level floor.
;
; Swung: the arm and head turned about the post's axis by any angle s in
; [0, 360] (the param swing). The mass is the same at every s (2747.2);
; the centre of mass circles the post, (45.57 cos s, 45.57 sin s, 67.52),
; so along each of the octagon's normals it is at most 45.57 < 69.29
; level, and 45.57 + 11.90 = 57.48 < 69.29 on a 10 degree slope: it
; stands whichever way the arm points. The long-armed lamp's centre
; circles at 92.97: at s = 180 it lies 92.97 along the normal facing -x.
; Each swung body is integrated in its own frame, where it does not move.
;
; Each claim is a weighted mean, decided without dividing bounds. With W
; the integral of density and A that of density (h - v), mean - v = A / W.
; W > 0, so the claim is A's sign (docs/IR.md, Claims).

(import "std")

(def steel 7.85g/cm3)
(def aluminium 2.7g/cm3)
(def pla 1.24g/cm3)

(body base (shape (cylinder-z 75 15)) (attr density steel) (attr color (vec 0.3 0.3 0.32)))
(body post (shape (cylinder-z 8 300)) (attr density aluminium) (attr color (vec 0.8 0.8 0.82)))
(fn arm-bar (length) (move (box length 16 10) (+ 8 (* 0.5 length)) 0 0))
(body arm (shape (arm-bar 300)) (attr density aluminium) (attr color (vec 0.8 0.8 0.82)))
(body long-arm (shape (arm-bar 550)) (attr density aluminium) (attr color (vec 0.8 0.8 0.82)))
(body head (shape (cylinder-z 40 60)) (attr density pla) (attr color (vec 0.95 0.9 0.7)))

; Each part is modelled centred on its own origin, and placed.
(assembly lamp
  (instance base base (at (move 0 0 7.5)))
  (instance post post (at (move 0 0 165)))
  (instance arm arm (at (move 0 0 310)))
  (instance head head (at (move 280 0 275))))

(assembly long-lamp
  (instance base base (at (move 0 0 7.5)))
  (instance post post (at (move 0 0 165)))
  (instance arm long-arm (at (move 0 0 310)))
  (instance head head (at (move 530 0 275))))

; The horizontal coordinate along the octagon's k-th normal, 45 k degrees
; from +x, and the octagon's apothem.
(fn along (k) (+ (* (cos (deg->rad (* 45 k))) x) (* (sin (deg->rad (* 45 k))) y)))
(def octagon (map along (range 8)))
(def apothem (* 75 (cos (deg->rad 22.5))))
(def tilt (tan (deg->rad 10)))

(assert mass-under-3kg (integral lamp (attr density) < 3000))
(assert centre-above-the-base (mean lamp z (weight (attr density)) > 0))
(assert centre-below-70 (mean lamp z (weight (attr density)) < 70))
; Level, the centre of mass is over the octagon; on a 10 degree slope, the
; vertical through it still meets the base's plane inside it: the mean of
; each of the eight fields (a vec) is below the apothem.
(assert stands-level (mean lamp octagon (weight (attr density)) < apothem))
(assert stands-on-10-degrees (mean lamp (map (lambda (n) (+ n (* tilt z))) octagon) (weight (attr density)) < apothem))

; The control: level, with the long arm, the centre of mass is past the
; base's edge, so the claim that it stands fails.
; expect-status: long-arm-stands fails estimate
(assert long-arm-stands (mean long-lamp octagon (weight (attr density)) < apothem))
(assert long-arm-tips-over (mean long-lamp x (weight (attr density)) > 75))

; Swung: the arm and head turned about the post through any angle.
(def swing (ranged 0 360))
(assembly swung-lamp
  (instance base base (at (move 0 0 7.5)))
  (instance post post (at (move 0 0 165)))
  (instance arm arm (at (move 0 0 310) (rot-z swing)))
  (instance head head (at (move 280 0 275) (rot-z swing))))
(assembly swung-long-lamp
  (instance base base (at (move 0 0 7.5)))
  (instance post post (at (move 0 0 165)))
  (instance arm long-arm (at (move 0 0 310) (rot-z swing)))
  (instance head head (at (move 530 0 275) (rot-z swing))))

(assert swung-mass-under-3kg (integral swung-lamp (attr density) < 3000))
(assert swung-stands-level (mean swung-lamp octagon (weight (attr density)) < apothem))
(assert swung-stands-on-10-degrees (mean swung-lamp (map (lambda (n) (+ n (* tilt z))) octagon) (weight (attr density)) < apothem))
; expect-status: swung-long-arm-stands fails estimate
(assert swung-long-arm-stands (mean swung-long-lamp octagon (weight (attr density)) < apothem))

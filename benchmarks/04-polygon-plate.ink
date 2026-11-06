; Benchmark 04: extrude an arbitrary 12-vertex non-convex polygon into a plate.
;
; Task: an 8 mm plate from a 12-vertex outline with a notch, a slanted edge
; and a V-cut, using the built-in `polygon` field.
;
; Reference: shoelace area 3352.5 mm^2, volume 26 820 mm^3.
;
; Pass criteria (ink assert --res 0.1): the volume estimate within 3.5% of
; 26 820 (plate-volume-lo and plate-volume-hi).

(import "std")

(param thick 8)            ; plate thickness (mm)

; Vertices (CCW): an irregular L-ish outline with a notch, a slanted edge and a
; V-cut.
(def outline
  (vec (vec 0 0) (vec 90 0) (vec 90 20) (vec 70 20) (vec 70 12) (vec 55 12)
       (vec 40 35) (vec 40 70) (vec 22 70) (vec 22 50) (vec 10 60) (vec 0 40)))

(def profile (polygon outline))
(def plate (extrude-z profile 0 thick))

(body plate
  (shape plate)
  (attr color "#4C8BD6")
  (attr density 2.7g/cm3))                     ; aluminium, g/mm^3

; Analytic: shoelace area 3352.5 mm^2, so volume = 3352.5 * 8 = 26820 mm^3.
; `ink mass --res 0.1` estimates 26820.1. The claims allow +-3.5%, measured at
; --res 0.1; a volume is an estimate at any resolution.
(assert plate-volume-lo (volume (shape-of plate) > 25880))
(assert plate-volume-hi (volume (shape-of plate) < 27760))
; The plate's faces lie on its envelope's: no search separates the faces at
; the model's floor within its budget, so it is decided at cells of 0.2 (to
; 0.4 mm).
(assert in-envelope (within (shape-of plate) (0 0 0) (90 70 8)) (res 0.2))
; The notch and the V-cut are empty.
(assert notch-empty (disjoint (shape-of plate) (move (box 13 6 20) 62.5 16.5 4)))
(assert vcut-empty (disjoint (shape-of plate) (move (box 4 4 20) 15 66 4)))

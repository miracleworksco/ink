; Benchmark 05: an L bracket with named-feature intent.
;
; Task: load "the top face"; support "the face where the upright meets the
; wall"; fillet "the edge where the arm meets the upright". The references
; must survive parameter changes. Every reference is derived from the same
; datum planes that build the geometry.
;
; Frames: wall is the plane y = 0 (+y out of the wall), z up, bracket centred
; on x = 0. Upright against the wall; arm along +y at the top.
;
; Pass criteria (ink require):
;   - each face reference selects material: over 100 mm^3 of the 1200 / 1500
;     expected (top-face-selects, wall-face-selects; proven);
;   - the bracket holds under the shelf load: safety factor > 2 (numerically).
;     On cells of 4, 2 and 1 mm it reads 4.472, 4.140 and 4.115, converging
;     faster than the elements can, so its grid error is three times that
;     movement, 1.07, inside its margin of 2.12; on 2, 1 and 0.5 mm, 4.497,
;     4.115 and 3.842, converging, 3.842 +- 0.86. (Not settled while the
;     coarser grids read 5.002 and 4.282, and on the finer three moving
;     apart: cells under 3% full were dropped, and a coarse grid split
;     where material that was not there looked thin.)
; expect-status: holds holds estimate

(import "std")

(param width 30)            ; x
(param thick 5)             ; plate thickness (upright and arm)
(param height 50)           ; upright height (z)
(param reach 40)            ; arm reach from the wall (y)
(param fillet 4)            ; inside-corner fillet radius
(param lip 0)               ; optional lip on the arm's tip (height above the arm); 0 = none

; ---------------------------------------------------------------- datums
; Named planes: the only "feature names" available. Each is a def (a number),
; and both the geometry and every face reference use it.
(def wall-plane 0)                              ; y of the upright's back face
(def upright-front (+ wall-plane thick))        ; y of the upright's front face
(def top-plane height)                          ; z of the arm's top face
(def arm-under (- height thick))                ; z of the arm's underside

; ---------------------------------------------------------------- geometry
(def slab-x (- (abs x) (* 0.5 width)))          ; |x| <= width/2
(def upright (max slab-x (- wall-plane y) (- y upright-front) (neg z) (- z top-plane)))
(def arm (max slab-x (- wall-plane y) (- y reach) (- arm-under z) (- z top-plane)))
(def lip-block (max slab-x (- (- reach 3) y) (- y reach) (- top-plane z) (- z (+ top-plane lip))))

; Fillet "the edge where the arm meets the upright": no edge exists to name.
; std corner-fillet needs the distances from the two walls that form the edge,
; i.e. the author restates both planes (upright-front, arm-under) and clips
; the fillet to the width by hand.
(def corner (max (corner-fillet (- y upright-front) (- arm-under z) fillet) slab-x))

; The optional lip. Params are load-time values now, so `if` builds it only
; when it exists. (A runtime `(select (neg lip) lip-block 1000)` also works,
; but then `check` rightly reports that lip-block forms no surface at lip 0.)
(def bracket-field (if (> lip 0) (min upright arm corner lip-block) (min upright arm corner)))

(body bracket
  (shape bracket-field)
  (attr color "#6A9F58")
  (attr material "petg")
  (attr density 1.27g/cm3))

; The obvious alternative, (smooth-union upright arm fillet), is wrong here:
; smooth-union blends wherever the two fields are both near zero, and the two
; boxes share the whole 5 x 5 overlap block's surfaces (back face, top face,
; sides). It adds material all around that block: eval at the outer corner
; (0, 0, height) gives -1.0 instead of 0, the back face bulges out to
; y = -0.5 near the seam, and the side faces get a visible ridge
; (ink render --target bracket-smooth). Kept for comparison.
(def bracket-smooth (smooth-union upright arm fillet))

; ---------------------------------------------------------------- face references
; A face reference = material within 1 mm of a datum plane, restricted to
; where the named primitive is. Written by hand for each face.
(def top-face-region (max (- (abs (- z top-plane)) 1) arm))
(def wall-face-region (max (- (abs (- y wall-plane)) 1) upright))

; Two cheaper forms that do NOT carry intent, for comparison (material each
; selects, ink mass --res 0.25, in mm^3; the arm top is 30 x 40 x 1 = 1200):
;                          default   height=70,reach=60   lip=5
;   top-face-region         1200         1800            1200
;   top-slab                1200         1800            1646 (includes the lip)
;   top-by-bounds            849         1677              74 (the lip's top only)
;   wall-face-region        1500         2100            1500
; (1) a bare slab by coordinate, "z > height - 1": survives param changes but
;     picks up anything else at that height (here the whole 5 mm lip).
(def top-slab (- (- top-plane 1) z))
; (2) from load-time proven bounds: "the top 1 mm of whatever is highest".
;     Needs no datum, but with --param lip=5 it moves the load onto the lip.
;     `bounds` is an outer bound padded by up to a cell (here 0.29 mm), so the
;     band is 0.71 mm, not 1 mm. It does follow --param (load-time values are
;     recomputed per run).
(def top-by-bounds (- (- (nth (nth (bounds bracket-field) 1) 2) 1) z))

; Material actually selected by each reference (ink mass --target NAME):
(def sel-top-face (max bracket-field top-face-region))
(def sel-top-slab (max bracket-field top-slab))
(def sel-top-bounds (max bracket-field top-by-bounds))
(def sel-wall-face (max bracket-field wall-face-region))

; ---------------------------------------------------------------- study
(study shelf-load
  (shape bracket-field)
  (elastic (E 2100MPa) (nu 0.38) (density 1270kg/m3) (yield 45MPa))
  (fix wall wall-face-region)
  (load top (force 0 0 -100N) (region top-face-region))
  (res 1))

; Guard the references: each band must hold a plausible amount of material
; (at least 100 mm^3 of the 1200/1500 expected). There is no claim form for
; "this region is (part of) the body's surface" or "this is one planar face".
(require top-face-selects (volume sel-top-face > 100mm3))
(require wall-face-selects (volume sel-wall-face > 100mm3))
(require holds (min-safety-factor shelf-load > 2))

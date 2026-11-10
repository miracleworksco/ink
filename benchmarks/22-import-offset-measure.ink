; Benchmark 22: imported STL, redistance / offset / shell, measured by rays
; and slices.
;
; Task: a supplier sends an STL of a counterbored plate (inputs/22-part.stl,
; 984 triangles: tools/22-part-source.ink meshed by ink mesh
; tools/22-part-source.ink --res 0.2 --out inputs/22-part.stl).
; Import it, measure it, and derive: a 1 mm offset (a clearance envelope),
; a 1.5 mm shell, and the same offset done on a non-distance field with and
; without redistance.
;
; Reference values (exact part): V = 60 x 40 x 10 - pi 8^2 6 - pi 12^2 4
;   = 20 984.07 mm^3; bore r 8, counterbore r 12, plate thickness 10,
;   counterbore floor 4 below the top.
;   Offset by 1 (Minkowski with a ball): the box grows to 62 x 42 x 12 with
;   rounded edges and the holes shrink by 1:
;   V = 60 40 10 + 2 (60 40 + 60 10 + 40 10) + pi (60 + 40 + 10) + 4/3 pi
;       - pi 7^2 6 - pi 11^2 4 + (hole-edge fillets, small).
;   Simpler exact checks used instead: the offset's bounds are
;   [-31, -21, -1] .. [31, 21, 11]; the bore of the offset part fits r = 7.
;   Non-distance control: g = 3 f (same zero set, gradient 3). (- g 1)
;   moves the surface by only 1/3 mm; (- (redistance g) 1) by 1.
;
; Pass criteria:
;   - imported volume within 0.3% of 20 984.07;
;   - bore radius 8 +- 0.02 (a slice), counterbore depth 4 +- 0.02 and
;     plate thickness 10 +- 0.02 (rays), on the imported field;
;   - offset: top face at z = 11 (+- 0.03); offset bore radius 7 +- 0.03;
;   - naive offset of 3f moves the top face by 1/3 (+- 0.02), redistanced
;     by 1 (+- 0.03);
;   - the shell is 1 solid + 1 cavity, watertight.
; Distances are taken with ink ray and ink slice.

; The offset, the shell and the two offsets of 3f are separate tasks,
; side by side.

; The shell offsets samples: approximate, to within their error and
; spacing (their trilinear interpolant is no distance), as ink check says.
; expect-warning: 22-import-offset-measure.ink:48): (offset-shape f h): an offset of sampled data is approximate
; expect-warning: 22-import-offset-measure.ink:48): (offset-shape f (neg h)): an offset of sampled data is approximate

(import "std")

(def imported (import-stl "inputs/22-part.stl" (res 0.25)))
(body part (shape imported) (attr color "#8FA3B8"))

(def grown (- imported 1))
(body grown (shape (move grown 100 0 0)))
(body hollow (shape (move (shell imported 1.5) 200 0 0)))

(def g3 (* 3 imported))
(def naive (- g3 1))
(def exact (- (redistance g3 (res 0.25)) 1))
(body naive (shape (move naive 0 100 0)))
(body exact (shape (move exact 0 200 0)))



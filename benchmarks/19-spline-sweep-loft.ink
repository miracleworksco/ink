; Benchmark 19: sweeps along a 3D spline and lofts.
;
; Task (a) duct: a 16 x 8 rounded-rectangle section (r 2) swept along a
; Catmull-Rom spline through 5 points in 3D (rotation-minimising frames,
; 12 segments per span). Checks: volume = section area x path length
; (Pappus: exact for a sweep whose bend radius exceeds the section), a
; single watertight shell, flat ends.
; Task (b) lofts, z 0..50:
;   frustum: circle r 20 -> circle r 10: exact frustum,
;     V = pi h (R^2 + R r + r^2) / 3 = 36 651.9 mm^3;
;   transition: circle r 20 -> square 40 x 40: ends must be exact
;     (area 1256.6 at z = 0, 1600 at z = 50); a ruled loft (points at the
;     same polar angle joined by straight lines) has a mid-height section
;     of 1419.3 mm^2 (computed numerically).
; Task (c) taper: a section narrowing from 16 to 8 wide along the duct's
;   route (std `sweep-varying`): volume = mean section area x path length.
;
; Hand values for (a): section area 16 x 8 - (4 - pi) 2^2 = 124.566 mm^2;
; path length from the polyline itself (`path-length` below, reported by
; ink check: 199.8046).
;
; Pass criteria (ink mass --res 0.25, ink mesh, ink slice):
;   - duct volume within 1% of area x length; watertight, one shell;
;   - frustum volume within 0.5% of 36 651.9;
;   - transition: end sections exact within 1%; mid section within 3% of
;     the ruled 1419.3;
;   - taper volume within 1% of mean area x length.

; The duct and the three lofts are separate tasks, side by side.

(import "std")

(def route (spline-points (vec (vec 0 0 0) (vec 40 10 20) (vec 70 50 30) (vec 110 60 10) (vec 150 40 40)) 12))
(def section (rounded-rect-2d 16 8 2))
(body duct (shape (sweep section route (vec 0 0 1))) (attr color "#6C8EBF"))
(def path-length (reduce + (map (lambda (a b) (norm (- b a))) (slice route 0 -1) (slice route 1))))
(def section-area (- 128 (* (- 4 pi) 4)))
(def duct-volume-expected (* section-area path-length))

(body frustum (shape (move (loft-z (circle-2d 20) (circle-2d 10) 0 50) 0 -80 0)) (attr color "#D79B00"))
(def transition-field (loft-z (circle-2d 20) (rect-2d 40 40) 0 50))
(body transition (shape (move transition-field 0 -150 0)) (attr color "#82B366"))
; the frustum within 0.5% of pi h (R^2 + R r + r^2) / 3 = 36 651.9 (measured)
(assert frustum-lo (volume (shape-of frustum) > 36468.6))
(assert frustum-hi (volume (shape-of frustum) < 36835.2))

; (c) a tapering duct on the same route: the section narrows from 16 to 8 wide
; along the length. Area is linear in s, so V = mean area x length.
(body taper (shape (move (sweep-varying (lambda (s) (rounded-rect-2d (mix 16 8 s) 8 2)) route (vec 0 0 1)) 0 120 0))
  (attr color "#B85450"))
(def taper-volume-expected (* (- 96 (* (- 4 pi) 4)) path-length))

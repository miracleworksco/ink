; Benchmark 01: swept pipe with end flanges.
;
; Task: a tube OD 20, wall 2, swept along a 3D cubic Bezier curve, with a
; flange at each end (OD 50, 6 thick) carrying 4 bolt holes (6.6 mm, PCD 38).
; The path is a 32-segment polyline (std `bezier3-points`), within 0.09 mm
; of the curve.
;
; Pass criteria (ink assert --res 0.1):
;   - the bore is clear: a 7.5 mm radius gauge along the path is disjoint
;     from the part (proven);
;   - the bolt holes go through: 3 mm pins through each flange are disjoint
;     from the part (proven, both ends).

(import "std")

(param od 20)               ; tube outer diameter
(param wall 2)              ; tube wall
(param fl-od 50)            ; flange outer diameter
(param fl-t 6)              ; flange thickness
(param bolt-d 6.6)          ; bolt hole diameter (M6 clearance)
(param pcd 38)              ; bolt pitch-circle diameter
(def nseg 32)               ; path segments (load-time number, not a param)

; Cubic Bezier control points: length 212.03 mm, min bend radius 46.5 mm.
(def route (bezier3-points (vec 0 0 0) (vec 90 0 0) (vec 40 110 30) (vec 130 110 90) nseg))

; End frames: origin at the end, direction into the pipe.
(def end-a (nth route 0))
(def dir-a (unit (- (nth route 1) (nth route 0))))
(def end-b (nth route -1))
(def dir-b (unit (- (nth route -2) (nth route -1))))

; Extrude a 2D profile from s0 to s1 along direction dir from point p.
(fn plate-at (profile p dir s0 s1) (sweep profile (vec (+ p (* dir s0)) (+ p (* dir s1))) (vec 0 0 1)))

(def bolt-holes-2d (pattern-polar (move (circle-2d (* 0.5 bolt-d)) (* 0.5 pcd) 0 0) 4))
(def flange-2d (subtract (circle-2d (* 0.5 fl-od)) bolt-holes-2d))
(def flanges (min (plate-at flange-2d end-a dir-a 0 fl-t) (plate-at flange-2d end-b dir-b 0 fl-t)))

(def bore-r (- (* 0.5 od) wall))
(def pipe-field (subtract (min (pipe route (* 0.5 od) wall) flanges) (tube route bore-r)))

(body swept-pipe
  (shape pipe-field)
  (attr color "#B87333")
  (attr material "copper")
  (attr density 8.96g/cm3))

; ---------------------------------------------------------------- checks
; The gaps are 0.3-0.5 mm; these hold to 0.05 at the automatic res.
; The bore is clear: a 7.5 mm radius gauge along the path never touches the part.
(assert bore-clear (disjoint (shape-of swept-pipe) (tube route 7.5)))
; Bolt holes go through: 3 mm pins through them (from 5 mm outside to 11 mm in)
; touch nothing.
(def pins-2d (pattern-polar (move (circle-2d 3) (* 0.5 pcd) 0 0) 4))
(assert holes-clear-a (disjoint (shape-of swept-pipe) (plate-at pins-2d end-a dir-a -5 11)))
(assert holes-clear-b (disjoint (shape-of swept-pipe) (plate-at pins-2d end-b dir-b -5 11)))
; Volume. Analytic for the true curve: tube pi (10^2 - 8^2) * 212.034 = 23981
; plus two flange rings outside the tube, 2 * 6 * (pi (25^2 - 10^2) - 4 pi 3.3^2)
; = 2 * 6 * (1649.3 - 136.8) = 18150, total 42131 mm^3.
; `ink mass --res 0.1` estimates 42124 (-0.02%). The volume claims are
; estimates too, with no guaranteed accuracy.
(assert pipe-volume-lo (volume (shape-of swept-pipe) > 38000))
(assert pipe-volume-hi (volume (shape-of swept-pipe) < 46500))

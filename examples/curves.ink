; Curves, profiles and sweeps: three shapes built from data with std, side
; by side.

(import "std")

; A pipe along a smooth spline through four points, with bolted flanges at both ends.
(def route (spline-points (vec (vec 0 0 0) (vec 60 10 20) (vec 90 60 30) (vec 140 70 10)) 12))
(def flange-2d (subtract (circle-2d 22) (pattern-polar (move (circle-2d 3) 16 0 0) 4)))

; Place a 2D profile as a plate of thickness th at path point p facing along dir.
(fn plate-at (profile p dir th)
  (sweep profile (vec (- p (* dir th)) p) (vec 0 0 1)))

(def route-pipe (pipe route 10 2))
(def flanges
  (min (plate-at flange-2d (nth route 0) (unit (- (nth route 1) (nth route 0))) -4)
       (plate-at flange-2d (nth route -1) (unit (- (nth route -1) (nth route -2))) 4)))
(body piping (shape (subtract (min route-pipe flanges) (tube route 8))) (attr color "#6C8EBF"))
; Its bore is open end to end: a 12 mm rod along the route touches nothing.
(assert bore-open (disjoint (shape-of piping) (tube route 6)))

; A rectangular bar swept along a helix (the profile stays square to the path).
(body coil (shape (move (sweep (rounded-rect-2d 6 3 0.8) (helix-points 25 15 2 120) (vec 0 0 1)) 0 -90 0))
  (attr color "#D79B00"))

; A notched star plate: an arbitrary non-convex polygon, extruded.
(def star (map (lambda (i) (let ((r (if (== (- i (* 2 (floor (/ i 2)))) 0) 30 14))
                                (a (/ (* 2 pi i) 12)))
                            (vec (* r (cos a)) (* r (sin a)))))
               (range 12)))
(body star-plate (shape (move (extrude-z (polygon star) 0 6) 100 -90 0)) (attr color "#82B366"))

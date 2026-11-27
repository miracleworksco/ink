; An outside read's case: tests/data/rated-beam.ink's bar under two studies, pressed with
; 1000 N and with 500 N, its yield rated 0.30 GPa. They are its two load cases.
(import "std")
(def rated (read-json "two-loads.json"))
(body bar (shape (move (box 210 10 10) 100 0 0)) (attr ratings (get rated "ratings")))

(fn edge (x0) (- (len2 (- x x0) (+ z 5)) 1))

(study press
  (shape (shape-of bar))
  (elastic (material aluminum-6061))
  (fix pin (edge 0) (pivot (at 0 0 -5) (axis 0 1 0)))
  (fix roller (edge 200) (dofs z) (pivot (at 200 0 -5) (axis 0 1 0)))
  (load hand (force 0 0 -1000N) (region (- (abs (- x 100)) 1)))
  (res 0.5))

(study fallback
  (shape (shape-of bar))
  (elastic (material aluminum-6061))
  (fix pin (edge 0) (pivot (at 0 0 -5) (axis 0 1 0)))
  (fix roller (edge 200) (dofs z) (pivot (at 200 0 -5) (axis 0 1 0)))
  (load hand (force 0 0 -500N) (region (- (abs (- x 100)) 1)))
  (res 0.5))

(assembly root
  (instance bar bar (at)))
(require strong-and-stiff (within-ratings root))

; examples/beam.ink, with the bar carrying its limits as a catalog part
; carries ratings (the kernel's ink-ratings: its yield, 0.24 GPa, and the
; 0.4 cm it may sag), placed in an assembly, and the design asking whether
; it stays within them. Nothing states its stress or its sag: the solver
; gives them. The study is beam.ink's, to the letter.
(import "std")
(def rated (read-json "rated-beam.json"))
(body bar (shape (move (box 210 10 10) 100 0 0)) (attr ratings (get rated "ratings")))

(fn edge (x0) (- (len2 (- x x0) (+ z 5)) 1))

(study press
  (shape (shape-of bar))
  (elastic (material aluminum-6061))
  (fix pin (edge 0) (pivot (at 0 0 -5) (axis 0 1 0)))
  (fix roller (edge 200) (dofs z) (pivot (at 200 0 -5) (axis 0 1 0)))
  (load hand (force 0 0 -1000N) (region (- (abs (- x 100)) 1)))
  (res 0.5))

(assembly root
  (instance bar bar (at)))
(require strong-and-stiff (within-ratings root))

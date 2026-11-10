; offset is a name a model may define; the built-in, offset-shape, keeps
; its meaning wherever it is read, passed as a value too.
(import "std")
(fn offset (f d) (- f (* 2 d)))
(def grow offset-shape)
(body ball (shape (grow (sphere 1) 1)))
(assert grown-by-one (disjoint (probe 2.5 0 0) (shape-of ball)))
(assert grown (subset (probe 1.9 0 0) (shape-of ball)))
(body plate (shape (move (shell (box 10 10 2) 1) 4 0 0)))
(assert hollow (disjoint (probe 4 0 0) (shape-of plate)))

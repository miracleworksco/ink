; A 12-tooth gear spinning up from rest as t^2, through 30 turns, and a
; shaft turning with it 2 mm past its tips. The pair reads t only through
; the drive's angle, and the gear is the same turned by a tooth, the shaft
; by any turn: through that angle the pair repeats every 30 degrees, so
; one tooth's turn, searched, stands for the 360 that pass the shaft. The
; gear's field folds at its pitch, so how fast it changes is not bounded
; where the turn crosses a fold: the search halves the angle first there.
(import "std")
(t-span 0 1)
(def drive (* 10800 (square t)))
(body gear (shape (extrude-z (involute-gear-2d 12 2 20 0) 0 5)))
(body shaft (shape (extrude-z (circle-2d 4) -5 10)))
(assembly root
  (instance g gear (at) (turn z drive))
  (instance s shaft (at (move 20 0 0)) (turn z drive)))
(assert clear (disjoint (shape-of root/s) (shape-of root/g)))

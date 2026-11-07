; A 10:20 gear pair spinning up from rest as t^2, the wheel at mesh-turn of
; the pinion's angle: the pair repeats every pinion tooth (the wheel turns
; a whole tooth with it), so one tooth, searched, stands for the 150 that
; pass. In phase, their flanks pass 0.1 mm apart: proven. Half a wheel
; tooth out of phase, the teeth clash at a moment.
(import "std")
(t-span 0 1)
(def drive (* 5400 (square t)))
(fn gear (n) (extrude-z (involute-gear-2d n 2 20 0.1) 0 5))
(body pinion (shape (gear 10)))
(body wheel (shape (gear 20)))
(assembly root
  (instance p pinion (at) (turn z drive))
  (instance w wheel (at (move (mesh-distance 10 20 2) 0 0)) (turn z (mesh-turn 10 20 0 drive)))
  (instance q pinion (at (move 0 0 20)) (turn z drive))
  (instance v wheel (at (move (mesh-distance 10 20 2) 0 20)) (turn z (+ 9 (mesh-turn 10 20 0 drive)))))
(assert meshes (disjoint (shape-of root/p) (shape-of root/w)))
(assert clashes (disjoint (shape-of root/q) (shape-of root/v)))
; expect-status: clashes fails

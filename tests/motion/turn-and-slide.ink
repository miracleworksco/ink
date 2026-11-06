; An arm turned once a second and a post slid toward it at 2 mm/s, both
; driven by t: the pair reads t through the slide too, so it has no
; period. The post reaches the arm's tip after t = 9.75; a cut to one
; turn would see it at 58..60 only.
(import "std")
(t-span 0 10.2)
(body arm (shape (move (box 30 2 2) 25 0 0)))
(body post (shape (cylinder-z 0.5 4)))
(assembly root
  (instance a arm (at) (turn z (* 360 t)))
  (instance p post (at (move 60 0 0)) (slide x (* -2 t))))
(assert hits (disjoint (shape-of root/a) (shape-of root/p)))
; expect-status: hits fails

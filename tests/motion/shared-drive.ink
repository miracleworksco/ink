; Two lines turned by one drive, as rollers on one belt are written: the
; pair of an arm on one line and a post reads both lines' turns (one node),
; and nothing but the sine and cosine of that one angle, so it repeats
; every turn of the drive: a post out of reach is proven clear over one
; turn, and so over the ten. A post in its reach at 90 degrees is hit at a
; moment it gets there.
(import "std")
(t-span 0 1)
(body arm (shape (box 30 2 2)))
(body post (shape (cylinder-z 0.5 4)))
(assembly line (instance arm arm (at (move 25 0 0))))
(assembly root
  (instance a line (at) (turn z (* 3600 t)))
  (instance b line (at (move 0 0 20)) (turn z (* 3600 t)))
  (instance p post (at (move 41 0 0)))
  (instance q post (at (move 0 35 0))))
(assert clear (disjoint (shape-of root/a/arm) (shape-of root/p)))
(assert hits (disjoint (shape-of root/a/arm) (shape-of root/q)))
; expect-status: hits fails
; expect-status: clear holds

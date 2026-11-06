; A roller as a CAD converter writes one (its tube round about its own z at
; (-50.4, 16.1), placed by a float turn), placed on its line by six digits
; that leave its axis 3e-4 mm off the line's and tilted about 5e-7 rad,
; spinning ten turns on an arm that swings 60 degrees past a post 0.5 mm
; from it (at 30 degrees, 113 from the arm's axis). Its axis is a float
; turn off its spin's, so no exact turn about its own axis places it, and
; t ranged leaves the clear pairs undecided: held at the moments measured,
; estimates. A post in its sweep fails. The same
; roller 0.02 mm off, the same.
(import "std")
(t-span 0 1)
(body roller (shape (place (subtract (move (cylinder-z 12 30) -50.4 16.1 0) (move (cylinder-z 10 31) -50.4 16.1 0)) (at (rot 0.57735 0.57735 0.57735 120)))))
(body arm (shape (move (box 100 6 4) 50 0 -20)))
(body post (shape (cylinder-z 0.5 40)))
(assembly line (instance tube roller (at (move 0 50.4003 -16.1) (rot 0.57735 0.57735 0.577351 240))))
(assembly far-line (instance tube roller (at (move 0 50.42 -16.1) (rot 0.57735 0.57735 0.577351 240))))
(assembly swing (instance a arm) (instance l line (at (move 100 0 0)) (turn z (* 3600 t))))
(assembly far-swing (instance a arm) (instance l far-line (at (move 100 0 0)) (turn z (* 3600 t))))
(assembly root
  (instance s swing (at) (turn z (* 60 t)))
  (instance f far-swing (at (move 0 0 100)) (turn z (* 60 t)))
  (instance clear post (at (move 97.860871 56.5 0)))
  (instance hit post (at (move 96.12882 55.5 0)))
  (instance far-clear post (at (move 97.860871 56.5 100))))
(assert clear (disjoint (shape-of root/s/l/tube) (shape-of root/clear)) (accuracy 0.1))
(assert hits (disjoint (shape-of root/s/l/tube) (shape-of root/hit)) (accuracy 0.1))
(assert far (disjoint (shape-of root/f/l/tube) (shape-of root/far-clear)) (accuracy 0.1))
; expect-status: clear holds estimate
; expect-status: hits fails
; expect-status: far holds estimate

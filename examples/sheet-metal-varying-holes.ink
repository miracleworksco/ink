; Holes that vary, in std's sheet. sheet-collisions pairs pieces by the
; boxes of the sheet without its holes, so a hole placed by a toleranced
; param, or one that moves, is checked over its range. Walls rise 90 at
; x = 10 and 50, with lips 30 past their zones. Uncut, the left lip's tip
; reaches 0.5 into the right lip's zone (examples/sheet-metal.ink,
; lip-meets-zone). A cut 3 wide over the tip's flat end shortens it.
; expect-status: tip-moving-cut fails
(import "std")

(def th 2)
(def k 0.4)
(def r 2)
(t-span 0 1)
(fn lips (cut)
  (sheet (polygon (vec (vec -58.6009 0) (vec 100 0) (vec 100 30) (vec -58.6009 30)))
         (vec (move (polygon (vec (vec -60 -1) (vec -57 -1) (vec -57 31) (vec -60 31))) cut 0 0)) th k
         (vec (sheet-bend (vec 10 0) (vec 10 30) 90 r) (sheet-bend (vec 50 30) (vec 50 0) 90 r)
              (sheet-bend (vec -18.3 0) (vec -18.3 30) 90 r) (sheet-bend (vec 80 30) (vec 80 0) -90 r))))
; the cut placed within 1 of its place: it takes 1.6 to 2.6 off the tip,
; which stays clear of the zone over the whole band
(param cut-at 0 (tol 0 1))
(assert tip-toleranced-cut (empty (sheet-collisions (lips cut-at))))
; the cut moving 3 over the motion: at t = 1 it leaves the tip whole
(assert tip-moving-cut (empty (sheet-collisions (lips (* 3 t)))))

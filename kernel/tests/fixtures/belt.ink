; A circulating band with rigid rectangular teeth, using only std paths.
(import "std")
(fn track (r c)
  (path-loop (vec
    (path-arc (vec 0 0 0) (vec 0 0 1) r (vec 0 1 0) pi)
    (path-line (vec 0 (neg r) 0) (vec c (neg r) 0))
    (path-arc (vec c 0 0) (vec 0 0 1) r (vec 0 -1 0) pi)
    (path-line (vec c r 0) (vec 0 r 0))) (vec 0 0 1)))
(fn belt (r c shift)
  (let ((P (track r c)))
    (min (path-sweep (move (rect-2d .6 6) -.3 0 0) P)
      (max (path-repeat (move (rect-2d 1 1.5) 0 .5 0) P 2 shift) (- (abs z) 3)))))

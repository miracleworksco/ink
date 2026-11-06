; Chordal path placements and two link bodies, independent of chain catalogs.
(import "std")
(def pitch 6.35)
(def r (/ pitch (* 2 (sin (/ pi 16)))))
(def rho (/ (* 16 pitch) (* 2 pi)))
(def center (* 8 pitch))
(def P (path-loop (vec
  (path-arc (vec 0 0 0) (vec 0 0 1) r (vec 0 1 0) pi rho)
  (path-line (vec 0 (neg r) 0) (vec center (neg r) 0))
  (path-arc (vec center 0 0) (vec 0 0 1) r (vec 0 -1 0) pi rho)
  (path-line (vec center r 0) (vec 0 r 0))) (vec 0 0 1)))
(t-span 0 2)
(body inner-link (shape (box 5 2 2)))
(body pin-link (shape (box 5 3 2)))
(body wheel (shape (cylinder-z 14 2)))
(assembly root
  (instance a wheel) (instance b wheel (at (move center 0 0)))
  (for i (range 0 32 2) (instance link-$i inner-link (path-link-at P (+ (* i pitch) (* 20 t)) pitch)))
  (for i (range 1 32 2) (instance link-$i pin-link (path-link-at P (+ (* i pitch) (* 20 t)) pitch))))

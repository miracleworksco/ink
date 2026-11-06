; Jumps across t and p, both ranged in one search. An arm jumping from 0
; to 37.3 degrees at t = p / 100 under a ball 1.1 above its path at 37.3:
; the jump is a line through the boxes of t and p (stepping followed it
; without end: 11678 boxes, none proven, in 900 s); it holds.
; A second arm is at 37.3 only in a band 0.002 wide past t = 0.3 + p / 100
; (two parallel jump lines), over a ball 0.97 above: it fails in the band.
; A third is at 37.3 only where two crossing bands meet (a window about 10
; ms by 1 of p near t 0.4, p 40): over a ball 1.1 above it holds; over one
; 0.97 above it fails there (t and p ranged in one search).
(import "std")
(t-span 0 1)
(param p 45 (tol -45 45))
(body arm (shape (move (cylinder-x 2 40) 25 0 0)))
(body ball (shape (sphere 3)))
(def bx (* 30 (cos (deg->rad 37.3))))
(def by (* 30 (sin (deg->rad 37.3))))
(def band (lambda (from w) (select (- t from) 0 (select (- t (+ from w)) 1 0))))
(assembly root
  (instance arm arm (at (move 0 0 0)) (turn z (select (- t (/ p 100)) 0 37.3)))
  (instance ball ball (at (move bx by 6.1)))
  (instance band-arm arm (at (move 0 100 0)) (turn z (* 37.3 (band (+ 0.3 (/ p 100)) 0.002))))
  (instance band-ball ball (at (move bx (+ 100 by) 5.97)))
  (instance cross-arm arm (at (move 0 200 0)) (turn z (* 37.3 (band (/ p 100) 0.01) (band (- 0.8 (/ p 100)) 0.01))))
  (instance cross-ball ball (at (move bx (+ 200 by) 6.1)))
  (instance cross-low ball (at (move (* 0.5 bx) (+ 200 (* 0.5 by)) 5.97))))
(assert clear (clearance (shape-of root/arm) (shape-of root/ball) 1))
(assert clear-band (clearance (shape-of root/band-arm) (shape-of root/band-ball) 1))
; expect-status: clear-band fails
(assert clear-cross (clearance (shape-of root/cross-arm) (shape-of root/cross-ball) 1))
(assert clear-cross-low (clearance (shape-of root/cross-arm) (shape-of root/cross-low) 1))
; expect-status: clear-cross-low fails

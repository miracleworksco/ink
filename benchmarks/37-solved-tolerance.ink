; Benchmark 37: a part placed by a solve, proven over its tolerances.
;
; Task: a link of length L = 50 +-0.1 is pinned at the origin and its far
; end rides on a rail at height h = 30 +-0.1. The slider (a 10 mm block)
; on that end sits where the constraints put it: the link's angle theta
; (degrees) and the slider's x, sx, are solved, L cos theta = sx and
; L sin theta = h, not written in closed form. A stop at x = 45.5 must clear the slider
; by 0.25 at every L and h in their bands; a tighter stop at x = 45.1 does
; not, though it clears at the nominal sizes.
;
; Hand calculation (worst case):
;   - x = sqrt(L^2 - h^2), rising with L and falling with h:
;     x in [sqrt(49.9^2 - 30.1^2), sqrt(50.1^2 - 29.9^2)]
;       = [sqrt(1584), sqrt(1616)] = [39.7995, 40.1995] (nominal 40);
;   - theta = asin(h / L) in [asin(29.9 / 50.1), asin(30.1 / 49.9)]
;       = [36.643, 37.100] deg (nominal 36.870);
;   - the slider's right face is at x + 5 <= 45.1995: the stop at 45.5
;     leaves at least 0.3005 (a 0.25 clearance holds); the tight stop at
;     45.1 is 0.1 clear at the nominal sizes and 0.0995 into the slider at
;     L = 50.1, h = 29.9 (it fails wherever L^2 - h^2 > 40.1^2 = 1608.01).
;
; Pass criteria (ink assert, bands applied by default):
;   - slider-clears-stop, x-high, x-low, theta-high, theta-low: hold to
;     their accuracy (not estimates: over the bands the solved values are
;     certified roots, docs/IR.md, Solving constraints);
;   - x-under and tight-apart fail (deliberate), each witness with L and h
;     in their bands and L^2 - h^2 above the hand value where it
;     fails (x-under: x > 40.15, so L^2 - h^2 > 1612.02; tight-apart: >
;     1608.01);
;   - with --nominal (L and h at their values) tight-apart holds
;     numerically: the tolerances are what fails it.

; The stop stands clear of the slider by design (slider-clears-stop). The
; rail and link holding both are not modelled.

(import "std")

(param L 50 (tol -0.1 0.1))        ; the link's length
(param h 30 (tol -0.1 0.1))        ; the rail's height

; Where the constraints put the link and the slider.
(solve place
  (unknowns (theta 37) (sx 40))
  (== (* L (cos (* theta (/ pi 180)))) sx)
  (== (* L (sin (* theta (/ pi 180)))) h))

(body slider (shape (box 10 10 10)))
(body stop (shape (box 2 20 20)))

(assembly root
  (instance slider slider (at (move sx h 0)))
  (instance stop stop (at (move 46.5 30 0))))

(def tight-stop (move (box 2 20 20) 46.1 30 0))

(assert slider-clears-stop (clearance (shape-of root/slider) (shape-of root/stop) 0.25))
(assert x-high (value sx <= 40.2))
(assert x-low (value sx >= 39.79))
(assert theta-high (value theta <= 37.11))
(assert theta-low (value theta >= 36.63))

; The tight stop overlaps the slider by up to 0.0995 (long link, low rail).
; The claim fails at accuracy 0.01 and is unknown at the default 0.05.
; expect-status: x-under fails
(assert x-under (value sx < 40.15))
; expect-status: tight-apart fails
(assert tight-apart (disjoint (shape-of root/slider) tight-stop) (accuracy 0.01))

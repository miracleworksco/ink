; Claims of 2D profiles (fields of x and y alone, the same at every
; height) are asked in the plane, and answered as of the prism each is.
; expect-status: inside-too-small fails
; expect-status: touching-gap fails
; expect-status: plate-overlaps-hole fails
; expect-status: fin fails estimate
; expect-status: walls holds estimate
(import "std")
(def plate (subtract (rect-2d 40 20) (circle-2d 4)))
(def fin (move (rect-2d 0.6 10) 0 14 0))   ; standing off the plate's edge
(assert inside (subset (move (circle-2d 3) -12 0 0) plate))
(assert inside-too-small (subset plate (rect-2d 39 20)))
(assert apart (clearance plate (move (circle-2d 2) 30 0 0) 7.9))
(assert touching-gap (clearance plate (move (circle-2d 2) 30 0 0) 8.1))
(assert hole-is-air (disjoint plate (circle-2d 3.9)))
(assert plate-overlaps-hole (disjoint plate (circle-2d 5)))
(assert walls (min-thickness plate > 5.9))
(assert fin (min-thickness (min plate fin) > 1))
(assert no-voids (cavities plate < 1))

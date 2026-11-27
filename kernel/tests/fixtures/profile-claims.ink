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
(require inside (subset (move (circle-2d 3) -12 0 0) plate))
(require inside-too-small (subset plate (rect-2d 39 20)))
(require apart (clearance plate (move (circle-2d 2) 30 0 0) 7.9))
(require touching-gap (clearance plate (move (circle-2d 2) 30 0 0) 8.1))
(require hole-is-air (disjoint plate (circle-2d 3.9)))
(require plate-overlaps-hole (disjoint plate (circle-2d 5)))
; walls on cells of 0.16 mm: on 1.25, 0.63 and 0.31 they read 6.049, 6.012 and
; 6.003, converging faster than their elements can, a grid error of 0.14 past
; the 0.1 margin; on 0.64, 0.32 and 0.16, 6.002, 6, 6.0005 +- 0.007 (exact 6).
(require walls (min-thickness plate > 5.9) (res 0.16))
(require fin (min-thickness (min plate fin) > 1))
(require no-voids (cavities plate < 1))

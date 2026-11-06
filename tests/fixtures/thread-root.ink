; The atan2 cut at an internal ISO thread root: a point there.
(import "std")
(body nut (shape (subtract (cylinder-z 4 6) (thread-iso 4.08 .7))))
(assert thread-minor-metal (subset (probe 0 1.693613886 1.575) (shape-of nut)))
; A DIN 1587 M4 cap nut (its standard's sizes: s 6.89, h 7.82, m 3.05,
; dk 6.5, t 5.5; tapped M4x0.7, 0.04 radial allowance), a point at its root.
(body din
  (shape (subtract (min (max (hexagon-2d 6.89) (slab-z 0 3.05))
                        (min (extrude-z (circle-2d 3.25) 2.85 4.57) (max (move (sphere 3.25) 0 0 4.57) (- 2.85 z))))
                   (max (thread-iso 4.08 .7) (- z 5.5)))))
(assert din-root (subset (probe 0 1.662113886 1.575) (shape-of din)))

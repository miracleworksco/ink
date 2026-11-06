; An FDM-printed cable guide, checked against Protolabs' FDM design rules
; (docs/IR.md, Manufacturing checks), printed upright on a bed at z = 0:
; walls at least 0.8 mm, overhangs to 45 degrees past vertical, and 0.5 mm
; between parts printed in place that must move.
;
; A base and an upright wall, with a teardrop bolt hole and a pointed cable
; window, their roofs 40 degrees past vertical. Above them a shelf stands
; out 5 mm on a 40-degree chamfer. Nothing needs support or a bridge. Each
; failing control changes one thing: a round hole (its crown faces straight
; down), a shelf with a flat underside, a 0.5 mm fin, and a print-in-place
; hinge pin 0.3 mm from its knuckle. Bridges are not checked: a bridge
; reads as an overhang.
;
; By hand: each roof is two planes 50 degrees from horizontal, tangent to
; the hole or rising from the window's corners, so 40 degrees past
; vertical, 5 under the limit. The apex is an upward edge. At step 1 and
; the default accuracy 0.05 (how far off an answer may be), a pass means no
; face leans past about 49 degrees (asin(2 0.05 cos 45 / 1) = 4.1 past 45;
; docs/IR.md, Manufacturing checks). A face under 42.7 never fails, so the
; roofs' 40 is decided.
(import "std")

(param limit 45)           ; the steepest overhang, degrees past vertical
(param step 1)             ; how coarsely the overhang check reads the part, mm
(param roof 50)            ; the hole's and window's roof, degrees from horizontal

(def base (move (box 40 20 4) 0 0 2))
(def wall (move (box 40 4 30) 0 0 15))
(def rs (sin (deg->rad roof)))
(def rc (cos (deg->rad roof)))
; Along y through the wall: the hole's axis at x = 10, z = 20, radius 4, a
; roof of two planes tangent to it meeting 4 / cos(roof) above the axis;
; the window 8 wide, z 9 to 15, its roof two planes rising from its top
; corners to meet over its middle.
(def hx (- x 10))
(def hz (- z 20))
(def round-hole (- (len2 hx hz) 4))
(def teardrop (max (min round-hole (max (+ (* rs hx) (* rc (- hz (/ 4 rc)))) (+ (* (neg rs) hx) (* rc (- hz (/ 4 rc))))
                                         (- (* 4 rc) hz)))
                   (- (abs y) 5)))
(def window (max (slab-x -14 -6) (- (abs y) 5) (- 9 z)
                 (- (* rc (- z 15)) (* rs (+ x 14))) (+ (* rc (- z 15)) (* rs (+ x 6)))))
; A shelf on the wall's +y face, its underside rising at the roof's angle
; out from z = 22 at the wall.
(def shelf (max (move (box 16 5 9) -10 4.5 25.5) (- (* rs (- y 2)) (* rc (- z 22)))))
(def flat-shelf (move (box 16 5 3) -10 4.5 28.5))
(fn guide (h s) (subtract (min base wall s) h window))

(def ok (guide teardrop shelf))
(def round (guide (max round-hole (- (abs y) 5)) shelf))
(def flat (guide teardrop flat-shelf))
(def finned (min ok (move (box 0.5 8 10) 19 -6 9)))
(body guide (shape ok) (attr material "petg"))
(body guide-round-hole (shape round) (attr material "petg"))
(body guide-flat-shelf (shape flat) (attr material "petg"))
(body guide-finned (shape finned) (attr material "petg"))

; A print-in-place hinge: a pin 3 mm across with a head, through a
; knuckle 10 mm across, both standing on the bed, the gap between them
; radial.
(fn hinge (gap) (min (extrude-z (circle-2d 1.5) 0 12) (extrude-z (circle-2d 4) 12 14)
                       (extrude-z (ring-2d (+ 1.5 gap) 5) 0 10)))
(def loose (hinge 0.6))
(def tight (hinge 0.3))
(body hinge (shape loose) (attr material "petg"))
(body hinge-tight (shape tight) (attr material "petg"))

(assert overhangs-ok (empty (overhangs ok 0 limit step)))
; expect-status: round-hole-overhangs fails
(assert round-hole-overhangs (empty (overhangs round 0 limit step)))
; expect-status: flat-shelf-overhangs fails
(assert flat-shelf-overhangs (empty (overhangs flat 0 limit step)))
(assert walls (min-thickness ok > 0.8))
; expect-status: fin-walls fails estimate
(assert fin-walls (min-thickness finned > 0.8))
(assert hinge-moves (min-thickness (air-around loose 2) > 0.5))
; expect-status: tight-hinge-moves fails estimate
(assert tight-hinge-moves (min-thickness (air-around tight 2) > 0.5))

; A bar leaning 60 degrees, round in section, where a slope's bounds do
; not narrow across a branch: its underside and end lean at most 30
; past vertical. Read by its slope the claim stays unknown at every res
; along the line where the underside meets the end cap; read by signs it
; holds.
(def bar (max (extrude-x (- (segment-distance-at here2 (vec 0 -5) (vec -20 30)) 3) -20 20) (neg z)))
(body bar (shape bar) (attr material "petg"))
(assert bar-overhangs (empty (overhangs bar 0 limit step)))

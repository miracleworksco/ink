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
; hinge pin 0.3 mm from its knuckle. Bridges are checked apart, below.
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

; A cable tunnel: a block with a channel through it along y, 8 mm wide,
; its flat roof a bridge from wall to wall. FDM bridges 10 mm (docs/IR.md,
; Manufacturing checks); walls under 0.8 mm fail the walls check, so the
; bridge check reads the run for walls at steps of 0.8. The roof fails as
; an overhang and holds as a bridge. Each failing control changes one
; thing: a channel 14 mm wide; the roof cut back 4 mm from one wall, a
; cantilever with walls 8 mm apart either side; the roof split by a 0.3 mm
; slit, two cantilevers that nearly meet; and the guide's flat shelf, with
; a wall on one side only.
(param span 10)            ; the longest bridge, mm
(fn tunnel (w) (subtract (move (box (+ w 8) 12 8) 0 0 4) (move (box w 20 5) 0 0 2.5)))
(def bridged (tunnel 8))
(def wide (tunnel 14))
(def cut-back (subtract bridged (move (box 4 20 4) 2 0 7)))
(def split (subtract bridged (move (box 0.3 20 4) 1 0 7)))
(body tunnel (shape bridged) (attr material "petg"))
(body tunnel-wide (shape wide) (attr material "petg"))
(body tunnel-cut-back (shape cut-back) (attr material "petg"))
(body tunnel-split (shape split) (attr material "petg"))
; expect-status: tunnel-overhangs fails
(assert tunnel-overhangs (empty (overhangs bridged 0 limit step)))
(assert tunnel-bridged (empty (unbridged bridged 0 limit step span 0.8)))
; expect-status: wide-tunnel-bridged fails
(assert wide-tunnel-bridged (empty (unbridged wide 0 limit step span 0.8)))
; expect-status: cut-back-tunnel-bridged fails
(assert cut-back-tunnel-bridged (empty (unbridged cut-back 0 limit step span 0.8)))
; expect-status: split-tunnel-bridged fails
(assert split-tunnel-bridged (empty (unbridged split 0 limit step span 0.8)))
; The roof cut by a slit 2 mm wide where the field reads 0, not above it:
; air all the same (a shape is where its field is below 0), so it fails.
(def zero-slit (max bridged (min 0 (- 1 (abs x)))))
(body tunnel-zero-slit (shape zero-slit) (attr material "petg"))
; expect-status: zero-slit-tunnel-bridged fails
(assert zero-slit-tunnel-bridged (empty (unbridged zero-slit 0 limit step span 0.8)))
; expect-status: flat-shelf-bridged fails
(assert flat-shelf-bridged (empty (unbridged flat 0 limit step span 0.8)))
(assert guide-bridged (empty (unbridged ok 0 limit step span 0.8)))

; Posts 2 mm thick and 8 mm apart under 0.8 mm roofs. A roof split by a
; 0.4 mm slit under a sound roof 0.2 mm above it fails: the lower roof's
; own bridge is broken. The sound roof alone holds. A roof on posts
; stopping 0.2 mm under it fails: no post touches it.
(def posts (min (move (box 2 12 7) -5 0 3.5) (move (box 2 12 7) 5 0 3.5)))
(def slit-roof (subtract (move (box 12 12 0.8) 0 0 5.4) (move (box 0.4 20 2) 0 0 5.4)))
(def top-roof (move (box 12 12 0.8) 0 0 6.4))
(def covered (min posts slit-roof top-roof))
(def roofed (min posts top-roof))
(def floating (min (move (box 2 12 4.8) -5 0 2.4) (move (box 2 12 4.8) 5 0 2.4) (move (box 12 12 0.8) 0 0 5.4)))
(body covered-slit (shape covered) (attr material "petg"))
(body roofed (shape roofed) (attr material "petg"))
(body floating-roof (shape floating) (attr material "petg"))
; expect-status: covered-slit-bridged fails
(assert covered-slit-bridged (empty (unbridged covered 0 limit step span 0.8)))
(assert roofed-bridged (empty (unbridged roofed 0 limit step span 0.8)))
; expect-status: floating-roof-bridged fails
(assert floating-roof-bridged (empty (unbridged floating 0 limit step span 0.8)))

; A ceiling over a gap is a bridge however close what is under it, once
; the gap is wider than 0.3 mm, where a support's two points stop reading
; it as contact: a slot 20 mm wide and 0.4 mm tall through a block, a 30
; mm plate 0.5 mm over a base (a print-in-place lid), and a recess 20 mm
; wide and 0.5 mm deep in an underside on the bed all fail at 10 mm.
(def slot (subtract (move (box 28 12 8) 0 0 4) (move (box 20 20 0.4) 0 0 4)))
(def lid (subtract (move (box 40 10 6) 0 0 3) (move (box 30 20 0.5) 0 0 3.25)))
(def recess (subtract (move (box 28 12 6) 0 0 3) (move (box 20 20 1) 0 0 0)))
(body slot (shape slot) (attr material "petg"))
(body lid (shape lid) (attr material "petg"))
(body recess (shape recess) (attr material "petg"))
; expect-status: slot-bridged fails
(assert slot-bridged (empty (unbridged slot 0 limit step span 0.8)))
; expect-status: lid-bridged fails
(assert lid-bridged (empty (unbridged lid 0 limit step span 0.8)))
; expect-status: recess-bridged fails
(assert recess-bridged (empty (unbridged recess 0 limit step span 0.8)))

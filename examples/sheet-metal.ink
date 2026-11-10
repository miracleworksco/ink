; Sheet metal designed flat first (std's sheet): generic parts whose folded
; sizes are checked against hand calculation. 2 mm sheet, inside radius 2,
; K-factor 0.4. So a 90 degree bend's allowance is
; BA = (pi / 2)(2 + 0.4 x 2) = 4.39823 mm (SOLIDWORKS Help, K-Factor).
; Its zone is 2.19911 either side of its line. A 60 mm flat bent at x = 30
; is flat for 27.80089 either side of the zone. Each leg's outside length
; adds the outside setback r + th = 4: 31.80089. Twice that, less the bend
; deduction 2 x 4 - 4.39823 = 3.60177, is the 60 mm cut.
;
; ink slice examples/sheet-metal.ink --target l-flat --out l-flat.dxf   (the cut file)
; expect-status: short-flange fails
; expect-status: hole-near-bend fails
; expect-status: unrelieved-tab fails
; expect-status: overlapping-flanges fails
; expect-status: overlapping-zones fails
; expect-status: short-web fails
; expect-status: closed-short-flange fails
; expect-status: hole-in-zone fails
; expect-status: shallow-relief fails
; expect-status: ear-folds fails
; expect-status: b-end-relief fails
; expect-status: zones-apart-bends fails
; expect-status: walls-cross fails
; expect-status: lip-tears fails
; expect-status: lip-meets-zone fails
; expect-status: zones-meet fails
; expect-status: zone-off-flange fails
; expect-status: across-acute fails
; expect-status: hole-near-angled-end fails
(import "std")

(def th 2)
(def k 0.4)
(def r 2)

; An L bracket: 60 x 40 flat, bent up 90 degrees along x = 30, a hole in
; each leg.
(def l (sheet (polygon (vec (vec 0 0) (vec 60 0) (vec 60 40) (vec 0 40)))
              (vec (move (circle-2d 2.5) 12 20 0) (move (circle-2d 2.5) 48 20 0)) th k
              (vec (sheet-bend (vec 30 40) (vec 30 0) 90 r))))
(def lb (sheet-folded l))
(body l-bracket (shape lb)
  (attr material "aluminum-5052")
  ; the vertical leg's outside face at its hole, z out of the part (+x)
  (frame upright (sheet-frame l 48 20 -1)))
(def lf (sheet-flat l))
(body l-flat (shape lf) (mark bend (sheet-bend-lines l)))   ; the cut file, its bend line dashed
(assert l-outside-x-in (subset (probe 31.79 20 10) lb))
(assert l-outside-x-out (disjoint lb (probe 31.81 20 10)))
(assert l-inside-x-in (subset (probe 29.81 20 10) lb))
(assert l-inside-x-out (disjoint lb (probe 29.79 20 10)))
(assert l-top-in (subset (probe 30.8 20 31.79) lb))
(assert l-top-out (disjoint lb (probe 30.8 20 31.81)))
; the bend's outside arc, radius 4 about (27.80089, 4), at 45 degrees
(assert l-arc-in (subset (probe 30.6224 20 1.1785) lb))
(assert l-arc-out (disjoint lb (probe 30.6366 20 1.1644)))
; the hole in the upright, 48 - 32.19911 = 15.80089 past its zone: centre z 19.80089
(assert l-hole-clear (disjoint lb (probe 30.8 20 19.8)))
; its frame on the upright's outside face, 5 along its x from the hole
(def upright (sheet-frame l 48 20 -1))
(assert l-frame-behind (subset (place (probe 5 0 -0.05) upright) lb))
(assert l-frame-ahead (disjoint lb (place (probe 5 0 0.05) upright)))
; where the flat meets the bend, mid-thickness: the flats reach into the
; zone, so the union is inside there, with no seam at 0
(assert l-seam-base (subset (probe (- 30 (* 0.5 (sheet-bend-allowance 90 r th k))) 20 1) lb))
(assert l-seam-upright (subset (probe 30.8 20 4) lb))
(assert l-flanges (empty (sheet-short-flanges l 4)))
(assert l-holes (empty (sheet-holes-near-bends l 2.5)))
(assert l-reliefs (empty (sheet-unrelieved l)))
(assert l-collisions (empty (sheet-collisions l)))
; past the zone's 90 degree end, 160 round its axis at radius 3: air
(assert l-past-zone-end (disjoint lb (probe 28.82699 20 6.81908)))
(assert l-radius (value (sheet-min-radius l) >= th))

; The same flat bent 135 degrees: BA = (3 pi / 4)(2.8) = 6.59734, so the
; zone starts at 26.70133 and its axis is (26.70133, 4); the flange leaves
; the arc along (-0.70711, 0.70711), its mid-surface 3 from the axis, and
; runs 26.70133 to (9.94196, 25.00201).
(def l135 (sheet (polygon (vec (vec 0 0) (vec 60 0) (vec 60 40) (vec 0 40))) (vec) th k
                 (vec (sheet-bend (vec 30 40) (vec 30 0) 135 r))))
(def l135b (sheet-folded l135))
(assert l135-inside-in (subset (probe 21.0515 20 12.4924) l135b))
(assert l135-inside-out (disjoint l135b (probe 21.0373 20 12.4782)))
(assert l135-outside-in (subset (probe 22.4515 20 13.8924) l135b))
(assert l135-outside-out (disjoint l135b (probe 22.4657 20 13.9066)))
(assert l135-end-in (subset (probe 9.94903 20 24.99494) l135b))
(assert l135-end-out (disjoint l135b (probe 9.93489 20 25.00908)))
(def leaning (sheet-frame l135 50 20 -1))
(assert l135-frame-behind (subset (place (probe 0 0 -0.05) leaning) l135b))
(assert l135-frame-ahead (disjoint l135b (place (probe 0 0 0.05) leaning)))
(assert l135-flanges (empty (sheet-short-flanges l135 4)))
(assert l135-collisions (empty (sheet-collisions l135)))

; A Z: 80 x 30 flat, up 90 along x = 25, then down 90 along x = 50 (on the
; first bend's flange, so it folds with it). The web is flat from 27.19911
; to 47.80089 (20.60177), upright at x 24.80089..26.80089 from z 4; the
; second axis is 2 past its outside face at the web's top, z 24.60177, and
; the top flange runs from x 28.80089 for 80 - 52.19911 = 27.80089, its
; inside (the flat's bottom) face at z 26.60177, its top at 28.60177.
(def zs (sheet (polygon (vec (vec 0 0) (vec 80 0) (vec 80 30) (vec 0 30))) (vec) th k
               (vec (sheet-bend (vec 25 30) (vec 25 0) 90 r)
                    (sheet-bend (vec 50 30) (vec 50 0) -90 r))))
(def zb (sheet-folded zs))
(body z-bracket (shape zb)
  (frame top (sheet-frame zs 70 15 1)))
(assert z-web-in (subset (probe 26.79 15 20) zb))
(assert z-web-out (disjoint zb (probe 26.81 15 20)))
(assert z-top-in (subset (probe 45 15 28.59) zb))
(assert z-top-out (disjoint zb (probe 45 15 28.61)))
(assert z-under-in (subset (probe 45 15 26.61) zb))
(assert z-under-out (disjoint zb (probe 45 15 26.59)))
(assert z-end-in (subset (probe 56.59 15 27.6) zb))
(assert z-end-out (disjoint zb (probe 56.61 15 27.6)))
(assert z-flanges (empty (sheet-short-flanges zs 4)))
(assert z-collisions (empty (sheet-collisions zs)))

; A tray: a 100 x 60 base with four walls 15 past their zones, bent up
; along its edges. Its corners are notched 6.5 square, each wall that short
; of them: the reliefs (BA / 2 + th + r = 6.19911 the least).
(def d 6.5)
(def wall (+ (* 0.5 (sheet-bend-allowance 90 r th k)) 15))
(def ts (sheet (polygon (vec (vec d (neg wall)) (vec (- 100 d) (neg wall)) (vec (- 100 d) d) (vec (+ 100 wall) d)
                             (vec (+ 100 wall) (- 60 d)) (vec (- 100 d) (- 60 d)) (vec (- 100 d) (+ 60 wall))
                             (vec d (+ 60 wall)) (vec d (- 60 d)) (vec (neg wall) (- 60 d)) (vec (neg wall) d) (vec d d)))
               (vec) th k
               (vec (sheet-bend (vec (- 100 d) 0) (vec d 0) 90 r) (sheet-bend (vec 100 (- 60 d)) (vec 100 d) 90 r)
                    (sheet-bend (vec d 60) (vec (- 100 d) 60) 90 r) (sheet-bend (vec 0 d) (vec 0 (- 60 d)) 90 r))))
(def tb (sheet-folded ts))
(body tray (shape tb)
  (frame floor (sheet-frame ts 50 30 1)))
(assert tray-wall-in (subset (probe 101.79 30 10) tb))
(assert tray-wall-out (disjoint tb (probe 101.81 30 10)))
(assert tray-wall-top-in (subset (probe 100.8 30 18.99) tb))
(assert tray-wall-top-out (disjoint tb (probe 100.8 30 19.01)))
(assert tray-front-in (subset (probe 50 -1.79 10) tb))
(assert tray-front-out (disjoint tb (probe 50 -1.81 10)))
(assert tray-corner-open (disjoint tb (probe 98 -1 10)))
(assert tray-flanges (empty (sheet-short-flanges ts 4)))
(assert tray-reliefs (empty (sheet-unrelieved ts)))
(assert tray-collisions (empty (sheet-collisions ts)))

; A bend at an angle to the strip's edges: 100 x 30, up 90 along
; x = 40 + y / 3 (18.4 degrees off y), its ends squared to the bend 8
; either side of its line and a to b between them. Its flange's flat
; at the foot (45, 15) runs 47.43416 along n = (0.94868, -0.31623) to
; the edge y = 0; folded, the wall's inside face is r - BA / 2 = -0.19911
; along n from the line and its outside 1.80089, and the flat's end is
; th + r + 47.43416 - BA / 2 = 49.23505 up. The hole at (80, 15), 33.20392
; along n, is 35.00480 up, its centre on the outside face at
; (50.20847, 24.93051).
(def skew (sheet (polygon (vec (vec 0 0) (vec 30 0) (vec 33.41053 5.52982) (vec 48.58947 0.47018) (vec 55 0) (vec 100 0)
                               (vec 100 30) (vec 65 30) (vec 56.58947 24.47018) (vec 41.41053 29.52982) (vec 35 30) (vec 0 30)))
                 (vec (move (circle-2d 2.5) 80 15 0)) th k
                 (vec (sheet-bend (vec 49 27) (vec 41 3) 90 r))))
(def sb (sheet-folded skew))
(body skewed (shape sb) (frame wall (sheet-frame skew 80 15 -1)))
(assert skew-outside-in (subset (probe 46.69814 14.43395 20) sb))
(assert skew-outside-out (disjoint sb (probe 46.71712 14.42763 20)))
(assert skew-inside-in (subset (probe 44.81975 15.06008 20) sb))
(assert skew-inside-out (disjoint sb (probe 44.80078 15.06641 20)))
(assert skew-top-in (subset (probe 45.75894 14.74684 49.2) sb))
(assert skew-top-out (disjoint sb (probe 45.75894 14.74684 49.27)))
(assert skew-hole (disjoint sb (probe 49.25979 25.24674 35.00480)))
(assert skew-hole-rim (subset (probe 49.25979 25.24674 38.0) sb))
(assert skew-flanges (empty (sheet-short-flanges skew 4)))
(assert skew-holes (empty (sheet-holes-near-bends skew 2.5)))
(assert skew-reliefs (empty (sheet-unrelieved skew)))
(assert skew-collisions (empty (sheet-collisions skew)))
; The same bend across the plain strip, a and b past its edges so its
; zone spans the material: past the edges is air, not a short flange,
; but where an edge meets the bend at 71.6 degrees its flange tapers
; under 4 th.
(def across (sheet (polygon (vec (vec 0 0) (vec 100 0) (vec 100 30) (vec 0 30))) (vec) th k
                   (vec (sheet-bend (vec 51 33) (vec 39 -3) 90 r))))
(assert across-reliefs (empty (sheet-unrelieved across)))
(assert across-air (disjoint (sheet-short-flanges across 4) (probe 35.92 -0.54 1)))
(assert across-acute (empty (sheet-short-flanges across 4)))
; A 45 degree bend across a plain plate, its line leaving the outline at
; (40, 0). Past there, level with the foot (39.5, -0.5), the zone's
; material is only its base side (BA / 2 = 2.19911 across from the line,
; the plate's edge 0.70711 from it), so the bend is there and its flange
; is not: air 3 past the line is a short flange (z 0.318 is the place
; across at -1.5). A hole of r 0.5 at (35, 3), its projection on the line
; past the outline too, is 2.96 from that material.
(def angled (sheet (polygon (vec (vec 0 0) (vec 100 0) (vec 100 30) (vec 0 30)))
                   (vec (move (circle-2d 0.5) 35 3 0)) th k
                   (vec (sheet-bend (vec 75 35) (vec 35 -5) 90 r))))
(assert angled-end-short (subset (probe 41.62132 -2.62132 0.318) (sheet-short-flanges angled 4)))
(assert hole-near-angled-end (empty (sheet-holes-near-bends angled 2.5)))

; An open hem: a 30 mm strip folded 180 degrees back over its top along
; x = 50, inside radius 1 (inside diameter th, Protolabs' least). BA =
; pi (1 + 0.8) = 5.65487, so the zone starts at 47.17257, the fold's
; outside reaches r + th past it, 50.17257, and the return, the 12 of
; flat past the zone, lies over the base from 47.17257 back to 35.17257,
; its inside face 2r above the base's top, z 4 to 6.
(def hem (sheet (polygon (vec (vec 0 0) (vec 64.82743 0) (vec 64.82743 30) (vec 0 30))) (vec) th k
                (vec (sheet-bend (vec 50 30) (vec 50 0) 180 1))))
(def hb (sheet-folded hem))
(body hemmed (shape hb))
(assert hem-gap (disjoint hb (probe 40 15 3)))
(assert hem-under-in (subset (probe 40 15 4.01) hb))
(assert hem-under-out (disjoint hb (probe 40 15 3.99)))
(assert hem-top-in (subset (probe 40 15 5.99) hb))
(assert hem-top-out (disjoint hb (probe 40 15 6.01)))
(assert hem-end-in (subset (probe 35.18 15 5) hb))
(assert hem-end-out (disjoint hb (probe 35.16 15 5)))
(assert hem-fold-in (subset (probe 50.16 15 3) hb))
(assert hem-fold-out (disjoint hb (probe 50.18 15 3)))
(assert hem-flanges (empty (sheet-short-flanges hem 4)))
(assert hem-collisions (empty (sheet-collisions hem)))

; Controls, each failing its rule: a flange 2.8 past its zone (4 th from
; the mold line needs 8 - 4 = 4); a hole 5 from the bend line (2.5 th + r
; = 7 needs the edge 7 away); a 20 mm tab bent from the middle of a 60 mm
; edge with no reliefs beside it
(def short (sheet (polygon (vec (vec 0 0) (vec 35 0) (vec 35 40) (vec 0 40))) (vec) th k
                  (vec (sheet-bend (vec 30 40) (vec 30 0) 90 r))))
(assert short-flange (empty (sheet-short-flanges short 4)))
(def holed (sheet (polygon (vec (vec 0 0) (vec 60 0) (vec 60 40) (vec 0 40))) (vec (move (circle-2d 2.5) 22.5 20 0)) th k
                  (vec (sheet-bend (vec 30 40) (vec 30 0) 90 r))))
(assert hole-near-bend (empty (sheet-holes-near-bends holed 2.5)))
(def tab (sheet (polygon (vec (vec 0 0) (vec 30 0) (vec 30 20) (vec 50 20) (vec 50 0) (vec 60 0)
                              (vec 60 -30) (vec 0 -30)))
                (vec) th k (vec (sheet-bend (vec 30 0) (vec 50 0) 90 r))))
(assert unrelieved-tab (empty (sheet-unrelieved tab)))
; and two bends on one flat whose flanges both take a flare of the first's
(def overlap (sheet (polygon (vec (vec 0 0) (vec 90 0) (vec 90 90) (vec 75 90) (vec 75 56) (vec 56 56)
                                  (vec 56 75) (vec 0 75)))
                    (vec) th k (vec (sheet-bend (vec 60 56) (vec 60 0) 90 r) (sheet-bend (vec 0 60) (vec 56 60) 90 r))))
(assert overlapping-flanges (empty (sheet-unrelieved overlap)))
; two bends 3 apart, up then down: their zones (4.39823 wide) share 1.39823
; of the flat, which the second, on the first's flange, would lose
; (crossing bends are not built)
(def close (sheet (polygon (vec (vec 0 0) (vec 80 0) (vec 80 30) (vec 0 30))) (vec) th k
                  (vec (sheet-bend (vec 30 30) (vec 30 0) 90 r) (sheet-bend (vec 33 30) (vec 33 0) -90 r))))
(assert overlapping-zones (empty (sheet-unrelieved close)))
; and a Z whose web is 2 flat between its zones: no zones overlap, but the
; web is under 4 th from either bend's mold line
(def web (sheet (polygon (vec (vec 0 0) (vec 80 0) (vec 80 30) (vec 0 30))) (vec) th k
                (vec (sheet-bend (vec 30 30) (vec 30 0) 90 r) (sheet-bend (vec 36.39823 30) (vec 36.39823 0) -90 r))))
(assert short-web (empty (sheet-short-flanges web 4)))
; each rule at other angles: a flange 0.5 past a 135 degree zone; a hole
; inside a 180 degree, r 10 zone (16.96 either side of its line), 16 from
; the line; beside that bend's end, a relief cut only th + r = 12 behind
; its line, leaving flat material in its zone
(def b135 (sheet (polygon (vec (vec 0 0) (vec (+ 30 (* 0.5 (sheet-bend-allowance 135 r th k)) 0.5) 0)
                               (vec (+ 30 (* 0.5 (sheet-bend-allowance 135 r th k)) 0.5) 40) (vec 0 40))) (vec) th k
                 (vec (sheet-bend (vec 30 40) (vec 30 0) 135 r))))
(assert closed-short-flange (empty (sheet-short-flanges b135 4)))
(def c180 (sheet (polygon (vec (vec 0 0) (vec 100 0) (vec 100 40) (vec 0 40))) (vec (move (circle-2d 0.5) 34 20 0)) th k
                 (vec (sheet-bend (vec 50 40) (vec 50 0) 180 10))))
(assert hole-in-zone (empty (sheet-holes-near-bends c180 2.5)))
(def r180 (sheet (polygon (vec (vec 0 0) (vec 100 0) (vec 100 30) (vec 38 30) (vec 38 40) (vec 0 40))) (vec) th k
                 (vec (sheet-bend (vec 50 30) (vec 50 0) 180 10))))
(assert shallow-relief (empty (sheet-unrelieved r180)))
; a tab bent up from an edge, relieved, beside an ear past its line meant
; to stay flat: the ear would fold with the tab
(def ear (sheet (polygon (vec (vec 0 0) (vec 80 0) (vec 80 10) (vec 53.5 10) (vec 53.5 12) (vec 60 12)
                              (vec 60 30) (vec 80 30) (vec 80 40) (vec 0 40)))
                (vec) th k (vec (sheet-bend (vec 60 10) (vec 60 0) 90 r))))
(assert ear-folds (empty (sheet-unrelieved ear)))
; a tab bent from the middle of an edge, relieved th + r past its zone at
; its a end (x 17.5..20 down to -6.5) but only to -5 at its b end
(def bend-ends (sheet (polygon (vec (vec 0 -30) (vec 60 -30) (vec 60 0) (vec 42.5 0) (vec 42.5 -5) (vec 40 -5)
                                    (vec 40 20) (vec 20 20) (vec 20 -6.5) (vec 17.5 -6.5) (vec 17.5 0) (vec 0 0)))
                      (vec) th k (vec (sheet-bend (vec 20 0) (vec 40 0) 90 r))))
(assert b-end-relief (empty (sheet-unrelieved bend-ends)))
; two unrelated bends 3 apart folding away from each other: their zones
; share 1.4 of the flat
(def apart (sheet (polygon (vec (vec 0 -20) (vec 50 -20) (vec 50 30) (vec 0 30))) (vec) th k
                  (vec (sheet-bend (vec 50 0) (vec 0 0) 90 r) (sheet-bend (vec 0 3) (vec 50 3) 90 r))))
(assert zones-apart-bends (empty (sheet-unrelieved apart)))
; a U whose walls lean in at 135 degrees, 37 past their zones: folded, they
; pass through each other above the base
(def leaning-u (sheet (polygon (vec (vec -30 0) (vec 70 0) (vec 70 30) (vec -30 30))) (vec) th k
                      (vec (sheet-bend (vec 10 0) (vec 10 30) 135 r) (sheet-bend (vec 30 30) (vec 30 0) 135 r))))
(assert walls-cross (empty (sheet-collisions leaning-u)))
; an L (r 6) with a lip on its upright, relieved beside it, whose flat past
; the lip's line runs back to the L's zone: it would fold with the lip and
; tear from the L's arc
(def lip (sheet (subtract (polygon (vec (vec 0 0) (vec 70 0) (vec 70 70) (vec 0 70)))
                          (polygon (vec (vec 35.29 40) (vec 39.9 40) (vec 39.9 26.5) (vec 42.1 26.5) (vec 42.1 45.45) (vec 35.29 45.45))))
                (vec) th k (vec (sheet-bend (vec 30 70) (vec 30 0) 90 6) (sheet-bend (vec 42 40) (vec 70 40) 90 6))))
(assert lip-tears (empty (sheet-unrelieved lip)))
; walls up 90 at x = 10 and 50 with lips 30 past their zones: the left
; lip's tip reaches 0.5 into the right wall's lip's zone (a flange meeting
; a zone it is not joined to)
(def lip-zone (sheet (polygon (vec (vec -58.6009 0) (vec 100 0) (vec 100 30) (vec -58.6009 30))) (vec) th k
                     (vec (sheet-bend (vec 10 0) (vec 10 30) 90 r) (sheet-bend (vec 50 30) (vec 50 0) 90 r)
                          (sheet-bend (vec -18.3 0) (vec -18.3 30) 90 r) (sheet-bend (vec 80 30) (vec 80 0) -90 r))))
(assert lip-meets-zone (empty (sheet-collisions lip-zone)))
; walls leaning in at 135 whose lips turn back 90: the lips' zones meet
(def zone-zone (sheet (polygon (vec (vec -34.6969 0) (vec 94.6969 0) (vec 94.6969 30) (vec -34.6969 30))) (vec) th k
                      (vec (sheet-bend (vec 10 0) (vec 10 30) 135 r) (sheet-bend (vec 50 30) (vec 50 0) 135 r)
                           (sheet-bend (vec -20.4978 0) (vec -20.4978 30) -90 r) (sheet-bend (vec 80.4978 30) (vec 80.4978 0) -90 r))))
(assert zones-meet (empty (sheet-collisions zone-zone)))
; a wall up from x = 20 and a bend on it (r 6) whose line runs past the
; wall's end across an arm of the base: its zone there is the base's, and
; would be folded with the wall too
(def arm (sheet (polygon (vec (vec 0 0) (vec 64 0) (vec 64 30) (vec 13 30) (vec 13 60) (vec 8 60) (vec 0 52))) (vec) th k
                (vec (sheet-bend (vec 20 30) (vec 20 0) 90 r) (sheet-bend (vec 6 58) (vec 64 0) 90 6))))
(assert zone-off-flange (empty (sheet-unrelieved arm)))

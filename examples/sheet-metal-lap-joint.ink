; REV Robotics REV-21-2359, Lap Joint Bracket - 1in, designed flat first
; (std's sheet): a U of 3 mm 5052 aluminium, two 90 degree bends from one
; web. The probes prove its folded sizes match REV's drawing.
;
; From the drawing REV-21-2359-DR (12/16/2024, revrobotics.com/1in-Brackets):
; 31.4 outside, 25.4 inside; flanges 57.2 (2.25 in) from the web's outside
; face and 73.0 (2.875 in) long; holes on a 12.7 grid 9.5 below the top
; (the drawing writes 9.5 and 9.7 both as 0.38 in).
; From REV's STEP (REV-21-2359.STEP, its flange outline and bends), what
; the drawing does not dimension: inside radius 0.6096 (0.024 in, 0.2 th,
; tighter than Geomiq's general r >= th); holes 4.9784 (0.196 in, #10
; clearance); the slanted edge's ends (12.7 up at the web's end, 57.15 up
; at 50.8 along), the top's inner corner at 12.7 and the notch at 6.35.
; REV publishes no K-factor: 0.33 is Geomiq's for air bending aluminium at
; a radius under th. It sets only the cut (136.29 wide here); the folded
; part is the drawing's at any K.
;
; ink slice examples/sheet-metal-lap-joint.ink --target flat --out REV-21-2359.dxf
; writes the outline and holes to cut, and the two bend lines, dashed on
; layer bend.
(import "std")

(def th 3)
(def r 0.6096)
(def k 0.33)
(def ba (sheet-bend-allowance 90 r th k))
; The flat: the web between its zones, 31.4 - 2 (r + th) wide, centred on
; x = 0; a flange's point h from the web's outside face (the drawing's
; dimensions) is h - (r + th) past its zone.
(def x0 (+ (* 0.5 (- 31.4 (* 2 (+ r th)))) ba))
(fn fx (h) (+ x0 (- h r th)))
(def xb (- x0 (* 0.5 ba)))
(def hole (circle-2d (/ 4.9784 2)))
(def S
  (sheet (polygon (vec (vec (neg (fx 12.7)) 0) (vec (fx 12.7) 0) (vec (fx 57.15) 50.8) (vec (fx 57.15) 66.675)
                       (vec (fx 50.8) 73.025) (vec (fx 12.7) 73.025) (vec (fx 6.35) 66.675) (vec (fx 6.35) 50.8)
                       (vec (neg (fx 6.35)) 50.8) (vec (neg (fx 6.35)) 66.675) (vec (neg (fx 12.7)) 73.025)
                       (vec (neg (fx 50.8)) 73.025) (vec (neg (fx 57.15)) 66.675) (vec (neg (fx 57.15)) 50.8)))
         (concat (map (lambda (y) (move hole 0 y 0)) (vec 12.7 25.4 38.1))
                 (flatten (map (lambda (h) (vec (move hole (fx h) 63.5 0) (move hole (neg (fx h)) 63.5 0)))
                               (vec 12.7 25.4 38.1 50.8))))
         th k
         (vec (sheet-bend (vec (neg xb) 0) (vec (neg xb) 50.8) 90 r)
              (sheet-bend (vec xb 50.8) (vec xb 0) 90 r))))
(def part-shape (sheet-folded S))
(def flat-shape (sheet-flat S))

; The web's outside face at z = 0, the flanges' outside faces at x = -15.7
; and 15.7; frames on outside faces at the web's middle hole and each
; flange's first hole, z out of the part.
(body part (shape part-shape)
  (attr material "aluminum-5052")
  (frame back (sheet-frame S 0 25.4 -1))
  (frame left (sheet-frame S (neg (fx 12.7)) 63.5 -1))
  (frame right (sheet-frame S (fx 12.7) 63.5 -1)))
(body flat (shape flat-shape) (attr material "aluminum-5052") (mark bend (sheet-bend-lines S)))

(assert outside-in (subset (probe 15.69 25 30) part-shape))
(assert outside-out (disjoint part-shape (probe 15.71 25 30)))
(assert inside-in (subset (probe -12.71 25 30) part-shape))
(assert inside-out (disjoint part-shape (probe -12.69 25 30)))
(assert height-in (subset (probe 14.2 60 57.14) part-shape))
(assert height-out (disjoint part-shape (probe 14.2 60 57.16) ))
(assert length-in (subset (probe -14.2 73.02 40) part-shape))
(assert length-out (disjoint part-shape (probe -14.2 73.03 40)))
(assert web-hole (disjoint part-shape (probe 0 25.4 1.5)))
(assert flange-hole (disjoint part-shape (probe 14.2 63.5 25.4)))
(assert cut-width-in (subset (probe 68.13 60 1.5) flat-shape))
(assert cut-width-out (disjoint flat-shape (probe 68.15 60 1.5)))
(assert flanges (empty (sheet-short-flanges S 4)))
(assert holes (empty (sheet-holes-near-bends S 2.5)))
(assert reliefs (empty (sheet-unrelieved S)))
(assert collisions (empty (sheet-collisions S)))

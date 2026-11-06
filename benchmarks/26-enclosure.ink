; Benchmark 26: injection-moulded electronics enclosure (base + lid).
;
; Task: an ABS box, 80 x 60 x 30 at the parting line (z = 30), 2 mm uniform
; wall, R 3 on every outer edge, 1 deg draft on the side walls (the base
; narrows toward z = 0; the cavity stays 2 mm normal to the wall). Four
; drafted screw bosses (r 3.5 at the top, M3 pilot r 1.25, 27 deep) tied into
; the corners, three 1.2 x 8 ribs on the floor, r 0.5 foot fillets (std's
; fillet-union). A 2 mm lid with a 3 mm locating lip 0.2 clear of the drafted wall,
; stopped square short of the bosses, and r 1.7 screw clearance holes.
;
; Hand calc (tools/enclosure_hand.py, pure Python, independent of the
; kernel): horizontal sections of the drafted, rounded cup are rounded
; rectangles, A = 4ab + 4 rho (a + b) + pi rho^2 (rho = R / cos(draft) on the
; walls, sqrt(R^2 - (R - z)^2) in the bottom fillet), integrated in z; the
; bosses are disks cut by the cavity section per slice.
;   base without the foot fillets: cup 23 794.8 + ribs 1756.1 + bosses
;   3537.7 - pilot holes 530.1 = 28 558.5 mm^3
;   lid: plate 9512.2 + lip 924.2 - stopped ends 156.7 = 10 279.7 mm^3
;   wall 2.0 (normal), ribs 1.2, pilot r 1.25, draft 1 deg (the side plane's
;   normal is 89 deg from the floor's), lip gap 0.2 at its lowest edge.
;
; Pass criteria:
;   - volumes within 0.3% of the hand values (base at foot-r = 0);
;   - wall thickness (the walls claim, by rays): its minimum 1.2 (the ribs), no false thin spots;
;   - one corner's screw holes, by probes (proven): the base's pilot holds
;     a pin of r 1.24 on the boss axis and is metal from r 1.26, the lid's
;     hole holds r 1.69 and is metal from r 1.71 (coaxial within 0.01);
;   - lip clear by 0.18 (proven), and not by 0.24 (disproven: the control);
;     lid and base do not interfere (they touch at z = 30);
;   - both meshes watertight, one shell each.
; expect-status: lip-not-looser fails

(import "std")

; ---------------------------------------------------------------- params
(param L 80)            ; outer length (x) at the parting line
(param W 60)            ; outer width (y) at the parting line
(param H 30)            ; base height; parting line at z = H
(param wall 2)          ; nominal wall
(param R 3)             ; outer edge radius (vertical and bottom edges)
(param draft 1)         ; side-wall draft, degrees (base narrows toward z = 0)
(param boss-r 3.5)      ; boss outer radius at the top
(param hole-r 1.25)     ; M3 self-tapping pilot (2.5 mm)
(param boss-x 35.5) (param boss-y 25.5)   ; boss covers the cavity corner (see Status)
(param rib-t 1.2)       ; ribs 0.6 wall
(param rib-h 8)         ; rib height above the floor
(param tl 2)            ; lid plate thickness
(param lip-h 3)         ; lid locating lip depth
(param lip-t 1.2)       ; lip wall
(param gap 0.2)         ; lip-to-wall clearance at the lip's bottom edge
(param screw-r 1.7)     ; lid clearance hole
(param foot-r 0.5)      ; boss and rib foot fillets

(def ta (tan (deg->rad draft)))
(def ca (cos (deg->rad draft)))

; ---------------------------------------------------------------- base
; A drafted core box (half-widths grow with z, planes tilted by `draft`,
; scaled by cos so each side is an exact plane distance), offset by R: every
; outer edge is rounded by R, the walls stay drafted. Its top runs past the
; parting line; the cut at z = H makes the rim.
(fn drafted-core (hx hy z0 z1)
  (emax (* ca (- (abs x) hx (* ta (- z H))))
        (* ca (- (abs y) hy (* ta (- z H))))
        (- (abs (- z (* 0.5 (+ z0 z1)))) (* 0.5 (- z1 z0)))))
(def outer (- (drafted-core (- (* 0.5 L) R) (- (* 0.5 W) R) R (+ H 10)) R))
(def cavity (+ outer wall))                              ; uniform wall: inset by t
(def cup (max outer (neg cavity) (- z H)))

; Screw bosses: drafted cones (narrower at the top, the core pulls up),
; tied into the corner walls, flush with the parting line.
(fn at-bosses (f) (min (move f boss-x boss-y 0) (move f (neg boss-x) boss-y 0)
                       (move f boss-x (neg boss-y) 0) (move f (neg boss-x) (neg boss-y) 0)))
(def boss (at-bosses (max (* ca (- (len2 x y) boss-r (* ta (- H z)))) (- 1 z) (- z H))))
(def hole (at-bosses (max (- (len2 x y) hole-r) (- (+ wall 1) z))))

; Ribs: two across y at x = +-12, one along x at y = 0, rib-h above the floor.
(def rib-slab (max (- z (+ wall rib-h)) (- 1 z)))
(def rib (min (max (- (abs (- (abs x) 12)) (* 0.5 rib-t)) (- (abs y) (- (* 0.5 W) 1)) rib-slab)
              (max (- (abs y) (* 0.5 rib-t)) (- (abs x) (- (* 0.5 L) 1)) rib-slab)))
;  Foot fillets r 0.5: a fillet-union of the cup with the bosses and ribs,
; which also swells material along seams where the two sides are coplanar
; (boss tops flush with the rim, bosses clipped to the outer skin): the lid
; then interferes (see Status). So the bosses and ribs are clipped to the
; outer skin and the parting plane after the fillet, not before.
(def inside-bits (min boss rib))
(def filleted (subtract (fillet-union cup inside-bits foot-r 90) hole))
(def base-field (max filleted outer (- z H)))

(body base (shape base-field) (attr material "abs") (attr color "#3A4250"))

; ---------------------------------------------------------------- lid
; Plate: the base's outline at the parting line, extruded up by tl.
(def lid-plate (max (warp outer x y H) (- H z) (- z (+ H tl))))
; Lip: a ring inside the opening, gap clear of the wall at its lowest point
; (the drafted wall only opens upward from there). It stops square short of
; the corner bosses: a round notch (at-bosses cylinder) first left wedge-thin
; lip ends (0.25 mm, found by ray thickness).
(def lip-out (+ (warp cavity x y (- H lip-h)) gap))
(def lip (max lip-out (neg (+ lip-out lip-t)) (- (- H lip-h) z) (- z H)))
(def notch (max (- (- boss-x (+ boss-r gap 0.15)) (abs x)) (- (- boss-y (+ boss-r gap 0.15)) (abs y))))
(def screw-hole (at-bosses (- (len2 x y) screw-r)))
(def notched-lip (subtract lip notch))
(def lid-field (subtract (min lid-plate notched-lip) screw-hole))

(body lid (shape lid-field) (attr material "abs") (attr color "#5A6A80"))

(assembly root
  (instance base base)
  (instance lid lid))

; ---------------------------------------------------------------- claims
(assert lip-clear (clearance (shape-of root/lid notched-lip) (shape-of root/base) 0.18))
(assert lip-not-looser (clearance (shape-of root/lid notched-lip) (shape-of root/base) 0.24) (accuracy 0.02))
; The lid seats on the base: no search separates the faces at the model's
; floor within its budget, so it is decided at cells of 0.05 (to 0.1 mm).
(assert lid-seated (disjoint (shape-of root/lid) (shape-of root/base)) (res 0.05))
(assert walls (min-thickness (shape-of base) > 1.1))
; One corner's screw holes on the boss axis (35.5, 25.5): open to a pin
; 0.01 under their radius, metal in a ring from 0.01 over it.
(fn ring (r0 r1 z0 z1) (move (subtract (cylinder-z r1 (- z1 z0)) (cylinder-z r0 (+ (- z1 z0) 1))) 35.5 25.5 (* 0.5 (+ z0 z1))))
(assert pilot-open (disjoint (shape-of root/base) (move (cylinder-z 1.24 24) 35.5 25.5 17)))
(assert pilot-metal (subset (ring 1.26 1.6 6 28) (shape-of root/base)))
(assert screw-hole-open (disjoint (shape-of root/lid) (move (cylinder-z 1.69 1.6) 35.5 25.5 31)))
(assert screw-hole-metal (subset (ring 1.71 2 30.2 31.8) (shape-of root/lid)))

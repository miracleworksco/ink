; Benchmark 31: a 6061 sheet-metal channel bracket with two bends.
;
; Task: a 2 mm sheet, 40 mm wide: a 40 mm tangent-to-tangent web, two 20 mm
; straight flanges, two 90 degree bends of inside radius 3 mm, two d5 holes
; in the right flange (y = +-10, z = 16), and one 3 mm wide bend relief from
; the +y free edge, across the right bend plus 1 mm of each flat. Formed
; section: web x = -20..20, z = 0..2; bend centres (+-20, 5), outer radius 5;
; flanges at x = +-24, z = 5..25. Density 0.0027 g/mm^3. Check the hole
; ligaments against shop rules, and the stock against the formed part.
;
; Hand calculation (development with K = 0.33, a process assumption):
;   - BA = (pi/2)(Ri + K t) = 5.749115 mm per bend; flat length
;     40 + 2 x 20 + 2 BA = 91.498229 mm.
;   - flat area 40 x length - 2 pi 2.5^2 - 3 (BA + 2) = 3597.4119 mm^2;
;     volume 7194.8238 mm^3, mass 19.42602 g.
;   - The formed solid has constant thickness, so it conserves volume at
;     the mid-radius (K = 0.5): BA_mid = (pi/2)(Ri + t/2) = 6.283185, giving
;     7277.0707 mm^3 and 19.64809 g, +1.143% over the K = 0.33 blank. That
;     is this model's constant-thickness approximation (real bends thin and
;     shift the neutral line), not an integration error.
;   - ligaments: top 6.5, sides 7.5, hole to bend tangent 8.5 mm. Shop rules
;     for this benchmark (not universal standards): hole edge to free edge
;     >= 2 r (centre >= 1.5 d), hole edge to bend >= Ri + 2 t = 7, hole edge
;     to relief >= 5.
;
; Pass criteria (the claims; ink volume, ink mass and ink assert --param, b31 in
; crates/ink-cli/tests/cli/benchmarks.rs):
;   - check clean; assert --res 0.1: the three ligament rules, the open core
;     through each hole (a 4.9 mm cylinder) proven; the deliberately
;     excessive 7 mm top-edge rule disproven. No tolerance-only results.
;   - volume and mass --res 0.1: estimate within 0.1% of the formed
;     analytic value; within 1.5% of
;     the K = 0.33 blank (the stated approximation; 0.5% would fail).
;   - Controls: --param hole-z=19 disproves hole-edge-distance, hole-z=11
;     disproves hole-bend-distance.

(import "std")
(param sheet 2)
(param ri 3)
(param K 0.33)
(param web 40)
(param flange 20)
(param width 40)
(param hole-r 2.5)
(param hole-z 16)
(param hole-y 10)
(param relief-width 3)
(def ro (+ ri sheet))
(def tangent-z ro)
(def top-z (+ ro flange))
(def x-tan (* 0.5 web))
(def x-wall (+ x-tan ri (* 0.5 sheet)))
(def base (move (box web width sheet) 0 0 (* 0.5 sheet)))
(fn corner (xc side)
  (max (- ri (len2 (- x xc) (- z ro)))
       (- (len2 (- x xc) (- z ro)) ro)
       (* (neg side) (- x xc)) (- z ro) (- (abs y) (* 0.5 width))))
(def bend-right (corner x-tan 1))
(def bend-left (corner (neg x-tan) -1))
(def wall-right (move (box sheet width flange) x-wall 0 (+ ro (* 0.5 flange))))
(def wall-left (move (box sheet width flange) (neg x-wall) 0 (+ ro (* 0.5 flange))))
; Holes only through right flange: finite cylinders, shared by cut and claims.
(def holes (min (move (cylinder-x hole-r (+ sheet 2)) x-wall hole-y hole-z)
                (move (cylinder-x hole-r (+ sheet 2)) x-wall (neg hole-y) hole-z)))
; Notch removes one full bend-width strip, plus 1 mm into each adjacent flat.
(def relief (max (- (- x-tan 1) x) (- z (+ ro 1)) (- (- (* 0.5 width) relief-width) y)))
(body bracket (shape (subtract (min base bend-left bend-right wall-left wall-right) (min holes relief)))
  (attr material "aluminum-6061") (attr density 2.7g/cm3))
; Explicit manual development: K is a process assumption, not derived by ink.
(def BA (* (/ pi 2) (+ ri (* K sheet))))
(def flat-length (+ web (* 2 flange) (* 2 BA)))
(def flat-area (- (* width flat-length) (* 2 pi (square hole-r)) (* relief-width (+ BA 2))))
(def flat-volume (* sheet flat-area))
(def flat-mass (* 0.0027 flat-volume))
; A uniform-thickness circular solid conserves area at K=0.5, not K=0.33.
(def geometric-BA (* (/ pi 2) (+ ri (* 0.5 sheet))))
(def geometric-volume (* sheet (- (* width (+ web (* 2 flange) (* 2 geometric-BA)))
  (* 2 pi (square hole-r)) (* relief-width (+ geometric-BA 2)))))
; Geometric keep-out sets: the outside of the flange's free top/side edges.
; Required hole-edge ligament >= 1.5*d - r = 2*r (centre distance >= 1.5*d).
; Negative control requests 7.0 mm vs the actual 6.5 mm top ligament.
(def free-edge-exterior (min (- top-z z) (- (* 0.5 width) (abs y))))
(assert hole-edge-distance (clearance holes free-edge-exterior (* 2 hole-r)))
(assert hole-bend-distance (clearance holes (min bend-left bend-right) (+ ri (* 2 sheet))))
(assert hole-relief-distance (clearance holes relief (* 2 hole-r)))
; expect-status: control-top-too-close-disproven fails
(assert control-top-too-close-disproven (clearance holes free-edge-exterior 7))
(assert hole-open-core (disjoint (- holes -0.05) (shape-of bracket)))

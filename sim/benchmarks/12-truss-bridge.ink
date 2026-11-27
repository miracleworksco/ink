; Benchmark 12: Pratt truss bridge under a midspan load.
;
; Task: a planar steel Pratt truss, 6 panels of 60 x 60 mm (span 360), round
; members r = 3 swept along their centre lines (tube = exact distance to a
; segment). Verticals are a linear pattern; the diagonals are one side's two
; members mirrored about midspan. Pin support at the left bottom node, roller
; (z only) at the right, each on a bearing: a steel pad 16 x 10 x 4 mm under
; the node (z -6 to -2, 1 mm into the chord, filleted to it as the joints
; are), held on its underside (z < -5), so the hold's singular edge lies in
; the pad, away from the members' joints. 600 N down at the centre bottom
; node b3.
;
; Hand calc (benchmarks/tools/truss.py, pin-jointed, direct stiffness):
;   top chord t2-t3: -900 N -> -31.83 MPa (compression, the worst member)
;   bottom chord b2-b3: +600 N -> +21.22 MPa
;   diagonals: 424 N -> 15.0 MPa; b3-t3 and b1-t1 carry nothing
;   midspan deflection: 0.0907 mm
;   top chord buckling: Euler, pinned, L = 60: pi^2 E I / L^2 = 34.9 kN (x39)
; The FEA has rigid welded joints with material blobs at the nodes, so it
; should be a little stiffer, with some secondary bending. Each joint's
; concave corners are filleted as a weld's would be, at a radius of about
; 1.5 mm (half a member's radius): members meet by a smooth union of width
; 1.75 mm, which fills a square corner as a 1.5 mm fillet does (its depth
; along the bisector, 0.25 k, is a radius-R fillet's for R = 0.85 k), and
; leaves each member's surface as it is wherever no other member is near.
;
; Pass criteria (ink solve --points, ink require):
;   - axial stress at mid-member (on the centre line): top chord t2-t3
;     within 10% of -31.83 MPa, bottom chord b2-b3 within 10% of +21.22;
;   - midspan deflection within 20% of 0.0907 mm;
;   - min safety factor (yield 250) > 3 (false as drawn: below).

(import "std")

(param panel 60)
(param h 60)
(param r 3)
(param load 600)
(param weld 1.75)   ; the joints' blend width: a fillet of about 1.5 mm
(def span (* 6 panel))

(fn join (a b) (smooth-union a b weld))
(fn node-b (i) (vec (* i panel) 0 0))
(fn node-t (i) (vec (* i panel) 0 h))
(fn member (a b) (tube (vec a b) r))

; mirror copy about midspan (std symmetric-x keeps only the + side, which here
; is empty: that first version silently dropped every diagonal)
(fn symmetric-about-mid (f) (join f (warp f (- span x) y z)))
(def bottom-chord (tube (vec (node-b 0) (node-b 6)) r))
(def top-chord (tube (vec (node-t 1) (node-t 5)) r))
(def end-posts (symmetric-about-mid (member (node-b 0) (node-t 1))))
(def verticals (pattern-linear (member (node-b 1) (node-t 1)) 5 panel 0 0))
(def diagonals (symmetric-about-mid (join (member (node-t 1) (node-b 2)) (member (node-t 2) (node-b 3)))))
; Bearing pads under the end nodes, and the region of each held: the pad's
; underside.
(fn pad (i) (move (box 16 10 4) (* i panel) 0 -4))
(fn pad-base (i) (max (+ z 5) (- (abs (- x (* i panel))) 8.5) (- (abs y) 5.5)))
(def truss (join (join (join (join (join bottom-chord top-chord) end-posts) verticals) diagonals) (min (pad 0) (pad 6))))
(body bridge (shape truss) (attr material "steel") (attr color "#8C6D46"))

(study midspan (shape (shape-of bridge))
  (elastic (material steel) (yield 250MPa))
  (fix pin (pad-base 0))
  (fix roller (pad-base 6) (dofs z) (support))
  (load deck (force 0 0 (* 1N (neg load))) (region (- (len3 (- x (* 3 panel)) y z) 5)))
  (res 0.5))

; The solve's linear buckling analysis finds the truss buckles at 8.5
; times the load: the top chord
; swaying out of the truss's plane (the roller leaves its end free across
; it), lower than the hand calc's in-plane 39. The loads do not buckle it,
; so no claim leaves buckling out; (min-buckling-factor midspan > 4) would
; claim a margin.
; Filleted, the joints are resolved: no gap between members is glued
; shut, and on cells of 0.5 mm (whose coarsest check, 2 mm, still holds a
; member three cells across) the deflection is decided. The strength is
; not, though it converges: its peak is in the crotch where an end post
; meets the bottom chord (its 45 degree corner blended to a fillet of
; about 0.7 mm), away from the bearings' holds, and it reads 5.86, 4.81,
; 4.11 and 3.69 on cells of 2, 1, 0.5 and 0.25 mm, settling at about 3.1:
; within its grid error of the claim's 3 on any cells this file can take.
; expect-status: strong unknown
(require strong (min-safety-factor midspan > 3))
(require stiff (max-displacement midspan < 0.1088))
(require not-too-stiff (max-displacement midspan > 0.0726))

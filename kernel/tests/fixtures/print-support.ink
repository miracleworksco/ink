; A lower floor cannot support a separate thin roof. The field of the
; second cavity is exactly zero throughout its air.
; expect-status: printable fails
; expect-status: smaller-step fails
; expect-status: flat-air-roof fails
; expect-status: flat-air-pull fails
(import "std")

(def floor (extrude-z (rect-2d 12 12) 0 2))
(def roof (extrude-z (rect-2d 12 12) 2.5 2.9))
(def part (min floor roof))
(require printable (empty (overhangs part 0 45 1)) (accuracy 0.01))
(require smaller-step (empty (overhangs part 0 45 0.25)) (accuracy 0.01))

(def stock (box 10 10 10))
(def pocket (box 8 8 6))
(def ordinary (subtract stock pocket))
(def flat-air (max stock (min ordinary 0)))
(def views (vec (vec 0 0 1) (vec 0 0 -1)))
(require flat-air-roof (empty (overhangs flat-air -5 45 1)) (accuracy 0.01))
(require flat-air-pull (empty (hidden flat-air stock (clear-air flat-air stock 0.5) views)) (accuracy 0.01))

; Both the roof and its gap are smaller than the requested accuracy.
; Even a clipped check has to prove support before it can pass.
(def foil (min floor (extrude-z (rect-2d 12 12) 2.003 2.005)))
; expect-status: clipped-foil fails
(require clipped-foil (empty (max (overhangs foil 0 45 1) (slab-x -1 1))) (accuracy 0.5))
(def foil-mask (overhangs foil 0 45 1))
; expect-status: foil-disjoint fails
; expect-status: foil-subset fails
(require foil-disjoint (disjoint (slab-x -1 1) foil-mask) (accuracy 0.5))
(require foil-subset (subset foil-mask 1) (accuracy 0.5))
(require on-bed (empty (overhangs floor 0 45 1)))

; The zero set and both signs survive this nonlinear re-encoding.
(def encoded (* part (+ 1 (square part))))
; expect-status: nonlinear-roof fails
(require nonlinear-roof (empty (overhangs encoded 0 45 1)) (accuracy 0.01))

; A slender column leans 35 degrees, halfway between two tested directions.
; Its straight interior support misses the rays. A failure must also rule
; out that support across the whole section.
(def stick (extrude-z (rect-2d 0.3 0.3) 0 8))
(def leaning (warp stick (- x (* 0.7 z (cos (deg->rad 22.5)))) (- y (* 0.7 z (sin (deg->rad 22.5)))) z))
; expect-status: between-directions unknown
(require between-directions (empty (overhangs leaning 0 45 1)) (accuracy 0.1))

; Both pieces have a zero-valued face at their shared height. A single
; plane of zeros separates material on both sides of the joint.
(def flush (min floor (extrude-z (rect-2d 12 12) 2 4)))
(require touching-roof (empty (overhangs flush 0 45 1)) (accuracy 0.01))

(def wall (move (extrude-z (rect-2d 1 12) 1.5 2.8) -5.5 0 0))
; expect-status: edge-attached fails
(require edge-attached (empty (overhangs (min floor wall roof) 0 45 1)) (accuracy 0.01))

(def neighbours (min (move (box 2 20 8) -1.3 0 4) (move (box 2 20 8) 1.3 0 4)))
(def strip (move (box 0.4 20 0.4) 0 0 5))
; expect-status: neighbouring-walls fails
(require neighbouring-walls (empty (overhangs (min neighbours strip) 0 45 1)) (accuracy 0.05))

(def blank (move (box 24 3 24) 0 0 12))
(def slot (max (slab-x -9 9) (- (abs (- z (+ 12 (* x (tan (deg->rad 40)))))) 0.875)))
; expect-status: oblique-slot fails
(require oblique-slot (empty (overhangs (rotate-z (subtract blank slot) 23) 0 45 1)) (accuracy 0.05))

; Neither a tiny air gap above the bed nor a thin unresolved layer is support.
; expect-status: floating fails
(require floating (empty (overhangs (move floor 0 0 5) 0 45 1)) (accuracy 0.05))
; expect-status: floating-near-bed fails
(require floating-near-bed (empty (overhangs (move floor 0 0 0.001) 0 45 1)) (accuracy 0.05))
; expect-status: foil fails
(require foil (empty (overhangs foil 0 45 1)) (accuracy 0.5))

; The vertical ray from (0,0,3) has air below the cap before it meets
; the post at 2.5. The +x descending ray is solid across that gap. It
; runs out of brace at 2.4, by which height the centre is already solid.
(def post (extrude-z (rect-2d 0.2 0.2) 0 2.5))
(def cap (extrude-z (rect-2d 0.4 0.4) 2.8 3.2))
(def brace (max (- (abs (- (+ x z) 3)) 0.06) (slab-y -0.08 0.08) (slab-z 2.4 3.1)))
(def unsupported (min post cap brace))
(def region (move (box 0.02 0.02 0.02) 0 0 3))
; expect-status: staggered unknown
(require staggered (empty (max region (overhangs unsupported 0 45 1))) (accuracy 0.001))

; Two halves of the same block meet on a tilted plane. Its height ranges
; across a quarter millimetre, and both halves use the same parameter.
(param seam-height 3 (tol -0.25 0.25))
(def seam-stock (max (rect-2d 12 12) (slab-z 0 7)))
(def tilted-face (- (+ z (* 0.125 x)) seam-height))
(def halves (min (max seam-stock tilted-face) (max seam-stock (neg tilted-face))))
(require tilted-joint (empty (overhangs halves 0 45 1)))

; Equal ranges on two independent dimensions can leave a real gap.
(param second-height 3 (tol -0.25 0.25))
; expect-status: independent-faces fails
(require independent-faces (empty (overhangs (min (max seam-stock (- z seam-height)) (max seam-stock (- second-height z))) 0 45 1)))

; A gap of one representable step at 5 mm is smaller than the witness
; search can resolve. The support proof must still see the distinct faces.
(def footprint (rect-2d 12 12))
(def lower (max footprint (slab-z 0 5)))
(fn cap-at (h) (max footprint (slab-z h 6)))
; expect-status: one-step-gap unknown
(require one-step-gap (empty (overhangs (min lower (cap-at 5.000000000000001)) 0 45 1)))
; expect-status: bearing-gap fails
(require bearing-gap (empty (overhangs (min lower (cap-at 5.15)) 0 45 1)))

; Both signs of a shared expression still exclude a volume where it is zero.
(def plateau (max (- z 5) (min 0 (- z 4))))
; expect-status: zero-plateau fails
(require zero-plateau (empty (overhangs (min (max seam-stock plateau) (max seam-stock (neg plateau))) 0 45 1)))

; A rectangular rib crosses a cut through the centre of its bearing.
(def missing (max (slab-x -2 2) (slab-y -2 2) (slab-z 0 8)))
(def rib (max (slab-x -5 5) (slab-y -3 3) (slab-z 4.8 7)))
; expect-status: incomplete-footprint fails
(require incomplete-footprint (empty (overhangs (min (subtract lower missing) rib) 0 45 1)))

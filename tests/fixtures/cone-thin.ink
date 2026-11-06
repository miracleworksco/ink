; Long thin cones (radius 1 at z = 0, 1 + 10^n at z = 1, n = 3..9): within
; the box each is the whole box's slab 0 < z < 1, so a ball at z = -0.55 of
; radius 0.1 is 0.45 from it. Their wedges' edges are within 1e-n of
; parallel; normals built through sqrt(1 - cs^2) once had slope 1.34 at
; 5e7 under a promise of 1, and a clearance of 1 held.
; Each 1 fails; each 0.4 holds, so the promise still speeds the proof.
(import "std")
(def b (move (sphere 0.1) 0 0 -0.55))
(def shifted (max (cone-z 1 50000001 0 1) (box 2 2 2)))
; expect-status: shifted-gap fails
(assert shifted-gap (clearance shifted b 1.0) (accuracy 0.01))
(assert shifted-in (subset (probe 0 0 0.46) shifted))
(def c3 (max (cone-z 1 1001 0 1) (box 2 2 2)))
; expect-status: far-3 fails
(assert far-3 (clearance c3 b 1.0) (accuracy 0.01))
(assert near-3 (clearance c3 b 0.4) (accuracy 0.01))
(def c4 (max (cone-z 1 10001 0 1) (box 2 2 2)))
; expect-status: far-4 fails
(assert far-4 (clearance c4 b 1.0) (accuracy 0.01))
(assert near-4 (clearance c4 b 0.4) (accuracy 0.01))
(def c5 (max (cone-z 1 100001 0 1) (box 2 2 2)))
; expect-status: far-5 fails
(assert far-5 (clearance c5 b 1.0) (accuracy 0.01))
(assert near-5 (clearance c5 b 0.4) (accuracy 0.01))
(def c6 (max (cone-z 1 1000001 0 1) (box 2 2 2)))
; expect-status: far-6 fails
(assert far-6 (clearance c6 b 1.0) (accuracy 0.01))
(assert near-6 (clearance c6 b 0.4) (accuracy 0.01))
(def c7 (max (cone-z 1 10000001 0 1) (box 2 2 2)))
; expect-status: far-7 fails
(assert far-7 (clearance c7 b 1.0) (accuracy 0.01))
(assert near-7 (clearance c7 b 0.4) (accuracy 0.01))
(def c8 (max (cone-z 1 100000001 0 1) (box 2 2 2)))
; expect-status: far-8 fails
(assert far-8 (clearance c8 b 1.0) (accuracy 0.01))
(assert near-8 (clearance c8 b 0.4) (accuracy 0.01))
(def c9 (max (cone-z 1 1000000001 0 1) (box 2 2 2)))
; expect-status: far-9 fails
(assert far-9 (clearance c9 b 1.0) (accuracy 0.01))
(assert near-9 (clearance c9 b 0.4) (accuracy 0.01))

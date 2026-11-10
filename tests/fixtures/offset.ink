; Growing a shape by d is everything within d of it, however the shape is
; written; where the field cannot give that, ink check says so.
(import "std")

; One cube written two ways grows alike: rounded edges, 0.5 out.
(body box-grown (shape (offset-shape (box 2 2 2) 0.5)))
(body slabs-grown (shape (offset-shape (max (slab-x -1 1) (slab-y -1 1) (slab-z -1 1)) 0.5)))
(assert edge-rounded (disjoint (probe 1.4 1.4 0) (shape-of slabs-grown)))
(assert face-grown (subset (probe 1.45 0 0) (shape-of slabs-grown)))
; The field less d is the field's own growth: the max of slabs mitres.
(body slabs-mitred (shape (- (max (slab-x -1 1) (slab-y -1 1) (slab-z -1 1)) 0.5)))
(assert mitred (subset (probe 1.4 1.4 0) (shape-of slabs-mitred)))

; A shell of a box of slabs is rounded outside, sharp inside.
(body slabs-shell (shape (shell (max (slab-x -1 1) (slab-y -1 1) (slab-z -1 1)) 0.2)))
(assert shell-rounded (disjoint (probe 1.09 1.09 0) (shape-of slabs-shell)))
(assert shell-sharp-inside (subset (probe 0.92 0.92 0) (shape-of slabs-shell)))

; A pocket's mouth is a 0 its shape keeps out: the field is not the
; distance outside, so a shell closes the mouth, and ink check says so.
; expect-warning: offset.ink:23): (offset-shape f h): the field is not shown to be the distance outside its shape
(def pocket (subtract (move (box 40 40 10) 20 20 5) (move (box 20 20 5) 20 20 7.5)))
(body pocket-shell (shape (shell pocket 1)))
; expect-status: mouth-is-air fails
(assert mouth-is-air (disjoint (probe 20 20 10) (shape-of pocket-shell)))
; Shrunk, a subtract is exact.
(body pocket-shrunk (shape (offset-shape pocket -1)))

; A region read as 0 is air: grown, the field fills it.
; expect-warning: offset.ink:31 (body zeros ...): (offset-shape (min 0 (sphere 1)) 0.5): the field is not shown to be the distance outside its shape
(body zeros (shape (max (box 4 4 4) (offset-shape (min 0 (sphere 1)) 0.5))))

; A coefficient read from a range is a real number, not its float there:
; (1 + p) - 1 is p, so these planes are a wedge, not a slab.
(def p (ranged 1e-16 2e-16))
; expect-warning: offset.ink:37 (body wedge ...): (offset-shape (max
(body wedge (shape (max (move (box 10 10 10) 0 4000000 0) (offset-shape (max (- x 1e-10) (- (+ (neg x) (* (- (+ 1 p) 1) y)) 1e-10)) 1))))

; A value the shape reads as real arithmetic (here exactly 0, its float -1)
; is not trusted as its float: min(x, r) is not shown the distance.
(param s 1)
(def r (- (warp (- (+ 1e16 x) 1e16) s 0 0) s))
; expect-warning: offset.ink:44 (body real-zero ...): (offset-shape (min x r) 1): the field is not shown to be the distance outside its shape
(body real-zero (shape (max (offset-shape (min x r) 1) (move (box 2 2 2) 1000 0 0))))
; An amount of real arithmetic is bounded, not its float: (1e16 + q) - 1e16
; is q, so this is no zero offset of a field that is no distance.
(def q (ranged 1 2))
; expect-warning: offset.ink:49 (body real-amount ...) (at q = 1): (offset-shape (* 2 x) (- (+ 10000000000000000 q) 10000000000000000)): the field is not shown to be the distance
(body real-amount (shape (max (box 4 4 4) (offset-shape (* 2 x) (- (+ 1e16 q) 1e16)))))
; Bounded over its whole range, read as real numbers: 2 - ((1e16 + q) - 1e16)
; is 2 - q, from 0 to 1, not the float 2 at q = 1, so no radius is promised.
; expect-warning: offset.ink:53 (body ranged-bound ...) (at q = 1): (offset-shape (* 2 x) (- 2 (- (+ 10000000000000000 q) 10000000000000000))): the field is not shown to be the distance to its shape, nor its slope bounded (or the amount's sign is not shown)
(body ranged-bound (shape (max (box 4 4 4) (offset-shape (* 2 x) (- 2 (- (+ 1e16 q) 1e16))))))
; Nor a side: 0.5 - ((1e16 + q1) - 1e16) is -0.5, a shrinking, which a
; union is not shown to give exactly, though its float is +0.5.
(def q1 (ranged 1 1))
; expect-warning: offset.ink:58 (body ranged-sign ...) (at q1 = 1): (offset-shape (min (move (sphere 2) -1 0 0) (move (sphere 2) 1 0 0)) (- 0.5 (- (+ 10000000000000000 q1) 10000000000000000))): the field is not shown to be the distance to its shape, nor its slope bounded (or the amount's sign is not shown)
(body ranged-sign (shape (offset-shape (min (move (sphere 2) -1 0 0) (move (sphere 2) 1 0 0)) (- 0.5 (- (+ 1e16 q1) 1e16)))))
; A solved value over a range is bounded as proofs bound it: s = q1 is 1
; in real numbers, so 0.5 - ((1e16 + s) - 1e16) is -0.5, a shrinking.
(solve follows (unknowns (s1 1)) (== s1 q1))
; expect-warning: offset.ink:63 (body solved-sign ...) (at s1 = 1): (offset-shape (min (move (sphere 2) -1 0 0) (move (sphere 2) 1 0 0)) (- 0.5 (- (+ 10000000000000000 s1) 10000000000000000))): the field is not shown to be the distance to its shape, nor its slope bounded
(body solved-sign (shape (offset-shape (min (move (sphere 2) -1 0 0) (move (sphere 2) 1 0 0)) (- 0.5 (- (+ 1e16 s1) 1e16)))))
; and its slope: 1 + ((1e16 + s) - 1e16) is 2, not the float 1.
; expect-warning: offset.ink:66 (body solved-slope ...): (offset-shape (* (+ 1 (- (+ 10000000000000000 s1) 10000000000000000)) x) 1): the field is not shown to be the distance outside its shape; its slope is at most 3
(body solved-slope (shape (max (box 4 4 4) (offset-shape (* (+ 1 (- (+ 1e16 s1) 1e16)) x) 1))))

; Samples of a distance (redistance's, an imported mesh's) are no
; distance between them: a shell of them is approximate, by how much
; ink check says. Scaled, they are not even that.
(def ball (redistance (* 3 (sphere 2)) (res 0.1)))
; expect-warning: offset.ink:74): (offset-shape f h): an offset of sampled data is approximate
; expect-warning: offset.ink:74): (offset-shape f (neg h)): an offset of sampled data is approximate
(body redistanced-shell (shape (move (shell (move ball 0 0 0.5) 0.4) 0 30 0)))
; expect-warning: offset.ink:77): (offset-shape f h): the field is not shown to be the distance
; expect-warning: offset.ink:77): (offset-shape f (neg h)): the field is not shown to be the distance
(body scaled-samples (shape (move (shell (warp ball (* 2 x) y z) 0.4) 0 40 0)))

; Benchmark 21: gyroid sheet lattice block: density, connectivity, stiffness.
;
; Task: a 36 mm cube filled with a gyroid sheet lattice, cell 12 mm, wall
; 1.5 mm, between two 2 mm solid face plates (z). Check the relative
; density against the gyroid's surface area, that the lattice is one
; connected solid with no trapped voids, and its effective compressive
; stiffness.
;
; Reference values:
;   - exact gyroid minimal surface: area per cubic cell = 3.0937 a^2, so a
;     thin sheet of thickness t fills phi = 3.0937 t / a = 0.3867 of the
;     lattice region (36 x 36 x 32);
;   - std `gyroid` is the level set sin cos + ... = 0 divided by its
;     greatest gradient, sqrt(3) k: never over the distance, so (shell
;     (gyroid a) t) is at least t thick, thicker where the gradient is less
;     (phi 0.441 at res 0.3); the redistanced version should match the
;     formula closely;
;   - stiffness: sheet TPMS are near stretch-dominated; published E*/Es
;     for sheet gyroids at phi ~ 0.35-0.4 are about 0.15-0.3 (e.g. Al-Ketan
;     et al. 2018). Loose band, used as a sanity check only.
;
; Pass criteria (ink mass --res 0.3, ink components, ink solve):
;   - redistanced lattice volume fraction within 5% of 0.3867;
;   - one solid, zero cavities (the two sheet labyrinths stay open at the
;     sides);
;   - compression: equilibrium, and E*/Es in [0.1, 0.35].

; The lattices alone stand beside the block they fill.

(import "std")

(param cell 12)
(param wall 1.5)
(def region (box 36 36 32))
(def g-naive (gyroid cell))
; redistance needs a bounded field: sample the gyroid's solid half inside a
; box 3 mm larger than the region, so within the region (and a wall from it)
; the distance is to the gyroid surface only (redistance refuses a field
; that is not bounded).
(def g-exact (redistance (max (box 42 42 38) (gyroid cell)) (res 0.3)))
(def lattice-naive (lattice-fill region g-naive wall))
(def lattice (lattice-fill region g-exact wall))
(def plates (max (box 36 36 36) (- 16 (abs z))))
(body block (shape (min lattice plates)) (attr material "nylon-pa12") (attr color "#D0D0C8"))
(body lattice-only (shape (move lattice 60 0 0)))
(body lattice-naive-only (shape (move lattice-naive 120 0 0)))

; Compression between the plates: bottom plate fixed, 100 N on the top.
; E* = (F / A) / (delta / H), H = 32: the lattice only (the plates are stiff).
(study squash (shape (shape-of block)) (elastic (material nylon-pa12))
  (fix bottom (- z -17))
  (load press (force 0 0 -100N) (region (- 17 z)))
  (res 0.5))

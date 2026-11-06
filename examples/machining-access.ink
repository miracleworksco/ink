; A machined aluminium motor mount, checked for what a 3-axis mill reaches
; from its setups, and for its walls and holes (docs/IR.md, Manufacturing
; checks). It is cut from the top (+z: a pocket, counterbored holes), the
; bottom (-z: the holes' exits) and one side (+x: a cross hole into the
; pocket). A control adds a T-slot under the pocket's rim and fails: no
; straight approach reaches under the rim.
;
; The stock is the block the mount is cut from. An access claim says every
; point of the stock's air, 0.5 mm clear of the part, sees out along a
; setup's direction; std's hidden looks along it, ranged over every
; distance to past the part's far side. One claim per setup, over the
; region it cuts: the cross hole from the side, the rest from top or
; bottom. Together they cover the stock's air; from top and bottom alone
; the cross hole is hidden. They read only where material is, so a sealed
; void is hidden however its field is written.
(import "std")

(def top (vec 0 0 1))
(def bottom (vec 0 0 -1))
(def side (vec 1 0 0))

(def block (move (box 60 40 16) 0 0 8))
(def pocket (extrude-z (rounded-rect-2d 30 20 3) 6 20))                ; floor at z = 6
(def bolt (min (cylinder-z 2.2 40) (extrude-z (circle-2d 4) 12 20)))  ; M4 clearance, counterbored 4 deep
(def bolts (at-frames bolt (point-frames (rect-points 24 13))))   ; on a 48 x 26 pattern
(def cross (move (cylinder-x 1.6 30) 22 0 11))                          ; from the +x face into the pocket
(def mount (subtract block pocket bolts cross))
(def t-slot (subtract mount (move (box 36 20 3) 0 0 7.5)))            ; 3 mm under the pocket's rim, each end
(body mount (shape mount) (attr material "aluminum-6061"))
(body mount-t-slot (shape t-slot) (attr material "aluminum-6061"))

(def air (clear-air mount block 0.5))
(assert top-and-bottom-reach (empty (hidden mount block (max air (neg cross)) (vec top bottom))))
(assert side-reaches-cross (empty (hidden mount block (max air cross) (vec side))))
; expect-status: cross-from-top-and-bottom fails
(assert cross-from-top-and-bottom (empty (hidden mount block air (vec top bottom))))
; expect-status: t-slot-reaches fails
(assert t-slot-reaches (empty (hidden t-slot block (max (clear-air t-slot block 0.5) (neg cross)) (vec top bottom))))
; A cube with a sealed void, its field written as a box less a box and as
; that capped at 0.2 (the same shape): the void is hidden either way.
(def sealed (subtract (box 20 20 20) (box 12 12 12)))
(def capped (min sealed 0.2))
; expect-status: sealed-reaches fails
(assert sealed-reaches (empty (hidden sealed (box 20 20 20) (clear-air sealed (box 20 20 20) 0.5) (vec top bottom))))
; expect-status: capped-reaches fails
(assert capped-reaches (empty (hidden capped (box 20 20 20) (clear-air capped (box 20 20 20) 0.5) (vec top bottom))))
; The T-slot's field times 30, the same shape (a band of its values would
; be 30 times thinner, too thin to hold a hidden point).
(def steep (* 30 t-slot))
; expect-status: steep-t-slot-reaches fails
(assert steep-t-slot-reaches (empty (hidden steep block (max (clear-air steep block 0.5) (neg cross)) (vec top bottom))))
; The T-slot with its height toleranced, 3 to 3.2 mm: checked over the
; band, it fails as at 3.
(param slot-h 3 (tol 0 0.2))
(def t-slot-tol (subtract mount (move (box 36 20 slot-h) 0 0 (+ 6 (* 0.5 slot-h)))))
; expect-status: toleranced-t-slot-reaches fails
(assert toleranced-t-slot-reaches (empty (hidden t-slot-tol block (max (clear-air t-slot-tol block 0.5) (neg cross)) (vec top bottom))))
; The stock written without its move (z -8 to 8, the part 0 to 16): the
; rim over the T-slot lies outside it, and still blocks.
; expect-status: unmoved-stock-t-slot-reaches fails
(assert unmoved-stock-t-slot-reaches (empty (hidden t-slot (box 60 40 16) (max (clear-air t-slot (box 60 40 16) 0.5) (neg cross)) (vec top bottom))))

; Walls at least 0.8 mm (Protolabs Network's recommendation for metals),
; holes and slots at least 2.5 mm across (their smallest cavity for
; standard tools): the air around the part, 5 mm past it, is nowhere
; thinner.
(assert walls (min-thickness mount > 0.8))
(assert holes (min-thickness (air-around mount 5) > 2.5))

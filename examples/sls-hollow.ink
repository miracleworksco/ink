; An SLS-printed knob, hollowed to save powder and weight, checked against
; Protolabs' SLS design rules (docs/IR.md, Manufacturing checks): walls at
; least 0.7 mm, and escape holes at least 5 mm across to empty the
; unsintered powder inside. Sealed, the hollow is an enclosed void (counted
; from a mesh, an estimate). Two 6 mm holes across it empty it. 3 mm holes
; empty it too, but are under the 5 mm the powder needs.
(import "std")

(def knob (extrude-z (circle-2d 15) 0 20))
(def inside (extrude-z (circle-2d 13) 2 18))                  ; walls 2 thick
(def sealed (subtract knob inside))
(fn drained (d) (subtract sealed (move (cylinder-x (* 0.5 d) 40) 0 0 10)))
(def holed (drained 6))
(def small (drained 3))
(body knob-sealed (shape sealed) (attr material "nylon-pa12"))
(body knob (shape holed) (attr material "nylon-pa12"))
(body knob-small-holes (shape small) (attr material "nylon-pa12"))

; expect-status: sealed-empties fails estimate
(assert sealed-empties (cavities sealed < 1))
(assert empties (cavities holed < 1))
(assert walls (min-thickness holed > 0.7))
; The air around the part, 6 mm past it, is nowhere thinner than the
; escape holes must be.
(assert escape-holes (min-thickness (air-around holed 6) > 5))
; expect-status: small-escape-holes fails estimate
(assert small-escape-holes (min-thickness (air-around small 6) > 5))

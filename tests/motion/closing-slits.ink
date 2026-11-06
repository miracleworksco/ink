; Three plates 0.06 thick, pitch 0.059 + 0.004 t apart, inside a solid:
; at t = 0 they overlap and make one slab 0.178 thick (its middle 0.089
; deep in the solid), from t = 0.25 three plates no deeper than 0.03. Each
; moves 0.004 at most, but how deep their union is jumps where the gaps
; close: a step proven at its middle with that travel taken out of the
; tolerance (once: proven with one step) says nothing of t = 0. Proven
; over every moment of a step at once, it is not proven: unknown, alone
; (slits) and with a block of the same shape beside them against a ball
; (slits-in-ball).
(import "std")
(t-span 0 1)
(def pitch (+ 0.059 (* 0.004 t)))
(def plates (min
  (place (box 0.06 1 1) (at (move (neg pitch) 0 0)))
  (box 0.06 1 1)
  (place (box 0.06 1 1) (at (move pitch 0 0)))))
(assert slits (disjoint plates (box 1 1 1)) (accuracy 0.05))
(assert slits-in-ball (disjoint (min plates (move (box 0.2 0.2 0.2) 0.4 0.4 0.4)) (sphere 0.5)) (accuracy 0.05))
; expect-status: slits unknown
; expect-status: slits-in-ball unknown

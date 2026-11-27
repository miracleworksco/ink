; examples/arm.ink, with the arm carrying its motor's rating as a catalog
; part does (the kernel's ink-ratings), and the design asking whether it
; stays within it: the torque its joint needs, at most the 150 N mm its
; effort-limits state. Nothing states that torque: the solver gives it.
(import "std")
(gravity 0 0 -9.81m/s2)
(t-span 0 2s)
(param limit 150N*mm)
(def rated (read-json "rated-arm.json"))
(body arm (shape (move (box 200 20 10) 100 0 0)) (attr material "aluminum-6061") (attr ratings (get rated "ratings")))
(assembly (joint (angle 0))
  (instance arm arm (at) (turn y angle) (effort-limits -150N*mm 150N*mm)))
(assembly root
  (instance j (joint (* 45deg (sin (* 180deg/s t))))))
(require motor-ok (within-ratings root))

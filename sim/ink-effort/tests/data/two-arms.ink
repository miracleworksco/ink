; rated-arm.ink's arm twice, in two assemblies that each move their own
; way: `swing` through level (its joint's torque peaks there), `raised`
; within 5 degrees of 80 up (a sixth of it at most). Each requirement is
; answered over its own assembly's motion.
(import "std")
(gravity 0 0 -9.81m/s2)
(t-span 0 2s)
(def rated (read-json "rated-arm.json"))
(body arm (shape (move (box 200 20 10) 100 0 0)) (attr material "aluminum-6061") (attr ratings (get rated "ratings")))
(assembly (joint (angle 0))
  (instance arm arm (at) (turn y angle) (effort-limits -150N*mm 150N*mm)))
(assembly swing
  (instance j (joint (* 45 (sin (* 180deg/s t))))))
(assembly raised
  (instance j (joint (+ 80 (* 5 (sin (* 180deg/s t)))))))
(require swing-ok (within-ratings swing))
(require raised-ok (within-ratings raised))

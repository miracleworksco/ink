; arm.ink's motor and arm (root/j), beside two more of them: one whose arm has no material and so no mass
; (root/k), and one whose arm carries such a body at its tip (root/c). No body of k's or c's is on j's arm.
(import "std")
(gravity 0 0 -9.81m/s2)                      ; m/s^2, z up
(t-span 0 2s)                             ; s

(body arm (shape (move (box 200 20 10) 100 0 0)) (attr material "aluminum-6061"))
(body bare (shape (move (box 200 20 10) 100 0 0)))
(assembly pair
  (instance a arm)
  (instance tip bare (at (move 200 0 0))))

(assembly (joint (angle 0))              ; angle: deg, an argument the motion passes
  (instance arm arm (at) (turn y angle) (effort-limits -150N*mm 150N*mm)))   ; the motor: N mm
(assembly (loose (angle 0))
  (instance arm bare (at) (turn y angle) (effort-limits -150N*mm 150N*mm)))
(assembly (carrying (angle 0))
  (instance arm pair (at) (turn y angle) (effort-limits -150N*mm 150N*mm)))
(assembly root
  (instance j (joint (* 45deg (sin (* 180deg/s t)))))
  (instance k (loose (* 30deg (sin (* 180deg/s t)))) (at (move 0 100 0)))
  (instance c (carrying (* 20deg (sin (* 180deg/s t)))) (at (move 0 200 0))))

; An outside read's case: examples/arm.ink, its arm with a frame `other-axis` and two ratings,
; one reading the torque on a frame the arm does not have (`missing`).
(import "std")
(gravity 0 0 -9.81m/s2)                      ; m/s^2, z up
(t-span 0 2s)                             ; s

(body arm (shape (move (box 200 20 10) 100 0 0)) (attr material "aluminum-6061") (frame other-axis (axes (origin 0 0 0) (x 0 1 0) (z 0 0 1))) (attr ratings (read-json "frame-torque.json")))

(assembly (joint (angle 0))              ; angle: deg, an argument the motion passes
  (instance arm arm (at) (turn y angle) (effort-limits -150N*mm 150N*mm)))   ; the motor: N mm
(assembly root
  (instance j (joint (* 45deg (sin (* 180deg/s t))))))
(require ok (within-ratings root))

; A motor swings an aluminium arm 45 degrees each way about a level shaft.
; What torque must it give? By hand: the arm (200 x 20 x 10 mm, 108 g) held
; level takes m g L / 2 = 0.108 x 9.81 x 100 = 105.9 N mm.
; ink effort examples/arm.ink: 105.9 as it passes level, and its peak,
; 106.9, just below level, where it also slows the arm's swing. ink require
; checks the motor's limits (+-150 N mm) along the whole motion: they hold.
(import "std")
(gravity 0 0 -9.81m/s2)                      ; m/s^2, z up
(t-span 0 2s)                             ; s

(body arm (shape (move (box 200 20 10) 100 0 0)) (attr material "aluminum-6061"))

(assembly (joint (angle 0))              ; angle: deg, an argument the motion passes
  (instance arm arm (at) (turn y angle) (effort-limits -150N*mm 150N*mm)))   ; the motor: N mm
(assembly root
  (instance j (joint (* 45deg (sin (* 180deg/s t))))))

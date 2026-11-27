; A stepper's moving pieces (17HS19-2004S1's numbers): (drive ANGLE) turns
; its rotor (82 g cm2) in its still stator, held at its holding stiffness
; (0.59 N m x 50 rotor teeth), its effort limits its 24 V pull-out curve,
; N mm at rpm. Below the curve's first speed its first torque; past its
; last speed (1200 rpm) the limits cross, so no effort lies within them.
(import "std")
(def curve (vec (vec 100 412) (vec 200 406.1) (vec 300 406.2) (vec 400 422.1) (vec 500 402) (vec 600 362.2)
                (vec 700 324) (vec 800 286) (vec 900 254.1) (vec 1000 232.1) (vec 1100 200) (vec 1200 172)))
(fn pull-out (rate) (let ((s (/ (+ 0deg/s (abs rate)) 1rpm))) (select (- 1200 s) (neg inf) (curve-at curve s))))
(body stator (shape (subtract (move (box 42 42 48) 0 0 -24) (cylinder-z 3 200))) (attr density 7.85g/cm3) (frame mount (at)))
(body rotor (shape (max (move (cylinder-z 2.5 64) 0 0 -8) (- y 2))) (attr density 7.85g/cm3) (attr spin-inertia 8200g*mm2)
  (frame shaft (at (move 0 0 24))))
(assembly (drive (angle 0))
  (instance stator stator)
  (instance rotor rotor (frame-of stator) (turn z angle) (stiffness 29500N*mm)
    (effort-limits (* 1N*mm (neg (pull-out (diff angle t)))) (* 1N*mm (pull-out (diff angle t)))))
  (frame mount (frame-of stator mount))
  (frame shaft (frame-of rotor shaft)))

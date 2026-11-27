; An outside read's case: the motor's effort limits cross from t = 1/32, so its effort is
; refused from there, while its speed, which a rating of its own limits to 1 deg/s, is known.
(import "std")
(gravity 0 0 0)
(t-span 0 1s)
(body arm
  (shape (move (box 200 20 10) 100 0 0))
  (attr material "aluminum-6061")
  (attr ratings (read-json "speed-crossed.json"))
  (frame shaft (axes (origin 0 0 0) (x 1 0 0) (z 0 1 0))))   ; the shaft, along the joint's y
(assembly (drive (angle 0deg))
  (instance arm arm (at) (turn y angle) (effort-limits (* 1N*mm (- (* 3200 (/ t 1s)) 100)) (* 1N*mm (- 100 (* 3200 (/ t 1s)))))))
(assembly root
  (instance j (drive (* 16deg/s2 t t))))
(require ok (within-ratings root))

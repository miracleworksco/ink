; An outside read's case: a made-up curve with a narrow dip at 2 deg/s, which the speed passes
; at t = 0.0625 s, between the moments the run reads.
(import "std")
(gravity 0 0 0)
(t-span 0 1s)
(body arm
  (shape (move (box 200 20 10) 100 0 0))
  (attr material "aluminum-6061")
  (attr ratings (read-json "speed-notch.json"))
  (frame shaft (axes (origin 0 0 0) (x 1 0 0) (z 0 1 0))))   ; the shaft, along the joint's y
(assembly (drive (angle 0deg))
  (instance arm arm (at) (turn y angle) (effort-limits -100N*mm 100N*mm)))
(assembly root
  (instance j (drive (* 16deg/s2 t t))))
(require ok (within-ratings root))

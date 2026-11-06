; Four Brotherhobby Avenger 2806.5 motors placed by frames at the ends of an
; X frame's arms, 150 mm from its centre. Each arm's pad has a frame on its
; top face, z up; each motor meets it by its mount frame (its base, z down),
; face to face.
;
; (at STEP ...) applies its steps in order, about the assembly's axes. A pad
; at radius L and angle a is (at (move L 0 4) (rot-z a)): out along x, then
; turned about the centre. (at (rot-z a) (move L 0 4)) turns the pad in
; place, then moves it along x: all four land on +x, where there is no arm.
(import "std")

; A Brotherhobby Avenger 2806.5 (its drawing: 32.9 across, 19.1 tall, a 5 mm
; shaft 14 above; four M3 taps on its base, 9.5 from its axis on its x and
; y): a base 2.4 thick, the bell above, the shaft. mount: its base, z down
; out of it; tap-I: tap I's mouth, z along it, out of the motor.
(def taps (vec (vec 9.5 0) (vec 0 9.5) (vec -9.5 0) (vec 0 -9.5)))
(body motor
  (shape (subtract (min (extrude-z (circle-2d 15.8) 0 2.4) (extrude-z (circle-2d 16.45) 2.4 19.1)
                        (extrude-z (circle-2d 2.5) 19.1 33.1))
                   (reduce min (map (lambda (q) (move (extrude-z (circle-2d 1.5) -1 5) (nth q 0) (nth q 1) 0)) taps))))
  (attr mass 41) (attr color (vec .16 .17 .18))
  (frame mount (at (rot-x 180)))
  (for i (range 4) (frame tap-$i (at (rot-x 180) (move (nth (nth taps i) 0) (nth (nth taps i) 1) 0)))))

(param L 150)                                ; centre to motor axis
(def th 4)                                   ; plate thickness
(fn arm-angle (i) (+ 45 (* 90 i)))
; One arm along +x: a bar, a round pad, the motor's four M3 clearance holes
; (its taps 9.5 from its axis on its x and y) and a shaft hole.
(def holes (min (cylinder-z 4 20)
                  (min (move (cylinder-z 1.6 20) 9.5 0 0) (move (cylinder-z 1.6 20) -9.5 0 0)
                         (move (cylinder-z 1.6 20) 0 9.5 0) (move (cylinder-z 1.6 20) 0 -9.5 0))))
(def arm (subtract (min (move (box L 16 th) (/ L 2) 0 (/ th 2)) (move (cylinder-z 18 th) L 0 (/ th 2)))
                   (move holes L 0 0)))
(body frame
  (shape (min (move (rounded-box 70 70 th 6) 0 0 (/ th 2))
                (min (rotate-z arm (arm-angle 0)) (rotate-z arm (arm-angle 1))
                       (rotate-z arm (arm-angle 2)) (rotate-z arm (arm-angle 3)))))
  (attr color "#3A3A3A")
  (for i (range 4) (frame pad-$i (at (move L 0 th) (rot-z (arm-angle i))))))

(assembly root
  (instance frame frame)
  (for i (range 4) (instance motor-$i motor (frame-of frame pad-$i) mount)))

; Each motor stands on its pad; its tap-0 screw hole is over a hole in the
; pad (a point 2 mm into the pad along the tap's axis is clear of it); the
; motors are spread round the centre, adjacent ones 2 L sin 45 apart.
(for i (range 4)
  (assert seated-$i (touches (shape-of root/motor-$i) (shape-of root/frame)))
  (assert screw-hole-$i (disjoint (shape-of root/frame) (place (probe 0 0 2) (frame-of root/motor-$i tap-0)))))
(assert spread (clearance (shape-of root/motor-0) (shape-of root/motor-1) 150))

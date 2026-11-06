; A NEO motor spinning up, driving a GT2 pulley on its shaft and a 15:45
; gear pair. The pinion rides the motor's shaft. The wheel sits and turns as
; std's mesh-distance and mesh-turn give, and the two are a gear mesh
; (gear-mesh): interfere checks their teeth over the motor's angle, assert
; their centre distance and axes. t is seconds: each drive spins up from
; rest to its output speed over the span.
;
; The NEO's record gives its free speed, 5676 rpm, which its rotor's joint
; states as a speed limit. 2000 rpm at the output needs 6000 at the motor,
; so that claim fails when the spin-up passes 5676. A limit no part's record
; gives is the design's own, a value claim: the belt, 40 mm a turn on the
; 20-tooth pulley, stays under 4000 mm/s (the slow drive ends at 4500 rpm,
; 3000 mm/s). A stepper's rotor has no speed limit: its torque falls with
; speed, by its driver and supply, so it has no single top speed.
; expect-status: speed-limits:root/fast/m/rotor fails

; The two drives side by side.

(import "std")

; A REV NEO (REV-21-1650), from its drawing and datasheet: a 60 mm body 58.3
; long behind its mounting face (z = 0), a 19.05 pilot 3.5 tall, its 8 mm
; shaft to 35 out (shaft: its tip); 425 g; at 12 V 5676 rpm free. (neo
; ANGLE) turns the rotor, within its free speed either way.
(body housing
  (shape (min (extrude-z (circle-2d 30) -58.3 0) (extrude-z (circle-2d 9.525) 0 3.5)))
  (attr color (vec .1 .1 .11)) (attr metallic .4))
(body rotor (shape (extrude-z (circle-2d 4) 3.5 35)) (attr color (vec .66 .68 .71)) (attr metallic .85) (attr material "steel")
  (frame shaft (at (move 0 0 35))))
(def free-speed 5676rpm)                     ; deg/s
(assembly (neo (angle 0)) (attr mass 425)
  (instance housing housing)
  (instance rotor rotor (frame-of housing) (turn z angle) (speed-limits (neg free-speed) free-speed))
  (frame shaft (frame-of rotor shaft)))
; A GT2 pulley, 20 teeth (std's GT2 groove, gt-pulley-teeth-2d), for a 6 mm
; belt: 16 long and 16 across its flange and hub, its belt's centre plane at
; z = 0; hub: its hub's end, 11.5 below, z out of it.
(fn gt2-pulley (bore)
  (subtract (min (extrude-z (gt-pulley-teeth-2d 20 2 0) -3.5 3.5)
                 (extrude-z (circle-2d 8) -11.5 -3.5) (extrude-z (circle-2d 8) 3.5 4.5))
            (cylinder-z (* 0.5 bore) 40)))
(body pulley (shape (gt2-pulley 8)) (attr color (vec .8 .81 .83)) (attr metallic 1) (attr material "aluminum-6061")
  (frame hub (at (rot-x 180) (move 0 0 -11.5))))

(t-span 0 1)
(def zp 15)                                     ; pinion teeth
(def zw 45)                                     ; wheel teeth
(def module 1)
(def backlash 0.1)                              ; each gear's share
(def face 8)

(fn gear (teeth bore)
  (subtract (extrude-z (involute-gear-2d teeth module 20 backlash) 0 face) (cylinder-z (* 0.5 bore) 20)))
(body pinion (shape (gear zp 8)) (attr color (vec 0.7 0.7 0.72)))
(body wheel (shape (gear zw 8)) (attr color (vec 0.75 0.6 0.3)))

; Degrees turned spinning up from rest to r rpm in the span's 1 s.
(fn spin-up (r) (* 0.5 (* r 1rpm) (square t)))
; The wheel is the output; the motor and the pinion on its shaft turn
; zw / zp times as fast, the other way.
(fn motor-turn (output-rpm) (* (neg (/ zw zp)) (spin-up output-rpm)))
(assembly (gearbox (output-rpm 0))
  (instance m (neo (motor-turn output-rpm)))
  (instance belt pulley (at (move 0 0 -16) (frame-of m shaft)) hub)   ; rides the shaft
  (instance p pinion (at (move 0 0 -30) (frame-of m shaft)))
  (instance w wheel (at (move (mesh-distance zp zw module) 0 5))           ; level with the pinion
    (turn z (mesh-turn zp zw 0 (motor-turn output-rpm))))
  (relate p w (gear-mesh zp zw module 20 (* 2 backlash))))

(assembly root
  (instance slow (gearbox 1500))
  (instance fast (gearbox 2000) (at (move 0 90 0))))
(def belt-per-turn (* 20 2))                    ; teeth x 2 mm pitch
(assert belt-speed (value (abs (* (diff (motor-turn 1500) t) (/ belt-per-turn 360))) < 4000))

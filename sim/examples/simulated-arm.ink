; Forward simulation: a link on a TowerPro MG996R servo, steered to 45 deg by a PD command
; the file writes (the servo supplies its command clamped to its 6 V line at each speed),
; and a gate on a plain hinge whose return spring its own assembly declares, let go at
; 80 deg. Neither motion is written: each joint's argument is (simulate X0 [V0]), its motion
; stepped forward from its forces (Radau IIA in the mechanism's coordinates).
; Claims, ink effort and a viewer read the simulated motions as any motion: the arm never
; passes 50 and settles by 0.3 s at 43.74 deg, 1.26 short of 45 (a PD under gravity holds
; the weight's 50 N mm with that much error: 50.4 / 40); the gate swings through closed and
; back, its spring and gravity trading its energy, never past 80. The gate hangs from its
; hinge's line alone, beside the servo: nothing joins them. Each motion is known to within
; what its bodies' measured masses move it: the gate's claim, 0.001 deg of room, is inside
; what its leaf's inertia (known to 1-2%) shifts its swing by in time, 0.4 ms, at an even
; step beside a turning point, so it is not settled.
; expect-status: gate-within unknown
(import "std")

; A TowerPro MG996R, from its A-F drawing: the case 40.9 x 20, 26.8 below its
; tabs' underside (z = 0; its mount frame there, z out of the case) and 10.2
; above; tabs 54 across, 2.4 thick, with M3 holes on 48 x 10; the output's
; axis 10.15 along x, its spline 6 across, its top (the output frame) 15.9
; up. Its record: 55 g, 159 deg of travel; at 6 V 66.6667 rpm free and
; 1.07873 N m at stall, 1.4 A. (drive ANGLE): its output is an actuator,
; the torque it can give falling from its stall torque at rest to none at
; its free speed.
(def ax 10.15)
(body case
  (shape (subtract (min (extrude-z (rounded-rect-2d 40.9 20 1.5) -26.8 10.2)
                        (extrude-z (rounded-rect-2d 54 20 1.7) 0 2.4)
                        (move (extrude-z (circle-2d 6.2) 10.2 12.4) ax 0 0))
                   (extrude-z (warp (circle-2d 1.7) (- (abs x) 24) (- (abs y) 5) z) -1 4)))
  (attr color "#2B2B2B")
  (frame mount (at (rot-x 180))))
(body spline (shape (extrude-z (circle-2d 3) 12.4 15.9)) (attr color (vec .8 .8 .82))
  (frame output (at (move 0 0 15.9))))
(def rated 66.6667rpm)                    ; deg/s
(fn torque-max (rate) (dc-torque-max 1078.73N*mm rated 1.4 1.4 (+ 0deg/s rate)))
(fn torque-min (rate) (dc-torque-min 1078.73N*mm rated 1.4 1.4 (+ 0deg/s rate)))
(assembly (servo (angle 0)) (attr mass 55g)
  (instance case case)
  (instance output spline (at (move ax 0 0)) (turn z angle) (limits -79.5 79.5) (speed-limits (neg rated) rated)
    (effort-limits (torque-min (diff angle t)) (torque-max (diff angle t))))
  (frame mount (frame-of case mount))
  (frame output (frame-of output output)))

(gravity 0 0 -9.81m/s2)
(t-span 0 1s)

; An aluminium link 120 mm long on the servo's output, a 50 g brass weight at its tip.
(body link (shape (move (rounded-box 136 16 4 4) 60 0 2)) (attr density 2.7g/cm3) (attr color (vec 0.85 0.45 0.1))
  (frame hub (at (rot-x 180))) (frame tip (at (move 120 0 4))))
(body weight (shape (move (cylinder-z 8 30) 0 0 15)) (attr mass 50g) (attr color "#B5A642"))
(assembly (arm (angle 0))
  (instance s (servo angle) (at (rot-x 90) (move 0 0 60)) mount)
  (instance link link (frame-of s output) hub)
  (instance weight weight (frame-of link tip)))

; A 60 x 80 mm aluminium leaf hanging from a hinge along y, its spring pulling it back to
; hanging, 2 N mm per degree about the hinge: a torque the gate assembly declares, so every
; model placing a gate gets it.
(body leaf (shape (move (box 60 4 80) 30 0 -40)) (attr density 2.7g/cm3) (attr color "#8FA3B8"))
(assembly (gate (angle 0))
  (instance leaf leaf (at (move 0 0 120)) (turn y angle))
  (load spring leaf (torque 0 (* 1N*mm (* -2 angle)) 0)))

(def swing (simulate 0))
(def flap (simulate 80))
(assembly root
  (instance a (arm swing))
  (instance g (gate flap) (at (move 0 200 0))))
(command root/a/s/output (- (* 40 (- 45 swing)) (* 1.2 (diff swing t))))

(require no-overshoot (value swing < 50))
(require settled (value (select (- t 0.3s) -1 (- (abs (- swing 45)) 1.5)) < 0))
(require gate-within (value (abs flap) < 80.001))

; A two-link arm in the vertical plane on two TowerPro MG996R servos, a 100 g steel block
; at its tip: is each servo enough for the move? Its output joint states the record's travel,
; 0-159 deg (+/-79.5 from centre), as limits and its rated no-load speed, 0.15 s per 60 deg
; at 6 V (400 deg/s), as speed limits: swung 120 deg the shoulder passes its travel, and every
; 0.3 s it passes its speed (a swing of A over T peaks at pi A / T deg/s). Each servo's output is an actuator
; (its effort-limits: the torque it can give at 6 V, from its stall torque at rest down to
; none at its rated speed). (gravity ...) makes their claims, effort-limits:PATH: the effort
; each supplies, found by inverse dynamics of what it moves, within its limits. The shoulder
; swings up 50 deg and back: once in 3 s it holds; every 0.4 s it passes 390 deg/s, near its
; 400 deg/s rated speed, where it can give 20 N mm, and holding the arm out takes far more.
; ink effort servo-arm-effort.ink gives each servo's peak, RMS and least margin, each joint's
; reaction (the shoulder's carries the links, the elbow servo and the block riding its output),
; and that the gate, a flap on a plain hinge moved by its argument, needs a torque nothing
; supplies (unpowered). Does the stand tip? Each arm's zero-moment point, where the floor's
; reaction acts (std's zmp), must stay on the 140 x 80 foot: the slow swing keeps it 25 mm
; inside; the fast one, starting at full acceleration, throws it past the foot's edge.
; expect-status: effort-limits:root/fast/s/output fails estimate
; expect-status: fast-stands fails estimate
; expect-status: limits:root/wide/s/output fails
; expect-status: speed-limits:root/quick/s/output fails
; expect-status: effort-limits:root/quick/s/output fails estimate

; The five setups side by side.
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
(t-span 0 3s)

(def H 90)                                   ; shoulder height
(def L1 120) (def L2 100)                    ; pivot to pivot, pivot to the block's centre
(def off 10.15)                              ; an MG996R's case centre to its output axis
(fn cutout () (box 41.3 20.4 60))            ; through which a case passes
(fn tab-holes () (warp (cylinder-z 1.7 60) (- (abs x) 24) (- (abs y) 5) z))

; An aluminium upright on a foot standing on the floor (z = 0); seat: on its -y face, z out of
; it (-y), so a servo's shaft is along -y.
(def seat (axes (origin (neg off) -3 H) (x 1 0 0) (z 0 -1 0)))
(def foot (max (- (abs x) 70) (- (abs y) 40)))   ; its footprint on the floor, in x and y
(body base
  (shape (min (extrude-z foot 0 6)
                (subtract (move (box 80 6 (+ H 30)) -10 0 (/ (+ H 30) 2))
                          (place (min (cutout) (tab-holes)) seat))))
  (attr density 2.7g/cm3) (attr color "#8FA3B8")
  (frame seat seat))

; Aluminium links, 4 thick, along +x from a hub on the output axis: hub, the hub's underside
; (z down); the upper link carries the elbow servo (its tabs on a pad, the case through a
; cutout), the forearm a block at its tip.
(fn link (l stand) (min (move (rounded-box (+ l 16) 16 4 4) (/ l 2) 0 (+ stand 2)) (move (cylinder-z 7 stand) 0 0 (/ stand 2))))
(body upper
  (shape (subtract (min (link L1 6) (move (rounded-box 62 26 4 4) (- L1 off) 0 8))
                   (move (min (cutout) (tab-holes)) (- L1 off) 0 0)))
  (attr density 2.7g/cm3) (attr color (vec 0.85 0.45 0.1))
  (frame hub (at (rot-x 180)))
  (frame seat (at (move (- L1 off) 0 10))))
(body fore
  (shape (link L2 0))
  (attr density 2.7g/cm3) (attr color (vec 0.85 0.45 0.1))
  (frame hub (at (rot-x 180)))
  (frame tip (at (move L2 0 4))))
(body block (shape (move (box 20 20 32) 0 0 16)) (attr density 7.85g/cm3) (attr color "#5A5F66"))

; The arm, its servos at `shoulder` and `elbow` degrees from the centre of their travel.
(assembly (arm (shoulder 0) (elbow 0))
  (instance base base)
  (instance s (servo shoulder) (frame-of base seat) mount)
  (instance upper upper (frame-of s output) hub)
  (instance e (servo elbow) (frame-of upper seat) mount)
  (instance fore fore (frame-of e output) hub)
  (instance block block (frame-of fore tip)))

; A gate on a plain hinge: a joint with no effort limits supplies nothing.
(body leaf (shape (move (box 60 4 80) 30 0 -40)) (attr density 2.7g/cm3))
(assembly (flap (angle 0)) (instance leaf leaf (at (move 0 0 120)) (turn y angle)))

; Up `swing` degrees and back every `period` seconds, smoothly.
(fn swing (swing period) (* 0.5deg swing (- 1 (cos (/ (* 360deg/s t) period)))))
(assembly root
  (instance slow (arm (swing 50 3) -30))
  (instance fast (arm (swing 50 0.4) -30) (at (move 0 200 0)))
  (instance wide (arm (swing 120 3) -30) (at (move 0 400 0)))
  (instance quick (arm (swing 50 0.3) -30) (at (move 0 600 0)))
  (instance gate (flap (swing 45 2)) (at (move 0 -200 0))))

; Each stand stays on its foot while the floor's reaction acts inside it: its zero-moment point
; in the floor (its own frame's z = 0), the footprint read there.
(def slow-p (zmp (reaction root/slow (frame-of root/slow))))
(def fast-p (zmp (reaction root/fast (frame-of root/fast))))
(require slow-stands (value (warp foot (nth slow-p 0) (nth slow-p 1) 0) < 0))
(require fast-stands (value (warp foot (nth fast-p 0) (nth fast-p 1) 0) < 0))

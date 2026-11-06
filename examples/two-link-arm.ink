; A two-link arm in the vertical plane, its pivots along y: two TowerPro
; SG90 servos, the shoulder on an upright, the elbow on the first link's
; end. Each frame is written in its own body's coordinates, and placement
; walks the chain: a link rides its servo's output frame, and the next
; servo stands on the link's seat. So turning the shoulder carries the
; elbow.
; A joint turns about z of the frame it is placed on. A servo's output z is
; its shaft, so a pivot along y needs a servo whose shaft points along y:
; here the upright's seat has z out of its -y face. Frames sit on the faces
; where parts touch, z out of their part: a seat on its face, a servo's
; mount on its tabs' underside (z out of the servo). Met face to face, z
; opposed, the servo stands on the seat. A rot-x 180 added to either would
; bury its tabs, and the claims shoulder-on-upright and elbow-on-link fail.
(import "std")

; A TowerPro SG90, from its A-F drawing: the case 22.7 x 12.2, 17 below its
; tabs' underside (z = 0; its mount frame there, z out of the case) and 10
; above; tabs 32.3 across, 1.8 thick, with 2.2 holes 27 apart; the output's
; axis 5.3 along x, its spline 4.8 across, its top (the output frame) 13.3
; up. Its record: 9 g, 150 deg of travel; at 4.8 V 83.3333 rpm free.
; (servo ANGLE) turns its output ANGLE degrees from the centre of its travel,
; within its free speed either way.
(def ax 5.3)
(body case
  (shape (subtract (min (extrude-z (rounded-rect-2d 22.7 12.2 1) -17 10)
                        (extrude-z (rounded-rect-2d 32.3 12.2 1.1) 0 1.8)
                        (move (extrude-z (circle-2d 4.45) 10 11.3) ax 0 0))
                   (extrude-z (warp (circle-2d 1.1) (- (abs x) 13.5) y z) -1 3)))
  (attr color "#2F6FB5")
  (frame mount (at (rot-x 180))))
(body spline (shape (extrude-z (circle-2d 2.4) 11.3 13.3)) (attr color "#F2F2F2")
  (frame output (at (move 0 0 13.3))))
(def rated 83.3333rpm)                    ; deg/s
(assembly (servo (angle 0)) (attr mass 9)
  (instance case case)
  (instance output spline (at (move ax 0 0)) (turn z angle) (limits -75 75) (speed-limits (neg rated) rated))
  (frame mount (frame-of case mount))
  (frame output (frame-of output output)))

(param shoulder 30)                          ; deg up from +x, each servo
(param elbow 45)                             ; from the centre of its travel
(def H 70)                                   ; shoulder height
(def L1 70) (def L2 60)                      ; pivot to pivot, pivot to tip
(def off 5.3)                                ; an SG90's case centre to its output axis
(fn cutout () (box 23.3 12.8 40))            ; through which a case passes
(fn tab-holes () (min (move (cylinder-z 1.1 40) 13.5 0 0) (move (cylinder-z 1.1 40) -13.5 0 0)))

; An upright on a foot; seat: on its -y face, z out of it (-y), x along x,
; the case centre 5.3 behind the shoulder axis at x = 0.
(def seat (axes (origin (neg off) -2 H) (x 1 0 0) (z 0 -1 0)))
(body base
  (shape (min (move (box 90 60 4) 0 0 2)
                (subtract (move (box 50 4 (+ H 20)) -5 0 (/ (+ H 20) 2))
                          (place (min (cutout) (tab-holes)) seat))))
  (attr color "#8FA3B8")
  (frame seat seat))

; A link along its +x, 3 thick on a boss `stand` tall (which holds the
; next servo's case clear of the upright), its hub at the origin: hub, its
; underside, z down; with a seat (the next servo's tabs on its top, the
; case through a cutout) or a tip at the far end.
(fn link (l stand) (min (move (rounded-box (+ l 10) 10 3 3) (/ l 2) 0 (+ stand 1.5)) (move (cylinder-z 5 stand) 0 0 (/ stand 2))))
(body upper
  (shape (subtract (min (link L1 5) (move (rounded-box 36 20 3 3) (- L1 off) 0 6.5))
                   (move (min (cutout) (tab-holes)) (- L1 off) 0 0)))
  (attr color (vec 0.85 0.45 0.1))
  (frame hub (at (rot-x 180)))
  (frame seat (at (move (- L1 off) 0 8))))
(body fore
  (shape (link L2 0))
  (attr color (vec 0.85 0.45 0.1))
  (frame hub (at (rot-x 180)))
  (frame tip (at (move L2 0 1.5))))

(assembly root
  (instance base base)
  (instance s (servo shoulder) (frame-of base seat) mount)
  (instance upper upper (frame-of s output) hub)
  (instance e (servo elbow) (frame-of upper seat) mount)
  (instance fore fore (frame-of e output) hub))

; The tip is where the angles put it; each servo stands on what holds it;
; the elbow servo clears the upright.
(def tip (get (frame-of root/fore tip) "origin"))
(def a1 (deg->rad shoulder))
(def a2 (deg->rad (+ shoulder elbow)))
(assert tip-x (value (abs (- (nth tip 0) (+ (* L1 (cos a1)) (* L2 (cos a2))))) < 0.001))
(assert tip-z (value (abs (- (nth tip 2) (+ H (* L1 (sin a1)) (* L2 (sin a2))))) < 0.001))
(assert shoulder-on-upright (touches (shape-of root/s) (shape-of root/base)))
(assert elbow-on-link (touches (shape-of root/e) (shape-of root/upper)))
(assert elbow-clear (clearance (shape-of root/e) (shape-of root/base) 2))

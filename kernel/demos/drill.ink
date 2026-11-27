; A cordless drill-driver, whole: a clamshell housing (motor barrel, raked
; pistol grip with finger grooves, battery foot) joined by large smooth
; blends into one organic shell, with its parting seam and screw bosses; a
; 1.6 mm soft overmold over the grip that is the housing grown by 1.6 mm
; there, its face textured by a gyroid; a 20-step clutch collar, a keyless
; chuck with its three jaws on a 1/4 in hex bit, a trigger in its recess,
; the forward/reverse button, a work light, and a slide-on battery with its
; release buttons and charge gauge. Every blend, offset and texture is a
; field: nothing here can fail to fillet or offset.
(import "std")

(param rake 15)        ; grip angle from vertical, deg
(param overmold 1.6)   ; soft layer, mm
(param groove 3.0)     ; finger groove depth, mm

; ------------------------------------------------------------ housing
; The motor barrel along x: 27 mm radius at the back to 23 at the gearbox,
; its back end rounded.
(def barrel (smooth-intersect (rotate-y (cone-z 27 23 -75 62) 90) (slab-x -75 62) 10))

; The grip: an ellipse section 40 mm fore and aft by 31 across, on an axis
; raked back from (10, 0, -20) downward; its front face carries four
; finger grooves 24 mm apart. (s, w) are the point's coordinates along the
; axis (down) and fore-aft (forward) in the grip's plane.
(def c (cos (deg->rad rake)))
(def k (sin (deg->rad rake)))
(def s (+ (* (- 10 x) k) (* (- -20 z) c)))
(def w (- (* (- x 10) c) (* (- z -20) k)))
(def bulge (* 2 (sin (* pi (clamp (/ (- s 10) 100) 0 1)))))
(def fingers (* groove (clamp (/ w 10) 0 1)
                (max 0 (sin (* 2 pi (/ (- s 30) 24))))
                (- 1 (clamp (/ (abs (- s 66)) 40) 0 1))))
(def half-a (- (+ 20 bulge) fingers))
(def half-b (+ 15.5 (* 0.5 bulge)))
(def ellipse (* (- (len2 (/ w half-a) (/ y half-b)) 1) 15))
(def grip (max ellipse (- (abs (- s 60)) 66)))

; The foot the battery slides onto, its front carrying the work light.
(def foot (move (rounded-box 100 70 16 6) -6 0 -143))
(def shell0 (smooth-union (smooth-union barrel grip 16) foot 12))

; The trigger, in a recess 1 mm round it; the forward/reverse button
; through the housing above it.
(def trigger-f (move (rotate-y (smooth-union (rounded-box 11 13 30 4.5) (move (rounded-box 6 9 8 2.5) -6 0 12) 2) (neg rake))
                     26 0 -50))
(def switch-f (move (rotate-x (capsule-z 46 4.2) 90) 13 23 -32))

; The clamshell's parting seam (0.7 mm wide, 0.6 deep, on y = 0), motor
; vents (20 slots round the back), and four screw bosses on the left half.
(fn seam-in (f) (max (- (abs y) 0.35) (neg (+ f 0.6))))
(def seam (seam-in shell0))
(def vents (move (warp (repeat-polar (move (rounded-box 8 3.2 18 1.2) 27 0 0) 20) y z x) -50 0 0))
(fn boss (bx bz) (move (cylinder-y 2.8 30) bx -30 bz))
(def bosses (min (boss -55 -14) (boss 35 -14) (boss -12 15) (boss -40 -140)))
(def housing-f (subtract shell0 (- trigger-f 1) seam vents bosses))
; The bosses' screws, heads 2 mm down.
(def screws (max bosses (+ (min (rotate-y (cone-z 27 23 -75 62) 90) foot) 2)))

; The overmold: the housing grown by overmold over the grip, between its
; blends (s from 24 to 106, edges rounded by a 0.5 mm blend: a 2 mm one
; would swallow the 1.6 mm layer's ends whole), its face
; textured by a gyroid 5 mm in period, 0.4 mm proud or sunk; and a band
; round the back of the barrel.
(def zone (min (max (- (abs (- s 65)) 41) (- ellipse 12))
               (max (- (abs (+ x 70)) 6) (- (len2 y z) 40))))
(def texture (* 0.6 (gyroid 5)))
(def upper (smooth-union barrel grip 16))   ; the shell where the overmold lies (the foot's blend is below it)
(def soft (subtract (smooth-intersect (max (- upper (+ overmold texture)) (neg upper)) zone 0.5)
                    (- trigger-f 1) (seam-in upper)))

; ------------------------------------------------------------ front end
; The clutch collar: 20 grip ribs; then the keyless chuck, a tapered
; sleeve with 24 flutes, its nose and three jaws closed on the bit.
(fn about-x (f n) (warp (repeat-polar (warp f z x y) n) y z x))
(def collar (subtract (fillet-intersect (cylinder-x 25 300) (slab-x 62 74) 1.5 90)
                      (about-x (move (cylinder-x 1.6 300) 0 25.4 0) 20)))
(def sleeve (subtract (smooth-intersect (rotate-y (cone-z 21 18 76 112) 90) (slab-x 76 112) 3)
                      (about-x (move (rotate-z (cylinder-x 1.6 60) -4.8) 94 20.3 0) 24)))
(def nose (subtract (fillet-intersect (rotate-y (cone-z 15 10 112 124) 90) (slab-x 112 124) 2 90)
                    (cylinder-x 6.5 400)))
(def jaw (max (- (abs z) 1.6) (- 3.4 y) (- y (- 6.6 (* 0.25 (- x 118)))) (slab-x 116 127)))
(def jaws (max (about-x jaw 3) (slab-x 115 128) (- (len2 y z) 8)))
(def chuck-f (min sleeve nose))

; A 1/4 in (6.35 mm) hex bit, pointed.
(def bit-f (min (extrude-x (hexagon-2d 6.35) 116 158)
                (fillet-intersect (rotate-y (cone-z 3.4 0.8 158 170) 90) (slab-x 158 170) 0.5 90)))

; ------------------------------------------------------------ battery
; A slide-on pack under the foot: rounded, a release button each side,
; a black base with ribs, and a three-step charge gauge at the back.
(def pack (move (rounded-box 104 78 36 7) -6 0 -170))
(def ribs (max (move (repeat-x (box 2 90 40) 7) 0 0 -170) (slab-z -188 -182) (slab-x -50 40)))
(def release (symmetric-y (move (rounded-box 26 6 10 2.5) -16 39 -162)))
(def pack-f (subtract (min pack release) (max (- (abs (+ z 158)) 0.4) (neg (+ pack 1))) ribs))
(def gauge (min (move (rounded-box 3 5 2 0.8) -58.5 -8 -160)
                (move (rounded-box 3 5 2 0.8) -58.5 0 -160)
                (move (rounded-box 3 5 2 0.8) -58.5 8 -160)))
(def light (move (rotate-y (cylinder-z 5 2) 90) 44.2 0 -143))

; ------------------------------------------------------------ bodies
(def graphite (vec 0.20 0.21 0.23))
(def accent (vec 0.05 0.45 0.62))
(body housing
  (shape (max housing-f (move (box 240 100 220) 15 0 -60)))
  (attr material "abs")
  ; graphite, with an accent stripe round the gearbox and on the foot's lip
  (attr color (select (min (max (- 52 x) (- x 61)) (- (abs (+ z 139)) 1.2)) accent graphite))
  (attr roughness 0.38))
(body overmold (shape soft) (attr material "tpu") (attr color (vec 0.07 0.07 0.075)) (attr roughness 0.9))
(body collar (shape collar) (attr material "abs") (attr color (vec 0.08 0.08 0.085)) (attr roughness 0.6))
(body chuck (shape chuck-f) (attr material "abs") (attr color (vec 0.06 0.06 0.065)) (attr roughness 0.5))
(body jaws (shape jaws) (attr material "steel") (attr color (vec 0.75 0.76 0.78)) (attr metallic 1) (attr roughness 0.25))
(body bit (shape bit-f) (attr material "steel") (attr color (vec 0.62 0.63 0.66)) (attr metallic 1) (attr roughness 0.3))
(body trigger (shape trigger-f) (attr material "pom") (attr color accent) (attr roughness 0.45))
(body switch (shape switch-f) (attr material "pom") (attr color (vec 0.10 0.10 0.11)) (attr roughness 0.5))
(body screws (shape screws) (attr material "steel") (attr color (vec 0.35 0.36 0.38)) (attr metallic 1) (attr roughness 0.4))
(body battery
  (shape (max pack-f (move (box 120 100 60) -6 0 -170)))
  (attr material "abs")
  (attr color (select (- z -180) (vec 0.06 0.06 0.065) graphite))
  (attr roughness 0.45))
(body lights
  (shape (min light gauge))
  (attr color (select (- x 0) (vec 0.25 0.95 0.35) (vec 0.97 0.97 0.92)))
  (attr roughness 0.15))

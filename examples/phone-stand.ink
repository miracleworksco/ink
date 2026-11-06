; A phone stand whose plate leans back at 65 degrees to the desk: a part on
; a tilted face of your own body (ink spec --section "A part on a tilted
; face"). The plate's front face is a frame written as ink prints one: an
; origin, an x along the desk, a z out of the face. The plate is built flat
; in that frame and put there. The phone's back meets the face 3 mm above
; the lip.
(import "std")

; An iPhone 16 Pro, from Apple's dimensional drawing: 71.45 x 149.61 x 8.25,
; its back at z = 0 (back: there, z out of it), screen up, top +y; its
; camera plateau 2.05 proud of the back from 28.26 to 74 above its centre,
; three lenses 16.2 across standing 4.28 proud.
(def phone-height 149.61)
(def lens (extrude-z (circle-2d 8.1) -4.28 0.05))
(body phone
  (shape (min (extrude-z (rounded-rect-2d 71.45 phone-height 12) 0 8.25)
              (move (extrude-z (rounded-rect-2d 44.42 45.74 9) -2.05 0.05) 12.71 51.13 0)
              (move lens 21.555 60.635 0) (move lens 21.555 41.395 0) (move lens 2.565 51.015 0)))
  (attr mass 199) (attr color (vec .37 .39 .42)) (attr metallic .35) (attr roughness .35)
  (frame back (at (rot-x 180)))
  (frame front (at (move 0 0 8.25))))

(def tilt 65)                                      ; the plate, from the desk (degrees)
(def a (deg->rad tilt))
(def face (axes (origin 0 0 8) (x 1 0 0) (z 0 (neg (sin a)) (cos a))))  ; its y runs up the plate
(def plate (min (move (box 80 100 5) 0 50 -2.5)    ; in the face's frame, the face at z = 0
                (move (box 80 6 12) 0 -3 1)))       ; a lip below the phone
(body stand (attr color "#8FA3B8")
  (shape (min (move (box 90 110 6) 0 15 3)        ; the base, on the desk at z = 0
              (place plate face)))
  ; The phone's back (its centre) 3 mm above the lip, up the face's own y.
  (frame seat (at (move 0 (+ (/ phone-height 2) 3) 0) face)))

(assembly root
  (instance stand stand)
  (instance phone phone (frame-of stand seat) back))

; The phone's back lies on the face, and its camera bump, 28 mm above its
; centre, clears the plate's top edge: a longer plate fails here.
(assert apart (disjoint (shape-of root/stand) (shape-of root/phone)))

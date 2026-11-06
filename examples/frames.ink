; Places parts by named frames instead of coordinates. A headed pin sits in
; a plate's bore: the plate's frame is on its top face at the bore, z out
; of the plate; the pin's is under its head, z out of the pin. The instance
; meets them face to face, so moving the hole moves the pin.
;
; A 608 bearing rests in a counterbore, its side face on the shoulder. It is
; modelled about its centre, as a catalog part is, so placing it by numbers,
; (at (move 0 0 4.5)), sinks it 3.5 into the shoulder. By frames it lands on
; the face: the housing's seat (on the shoulder, z up) meets the bearing's
; face-b (its -z side face, z out of the bearing), face to face.
(import "std")

(param hole-x 12)
(param hole-y 6)

(def slab (box 60 40 8))
(def bore (move (cylinder-z 5 20) hole-x hole-y 0))
(body plate (shape (max slab (- 0 bore))) (attr color "#8FA3B8")
  (frame bore-top (at (move hole-x hole-y 4))))

(def shaft (cylinder-z 4.8 24))
(def head (move (cylinder-z 8 3) 0 0 13.5))
(body pin (shape (min shaft head)) (attr color "#D79B00")
  (frame seat (at (rot-x 180) (move 0 0 12))))   ; under the head, z down

; A 608 (8 bore, 22 outside, 7 wide): its two rings, the gap between them
; (radii 6 to 9.5) recessed 0.3 from its side faces.
(body bearing
  (shape (subtract (cylinder-z 11 7) (cylinder-z 4 8) (max (ring-2d 6 9.5) (- 3.2 (abs z)))))
  (attr color (vec .78 .78 .8)) (attr metallic .9)
  (frame face-a (at (move 0 0 3.5)))                  ; +z side face, z out
  (frame face-b (at (rot-x 180) (move 0 0 -3.5))))    ; -z side face

(def ledge 4.5)                              ; the shoulder's height
(def at-y 60)                                ; the housing beside the plate
(body housing
  (shape (subtract (move (rounded-box 40 40 12 3) 0 0 6)
                   (move (cylinder-z 11.02 20) 0 0 (+ ledge 10))  ; the bore, 0.02 clear
                   (cylinder-z 8 40)))                          ; under the outer ring only
  (attr color "#8FA3B8")
  (frame seat (at (move 0 0 ledge))))

(param spin 0)
(assembly root
  (instance base plate)
  (instance peg pin (frame-of base bore-top) seat (turn z spin))
  (instance housing housing (at (move 0 at-y 0)))
  (instance bearing bearing (frame-of housing seat) face-b))

; Seated: the head is material at a point just above the plate's top (z = 4)
; beside the bore, so a pin raised off the plate fails.
; The shaft clears the bore by 0.2 all round.
(assert seated (subset (probe (+ hole-x 6.5) hole-y 4.02) (shape-of root/peg)))
(assert clear (clearance (shape-of root/base) (shape-of root/peg shaft) 0.15))

; The bearing rests on the shoulder: it touches the housing below the seat,
; neither deeper in nor farther off than the accuracy. Only the housing
; below the seat counts: the bore is 0.02 clear, closer than the accuracy,
; so a bearing afloat against the bore's wall would touch the whole
; housing. And it sits down in the bore, below the housing's top.
(assert rests (touches (shape-of root/bearing) (max (shape-of root/housing) (- z ledge))))
(assert in-bore (within (shape-of root/bearing) (-11.1 (- at-y 11.1) ledge) (11.1 (+ at-y 11.1) 12)))

; A sliding drawer: a box on a pair of Accuride 3832 slides in a cabinet.
; Each slide's cabinet member meets, by its mount frame, a seat on a wall's
; inner face (z out of the wall, into the cabinet). The box meets the right
; slide's drawer-mount, which (drawer-slide EXTENSION) moves with the
; drawer member, so the box follows the slide. The left slide is the same
; part met to the left wall's seat, whose z points the other way: the
; frames turn it over, so nothing is mirrored by hand.
(import "std")

; An Accuride 3832-250 (its datasheet: 250 long, 45.7 tall, 12.7 side-space,
; 254 travel), three steel channels: the cabinet member, its web on the
; wall (y = 0 to 1.2) and its flanges out to y = 6; the intermediate inside
; it; the drawer member, its face at y = 12.7. mount: the cabinet member's
; face on the wall, z out of the slide (-y); drawer-mount: the drawer
; member's face, z out of it (+y). (drawer-slide EXTENSION): the drawer
; member out EXTENSION along -x, the intermediate half as far.
(fn channel (w0 w1 f0 f1 h)                  ; its web from y w0 to w1, flanges f0 to f1, h tall
  (max (min (slab-y w0 w1) (max (slab-y f0 f1) (- (- (* 0.5 h) 1.2) (abs z)))) (- (abs z) (* 0.5 h)) (slab-x 0 250)))
(def steel (vec .64 .67 .7))
(body slide-cabinet (shape (channel 0 1.2 0 6 45.7)) (attr color steel) (attr metallic .65)
  (frame mount (at (rot-x 90))))
(body slide-middle (shape (channel 6.8 8 2 8 42.5)) (attr color steel) (attr metallic .65))
(body slide-drawer (shape (channel 11.5 12.7 8.5 12.7 36)) (attr color steel) (attr metallic .65)
  (frame mount (at (rot-x -90) (move 0 12.7 0))))
(assembly (drawer-slide (extension 0)) (attr mass 760)
  (instance cabinet slide-cabinet)
  (instance intermediate slide-middle (frame-of cabinet) (slide x (* -0.5 extension)))
  (instance drawer slide-drawer (frame-of cabinet) (slide x (neg extension)) (limits -254 0))
  (frame mount (frame-of cabinet mount))
  (frame drawer-mount (frame-of drawer mount)))

(param pull 200)                             ; how far the drawer is out, mm
(def w 200) (def deep 250) (def high 80)     ; the box, outside
(def gap 12.7)                               ; slide side-space, wall to box
(def mid 40)                                 ; the slides' centre height
(def wall-y (+ (/ w 2) gap))
(body cabinet
  (shape (min (move (box 280 15 100) 130 (neg (+ wall-y 7.5)) 40)     ; right wall
                (move (box 280 15 100) 130 (+ wall-y 7.5) 40)           ; left wall
                (move (box 280 (* 2 (+ wall-y 15)) 15) 130 0 -12.5)     ; floor, 5 below the box
                (move (box 15 (* 2 (+ wall-y 15)) 100) 262.5 0 40)))    ; back
  (attr color "#B79162")
  (frame seat-r (at (rot-x -90) (move 0 (neg wall-y) mid)))  ; z +y, out of the right wall
  (frame seat-l (at (rot-x 90) (move 0 wall-y mid))))        ; z -y, out of the left wall
(body drawer-box
  (shape (subtract (move (box deep w high) (/ deep 2) 0 (/ high 2))
                   (move (box (- deep 16) (- w 16) high) (/ deep 2) 0 (+ (/ high 2) 8))))
  (attr color "#D8C3A0")
  (frame side-r (at (rot-x 90) (move 0 (/ w -2) mid))))      ; right outer face, z -y

(assembly root
  (instance cabinet cabinet)
  (instance right (drawer-slide pull) (frame-of cabinet seat-r) mount)
  (instance left (drawer-slide pull) (frame-of cabinet seat-l) mount)
  (instance box drawer-box (frame-of right drawer-mount) side-r))

; The box rides both drawer members, each slide sits on its wall, the box
; clears the cabinet, and its front is `pull` out (the slides' travel is
; their limits claim). `--param pull=0` checks it shut.
(assert on-right (touches (shape-of root/box) (shape-of root/right/drawer)))
(assert on-left (touches (shape-of root/box) (shape-of root/left/drawer)))
(assert right-on-wall (touches (shape-of root/right/cabinet) (shape-of root/cabinet)))
(assert left-on-wall (touches (shape-of root/left/cabinet) (shape-of root/cabinet)))
(assert box-clear (clearance (shape-of root/box) (shape-of root/cabinet) 4))
(assert pulled-out (subset (probe (- 4 pull) 0 4) (shape-of root/box)))

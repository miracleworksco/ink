; An injection-moulded enclosure half, checked for undercuts along its
; mould's pull (docs/IR.md, Manufacturing checks). Open at the bottom, it
; leaves a two-part mould along z: the cavity pulls up, the core down. A USB
; opening in a side wall is undercut: no line of sight from its air runs
; out along +z or -z. A side action pulling along +x reaches it, so the
; part moulds with one: every point of the air sees out along +z, -z or +x.
;
; The band checked is the mould's air 1 mm clear of the part. So an
; undercut shallower than about 1.1 mm passes here, such as a snap bead:
; 1 mm plus twice the 0.05 accuracy (how far off an answer may be). Less
; clearance sees shallower, at more cost.
; Walls 2 mm. Draft is not modelled; a draft check reads a slope (Claims,
; Slopes and normals).
(import "std")

(def up (vec 0 0 1))
(def down (vec 0 0 -1))
(def side (vec 1 0 0))

(def shell (subtract (extrude-z (rounded-rect-2d 50 30 4) 0 20) (extrude-z (rounded-rect-2d 46 26 2) -1 18)))
(def usb (move (box 10 8 4) 24 0 8))                    ; through the +x wall
(def plain shell)
(def mould (move (box 54 34 24) 0 0 10))                ; the mould's box, 2 mm past the part
(def ported (subtract shell usb))
(body enclosure (shape plain) (attr material "abs"))
(body enclosure-usb (shape ported) (attr material "abs"))

(assert pulls (empty (hidden plain mould (clear-air plain mould 1) (vec up down))))
; expect-status: usb-pulls fails
(assert usb-pulls (empty (hidden ported mould (clear-air ported mould 1) (vec up down))))
(assert pulls-with-side-action (empty (hidden ported mould (clear-air ported mould 1) (vec up down side))))

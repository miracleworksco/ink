; expect-warning: body 'barrel': a fillet-union joins two shapes whose proven bounds both end near x = 150.146
;   (a false alarm: x = 150 is the skin's outer radius, a curve, and the
;   ribs end at r 149, inside it; their proven bounds only reach it.)
; An isogrid barrel section, the stiffened shell of rocket stages and
; satellite buses: a 300 mm diameter, 240 mm tall aluminium cylinder,
; machined from the inside into a 1.5 mm skin and a triangular grid of
; ribs 8 mm deep, every rib and node filleted into the skin, between two
; end rings drilled for 48 bolts each. The grid is three families of
; parallel ribs, each one field repeated round the cylinder in unrolled
; coordinates (arc, height), so its hundreds of pockets cost one, and the
; fillets at every rib, node and ring come from two blends. A quarter is
; cut away to show the inside, primer-green as such parts are.
(import "std")

(param radius 150)     ; outer radius, mm
(param height 240)     ; mm
(param skin 1.5)       ; mm
(param rib-h 8)        ; rib depth inward from the skin, mm
(param rib-w 2.2)      ; rib width, mm
(param around 30)      ; inclined ribs round the circumference
(param ring-h 12)      ; end ring height, mm

(def half (* 0.5 height))
(def r (len2 x y))
(def c30 (* 0.5 (sqrt 3)))
; The grid spacing that closes round the cylinder: the inclined families
; repeat every h / cos 30 along the arc.
(def h (/ (* 2 pi radius c30) around))

; Unrolled coordinates at the skin: a along the circumference, z up.
(def a (* radius (atan2 y x)))
; The grid's nodes sit at z = z0 + 2 k h (a ring rib at mid-height).
(def z0 0)
(fn family (nx nz shift)
  (let ((p (+ (* nx a) (* nz (- z z0)) shift)))
    (- (abs (- (mod (+ p (* 0.5 h)) h) (* 0.5 h))) (* 0.5 rib-w))))
(def rings-family (family 0 1 0))
(def grid (smooth-union (smooth-union (family c30 0.5 0) (family (neg c30) 0.5 0) 3)
                        rings-family 3))

(def inner (- radius skin))
(def shell (max (- r radius) (- inner r) (- (abs z) half)))
(def ribs (max grid (- r (+ inner 0.5)) (- (- inner rib-h) r) (- (abs z) half)))
(def rings (max (- r radius) (- (- inner 10) r) (- (abs z) half) (- (- half ring-h) (abs z))))
; 6.6 mm (M6 clearance) holes through each ring, 48 round it.
(def bolt (move (cylinder-z 3.3 (* 2.2 ring-h)) (- radius 6.5) 0 half))
(def bolts (symmetric-z (repeat-polar bolt 48)))
(def barrel (subtract (fillet-union (min shell rings) ribs 2.5 90) bolts))

; A quarter cut away (x > 0, y < 0).
(def cut (min x (neg y)))

(body barrel
  (shape (max barrel cut (cylinder-z (+ radius 1) (+ height 2))))
  (attr material "aluminum-7075")
  (attr color (mix (vec 0.80 0.81 0.83) (vec 0.33 0.38 0.12) (clamp (+ 0.5 (/ (- (+ inner 0.3) r) 0.4)) 0 1)))
  (attr metallic 0.5)
  (attr roughness 0.4))

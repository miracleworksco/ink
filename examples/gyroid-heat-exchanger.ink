; A gyroid heat exchanger built as one field, with no faces, edges or
; surface patches to keep consistent. One gyroid wall fills a barrel and
; splits it into two interlocking channels that never touch, so a hot and a
; cold fluid meet across a wall everywhere at once. Each face of the wall is
; coloured by the gyroid field's sign: warm one side, cool the other. A
; quarter of the barrel is cut away, and the lattice is filleted into its
; wall.
(import "std")

(param period 14)      ; gyroid cell, mm
(param wall 1.2)       ; gyroid wall, mm
(param skin 2)         ; barrel wall, mm
(param length 120)     ; barrel length, mm
(param twist 180)      ; the lattice's turn from one end to the other, deg

; The barrel: a revolved profile 30 mm in radius at the ends and 36 at the
; middle.
(def half (* 0.5 length))
(def radius (+ 30 (* 6 (cos (/ (* pi y) length)))))
(def ends (- (abs y) half))
(def barrel (revolve-z (max (- x radius) ends)))
(def bore (revolve-z (max (- x (- radius skin)) ends)))  ; open at both ends

; The gyroid wall, inside the barrel, filleted into its skin, its cells
; turning about the barrel's axis along its length. The warp accounts for
; the twist when bounding the gyroid's slope.
(def turn (* (deg->rad twist) (/ z length)))
(def g (warp (gyroid period)
             (+ (* x (cos turn)) (* y (sin turn)))
             (- (* y (cos turn)) (* x (sin turn)))
             z))
(def core (lattice-fill bore g wall))
(def skin-wall (revolve-z (max (- (abs (- x (- radius (* 0.5 skin)))) (* 0.5 skin)) ends)))
(def body-field (fillet-union skin-wall core 1.5 90))

; A flange at each end, with eight bolt holes.
(def flange-ring (max (- (len2 x y) 44) (neg (- (len2 x y) 28)) (- (abs z) 4)))
(def bolts (pattern-polar (move (cylinder-z 3 20) 38 0 0) 8))
(def flange (max flange-ring (neg bolts)))
(def flanges (min (move flange 0 0 (- half 4)) (move flange 0 0 (- 4 half))))

; Where the wall is painted: inside the bore, short of the flanges' end
; faces by 1 mm (a flange's face at the bore's end would take the paint).
(def paint-bore (revolve-z (max (- x (- radius skin)) (+ ends 1))))

; A quarter cut away (x > 0 and y < 0) to show the inside.
(def cut (min x (neg y)))

(body exchanger
  (shape (max (min body-field flanges) cut))
  ; The wall warm on the face toward one channel and cool toward the other
  ; (the gyroid's sign); the barrel and flanges bare metal.
  (attr color (mix (vec 0.72 0.74 0.77)
                   (mix (vec 0.15 0.55 0.95) (vec 1.0 0.45 0.12) (clamp (+ 0.5 (/ g wall)) 0 1))
                   (clamp (- 0.5 paint-bore) 0 1)))
  (attr metallic 0.6)
  (attr roughness 0.35))

; The two channels, the bore on either side of the wall, never meet: the
; wall keeps them at least 0.5 mm apart everywhere (the twist thins its 1.2
; to about 0.9).
(def hot (max bore (- (* 0.5 wall) g)))
(def cold (max bore (+ g (* 0.5 wall))))
(assert channels-apart (clearance hot cold 0.5))

; A pin-fin heatsink whose pin density follows the heat: 80 x 80 mm, under
; an 80 mm fan, for a 15 x 15 mm processor die at its centre. The pins are
; packed tightest over the die, where the heat flux is highest, and spread
; out toward the edges, where the air enters: 829 tapered pins, 1.1 mm
; across over the die and up to 2.3 mm toward the corners, filleted to the
; base. They are one field: one pin, written once, repeated over a grid
; stretched by a sinh map, so no pin is placed, booleaned or filleted one
; by one.
(import "std")

(param side 80)        ; base, mm square
(param base-t 4)       ; base thickness, mm
(param pin-h 20)       ; pin height above the base, mm
(param pitch0 1.8)     ; pin pitch over the die's centre, mm
(param grade 0.06)     ; how fast the pitch grows away from the centre, 1/mm
(param fill 0.30)      ; pin radius at the root, as a fraction of the local pitch
(param root-r 0.6)     ; fillet where each pin meets the base, mm

(def half (* 0.5 side))

; The grid: index space u is uniform; the plate's coordinate is
; X(u) = sinh(a p0 u) / a, so the pitch dX/du = p0 cosh(a p0 u) is p0 at the
; centre and grows outward. Its inverse is asinh(a x) / (a p0).
(fn asinh (v) (ln (+ v (sqrt (+ 1 (square v))))))
(fn sinh (v) (* 0.5 (- (exp v) (exp (neg v)))))
(fn to-u (v) (/ (asinh (* grade v)) (* grade pitch0)))
(fn to-x (u) (/ (sinh (* grade pitch0 u)) grade))

; The pin pitch here, continuous: dX/du at u(x), p0 sqrt(1 + (a x)^2) along
; each axis, the closer of the two.
(def pitch (* pitch0 (sqrt (+ 1 (square (* grade (min (abs x) (abs y))))))))
(def r0 (* fill pitch))

; One pin at (cx, cy): a cone from r0 at the base to 60 % of it at the tip,
; capped by a sphere. None where the centre is within 4.5 mm of a corner
; hole (the screwdriver's way) or 2.5 mm of the edge.
(def hole-at (- half 4.5))
(fn pin-at (cx cy)
  (let ((rr (len2 (- x cx) (- y cy)))
        (away (max (- 4.5 (len2 (- (abs cx) hole-at) (- (abs cy) hole-at)))
                   (- (max (abs cx) (abs cy)) (- half 2.5))))
        (cone (max (- rr (mix r0 (* 0.6 r0) (clamp (/ z pin-h) 0 1))) (- z pin-h)))
        (tip (- (len3 (- x cx) (- y cy) (- z pin-h)) (* 0.6 r0))))
    (select away (min cone tip) (max (min cone tip) 10))))

; The nearest pin: of the grid points either side of the point in u, along
; each axis, whichever is nearest in the plate (so the field is continuous
; across the cells' edges, which a stretched grid puts off the midpoints).
(def iu (floor (to-u x)))
(def iv (floor (to-u y)))
(def us (vec (to-x iu) (to-x (+ iu 1))))
(def vs (vec (to-x iv) (to-x (+ iv 1))))
(def pins
  (max (min (pin-at (nth us 0) (nth vs 0)) (pin-at (nth us 1) (nth vs 0))
            (pin-at (nth us 0) (nth vs 1)) (pin-at (nth us 1) (nth vs 1)))
       (- -1 z) (- (max (abs x) (abs y)) half)))

; The base, with an M3 clearance hole near each corner for the clip.
(def plate (move (rounded-box side side base-t 2) 0 0 (* -0.5 base-t)))
(def holes (pattern-polar (move (cylinder-z 1.7 100) hole-at hole-at 0) 4))

; The pins end at the plate's sides, where the plate's edge is rounded: the
; root fillet is asked within that square only (its cap).
(def sink (subtract (fillet-union plate pins root-r 90 (- (max (abs x) (abs y)) half)) holes))

; The heat flux from the die: a Gaussian 15 mm across, peak 1 at the centre.
(def heat (exp (neg (/ (+ (square x) (square y)) (* 2 (square 9))))))

(body heatsink
  (shape (max sink (box (+ side 2) (+ side 2) 80)))
  (attr material "aluminum-6063")
  ; Bare aluminium, glowing where the die heats it (warmest at the base).
  (attr color (mix (vec 0.70 0.72 0.75) (vec 1.0 0.32 0.05)
                   (clamp (* 1.2 heat (- 1 (* 0.5 (clamp (/ z pin-h) 0 1)))) 0 1)))
  (attr metallic 0.7)
  (attr roughness 0.35)
  (frame die (at (rot-x 180) (move 0 0 (neg base-t)))))

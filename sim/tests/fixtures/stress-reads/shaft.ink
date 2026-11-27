; Benchmark 28's stepped shaft: its stress at the left shoulder's r3 fillet,
; CalculiX (C3D10, 3 meshes) 151.9, Peterson ~150. The keyway's floor meets
; its walls in sharp inside corners, where the elastic stress is unbounded:
; over the whole part, `shoulder-peak-over-155` is true and
; `shoulder-peak-band` (< 1e6) false.
(import "std")

(param fillet-r 3)           ; shoulder fillet radius
(param hub-bore-r 15.01)     ; pulley hub bore (shaft 15.000: 0.01 radial clearance)

(def seat-l (move (cylinder-z 10 50) 0 0 5))          ; z -20 .. 30
(def body-mid (move (cylinder-z 15 160) 0 0 100))     ; z 20 .. 180
(def seat-r (move (cylinder-z 10 30) 0 0 185))        ; z 170 .. 200
(def filleted (fillet-union (fillet-union seat-l body-mid fillet-r 90) seat-r fillet-r 90))
(def keyway (max (- (len2 x (- z (clamp z 89 111))) 4) (- 11 y)))
(def shaft-field (subtract filleted keyway))
(body shaft (shape shaft-field) (attr material "steel") (attr density 7.85g/cm3) (attr color "#8A8F96"))

(def hub-slot (max (- (abs x) 4.02) (neg y) (- y 18.3)))
(def hub-bore (- (len2 x y) hub-bore-r))
(def hub-field (subtract (move (cylinder-z 35 24) 0 0 100) (min hub-bore hub-slot)))
(body hub (shape hub-field) (attr material "aluminum-6061") (attr color "#C8A06A"))
(def key-field (max (- (len2 x (- z (clamp z 90 110))) 3.99) (- 11.005 y) (- y 18.005)))
(body key (shape key-field) (attr material "steel") (attr color "#404040"))

(def P 3000)
(def T 100000)
(def seg-c 12)
(def xbar 13.210919583)
(def Ft (/ T (* 2 xbar)))
(fn seat (zc) (max (- (abs (- z zc)) 1) (- 8 (len2 x y))))
(study drive (shape shaft-field) (elastic (material steel))
  (fix brg-l (seat 10) (pivot (at 0 0 10)))
  (fix brg-r (seat 190) (pivot (at 0 0 190)) (slide 0 0 1))
  (fix coupling (- z -19) (dofs y) (support))
  (load belt (force (* 1N (neg P)) 0 0) (region (- (abs (- z 100)) 5)))
  (load torque-a (force 0 (* 1N Ft) 0) (region (max (- (abs (- z 100)) 5) (- seg-c x))))
  (load torque-b (force 0 (* 1N (neg Ft)) 0) (region (max (- (abs (- z 100)) 5) (+ x seg-c))))
  (res 1))

(require shoulder-peak-over-155 (max-stress drive > 155MPa))
(require shoulder-peak-band (max-stress drive < 1e6MPa))

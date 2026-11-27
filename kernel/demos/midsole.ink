; A 3D-printed running-shoe midsole (men's US 9: 285 mm long, 104 across
; the forefoot), its core three lattices that turn into one another along
; the foot: a diamond lattice under the heel, where the foot strikes; a
; gyroid under the arch, stiffer in twist; a Schwarz P under the forefoot,
; softer, to flex and spring. Each hand-off is a crossfade of two fields
; over 30 mm, so the cells morph from one kind to the next with no seam.
; The walls thicken under the heel and the ball of the foot (a pressure
; map: two pads, a design heuristic, not a gait measurement), and the cells
; follow the sole: their layers bend up with the toe spring and the heel's
; bevel (a conformal lattice: its height coordinate is height above the
; bottom of the sole). Heel stack 32 mm, forefoot 24 (8 mm drop).
(import "std")

(param drop 8)          ; heel stack minus forefoot stack, mm
(param fore-stack 24)   ; forefoot stack height, mm
(param spring 14)       ; toe spring: how far the toe lifts, mm

; ---------------------------------------------------------------- outline
; The footprint (a right foot, heel at x = 0, toe at x = 285, the arch on
; +y): a smooth closed curve (Catmull-Rom) through 20 points, 78 mm across
; the heel, 70 at the waist and 104 across the ball.
(def foot-pts (vec (vec 0 0) (vec 6 -26) (vec 30 -40) (vec 70 -38) (vec 120 -40) (vec 170 -46)
                   (vec 205 -50) (vec 240 -44) (vec 268 -28) (vec 283 -8) (vec 285 6) (vec 276 26)
                   (vec 250 46) (vec 205 54) (vec 165 44) (vec 125 30) (vec 85 32) (vec 45 38)
                   (vec 15 32) (vec 2 16)))
(def foot-n (len foot-pts))
(def foot-loop
  (apply concat (map (lambda (i)
    (map (lambda (s) (catmull-rom (nth foot-pts (- i 1)) (nth foot-pts i)
                                  (nth foot-pts (mod (+ i 1) foot-n)) (nth foot-pts (mod (+ i 2) foot-n)) s))
         (vec 0 0.2 0.4 0.6 0.8)))
    (range foot-n))))
(def outline (polygon foot-loop))

; ---------------------------------------------------------------- profile
; The sole's bottom: flat from 45 to 190 mm, rising in a rocker to the toe
; (toe spring) and a small bevel at the heel; the stack on it falls from
; heel to forefoot and thins at the toe.
(def toe-u (clamp (/ (- x 190) 95) 0 1))
(def heel-u (clamp (/ (- 45 x) 45) 0 1))
(def bottom (+ (* spring (square toe-u)) (* 5 (square heel-u))))
(def stack (- (+ fore-stack (* drop (clamp (/ (- 170 x) 110) 0 1))) (* 8 (square toe-u))))
; The footbed is cupped: it rises 5 mm toward the rim at the heel and arch.
(def rim (clamp (/ (+ outline 14) 14) 0 1))
(def cup (* 5 (square rim) (clamp (/ (- 170 x) 60) 0 1)))
(def top (+ bottom stack cup))
(def solid (smooth-intersect outline (max (- bottom z) (- z top)) 3))

; ---------------------------------------------------------------- lattices
; Conformal coordinates: height above the sole's bottom, so cells bend with
; the rocker. Each TPMS divided by its greatest gradient.
(def wz (- z bottom))
(fn tpms-k (p) (/ (* 2 pi) p))
(def diamond-p 11)
(def gyroid-p 10)
(def schwarz-p-p 9)
(def lat-diamond
  (let ((k (tpms-k diamond-p))
        (sx (sin (* k x))) (sy (sin (* k y))) (sz (sin (* k wz)))
        (cx (cos (* k x))) (cy (cos (* k y))) (cz (cos (* k wz))))
    (/ (+ (* sx sy sz) (* sx cy cz) (* cx sy cz) (* cx cy sz)) (* k (sqrt 3)))))
(def lat-gyroid (warp (gyroid gyroid-p) x y wz))
(def lat-schwarz
  (let ((k (tpms-k schwarz-p-p)))
    (/ (+ (cos (* k x)) (cos (* k y)) (cos (* k wz))) (* k (sqrt 3)))))

; The zones: heel to arch over 75..105 mm, arch to forefoot over 140..170.
(fn ease (u) (let ((c (clamp u 0 1))) (* c c (- 3 (* 2 c)))))
(def s1 (ease (/ (- x 75) 30)))
(def s2 (ease (/ (- x 140) 30)))
(def lattice (blend (blend lat-diamond lat-gyroid s1) lat-schwarz s2))

; Wall thickness: 0.6 mm, up to 1.6 mm on the two pressure pads.
(def pad (max (exp (neg (/ (+ (square (- x 44)) (square y)) (square 30))))
              (exp (neg (/ (+ (square (- x 198)) (square (- y 10))) (square 36))))))
(def wall (+ 0.6 (* 1.0 pad)))
(def core (- (abs lattice) (* 0.5 wall)))

; ---------------------------------------------------------------- skins
; A 1.6 mm footbed on top and a 2 mm outsole beneath, its forefoot cut by
; six 1 mm flex grooves, the lattice filleted into both. The outsole's
; edge is rounded 0.5 mm: a 1.5 mm round top and bottom would leave the
; 2 mm sole no side wall.
(def footbed (max solid (- (- top 1.6) z)))
(def grooves (max (- (abs (- (mod (- x 160) 18) 9)) 0.9) (- 160 x) (- x 270) (- z (+ bottom 1))))
(def outsole (subtract (smooth-intersect outline (max (- bottom z) (- z (+ bottom 2))) 0.5) grooves))
(def midsole-field (max solid (smooth-union (smooth-union core footbed 1.2) (max solid (- z (+ bottom 2.4))) 1.2)))

; The lattice coloured by zone (heel teal, arch violet, forefoot coral), the
; footbed white.
(def zone-color (mix (mix (vec 0.1 0.62 0.62) (vec 0.45 0.32 0.85) s1) (vec 0.98 0.42 0.35) s2))
; A step, which the mesh cuts as a sharp line.
(def paint (select (- z (- top 1.8)) zone-color (vec 0.94 0.94 0.92)))

(body midsole
  (shape (max midsole-field (- (+ bottom 2) z)))
  (attr material-note "printed elastomer (EPU / TPU)")
  (attr color paint)
  (attr metallic 0)
  (attr roughness 0.6))

(body outsole
  (shape outsole)
  (attr material-note "rubber")
  (attr color (vec 0.12 0.12 0.13))
  (attr metallic 0)
  (attr roughness 0.85))

; The same midsole without its footbed (for pictures): the three lattices
; seen from above, turning into one another from heel to toe.
(body midsole-open
  (shape (max solid (- z (- top 1.6)) (- (+ bottom 2) z)
              (smooth-union core (max solid (- z (+ bottom 2.4))) 1.2)))
  (attr color zone-color)
  (attr metallic 0)
  (attr roughness 0.6))

(assembly shoe
  (instance midsole midsole (at))
  (instance outsole outsole (at)))

; Cut along its length (the medial half, y > 6, removed) to show the cells
; bending with the sole and thickening on the pads (for pictures).
(def cut (- 6 y))
(body midsole-cut
  (shape (max (max midsole-field (- (+ bottom 2) z)) cut))
  (attr color paint)
  (attr metallic 0)
  (attr roughness 0.6))
(body outsole-cut
  (shape (max outsole cut))
  (attr color (vec 0.12 0.12 0.13))
  (attr metallic 0)
  (attr roughness 0.85))
(assembly shoe-cut
  (instance midsole midsole-cut (at))
  (instance outsole outsole-cut (at)))

(assembly shoe-open
  (instance midsole midsole-open (at))
  (instance outsole outsole (at)))

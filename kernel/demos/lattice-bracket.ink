; A lightweight cantilever bracket for metal printing: four M8 bolts at
; its foot, a 16 mm clevis pin 85 mm above, and between them a smooth,
; organic body that is a 1.6 mm skin over a gyroid core. The core's density
; is graded by distance to where the loads come in: 70 % solid around the
; bolts and the pin, falling to 20 % far from them, with every wall at
; least 0.55 mm thick. The grading is one field (a wall thickness that
; reads distance), so the lattice has no cells to place, trim or join.
(import "std")

(param period 10)      ; gyroid cell, mm
(param skin 1.6)       ; outer skin, mm
(param pin-d 16)       ; clevis pin, mm
(param bolt-d 9)       ; M8 clearance, mm
(param dense 0.7)      ; core density at the load points
(param sparse 0.2)     ; core density far from them
(param reach 28)       ; distance over which the density falls, mm

(def pin-z 85)
(def bolt-x 46)
(def bolt-y 24)

; The envelope: a base plate, a tapering web and the pin's boss, joined by
; large smooth blends.
(def plate (move (rounded-box 120 70 12 4) 0 0 6))
; The web's profile is in (x, z), extruded 32 mm along y.
(def web (warp (extrude-z (polygon (vec (vec -54 6) (vec 54 6) (vec 16 pin-z) (vec -16 pin-z))) -16 16) x z y))
(def boss (move (cylinder-y 19 52) 0 0 pin-z))
(def envelope (smooth-union (smooth-union plate web 14) boss 12))

; Distance to the load paths: the bolts' axes and the pin's.
(def to-bolts (len2 (- (abs x) bolt-x) (- (abs y) bolt-y)))
(def to-pin (len2 x (- z pin-z)))
(def to-load (min (- to-bolts (* 0.5 bolt-d)) (- to-pin (* 0.5 pin-d))))

; The core: a gyroid sheet whose thickness gives the density there; solid
; within 3 mm of a hole, where the bolt heads and the pin bear.
(def rho (mix dense sparse (clamp (/ to-load reach) 0 1)))
(def th (gyroid-thickness rho period))
(def g (gyroid period))
(def bearing (max envelope (- to-load 3)))
(def outer (max envelope (neg (+ envelope skin))))
; The holes: bolt clearances through the plate (spot-faced 16 mm on top),
; and the pin's bore.
(def bolts (min (- to-bolts (* 0.5 bolt-d)) (max (- to-bolts 8) (- 12 z))))
; The bracket with core walls `wall` thick (th, or any other field).
(fn bracket-with (wall)
  (let ((core (max (+ envelope skin) (- (abs g) (* 0.5 wall))))
        (solid (min outer bearing (fillet-union core outer 0.8 90))))
    (subtract solid bolts (move (cylinder-y (* 0.5 pin-d) 80) 0 0 pin-z))))
(def bracket-field (bracket-with th))

; Skin bare metal; the core coloured by its density (blue sparse, red dense).
(def paint-core (color-map rho (- sparse 0.12) dense))
(def paint-skin (vec 0.72 0.73 0.75))
(def colour (mix paint-core paint-skin (clamp (+ 0.5 (/ (+ envelope skin) 0.4)) 0 1)))

(body bracket
  (shape bracket-field)
  (attr material-note "AlSi10Mg, laser powder bed fusion")
  (attr density 2.67g/cm3)
  (attr color colour)
  (attr metallic 0.5)
  (attr roughness 0.4))

; The same, cut open along y = 0 to show the graded core (for pictures).
(body bracket-cut
  (shape (max bracket-field (neg y)))
  (attr color colour)
  (attr metallic 0.5)
  (attr roughness 0.4))

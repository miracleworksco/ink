; The lattice bracket as a morph in time: a solid bracket hollowed into
; its graded gyroid core by a front that sweeps down it and back, cut open
; along y = 0. One field of x, y, z and t: the core's wall thickness is
; blended from solid to graded across the front, so as it moves the part
; opens and closes hundreds of tunnels, a change of topology at every
; step that a boundary model would have to rebuild from scratch.
(import "std")
(import "lattice-bracket.ink" as b)

(t-span 0 4s)

; The front: from above the pin (z 125) down below the plate (z -15) over
; 2 s and back, eased; walls go from solid (3.2 mm, more than the gyroid's
; whole swing) to graded over 12 mm behind it.
(def phase (* 0.5 (- 1 (cos (* 90deg/s t)))))
(def front (mix 125 -15 phase))
(def k (clamp (/ (- z front) 12) 0 1))
(def wall (mix 3.2 b.th k))

(body bracket-morph
  (shape (max (b.bracket-with wall) (neg y)))
  (attr color (mix b.paint-skin b.colour k))
  (attr metallic 0.5)
  (attr roughness 0.4))

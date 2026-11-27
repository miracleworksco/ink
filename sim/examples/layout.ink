; The stiffest bracket that uses 40% of a PETG plate, held at a wall and
; loaded at the bottom of its free end.
; ink topology examples/layout.ink: the layout found and its compliance
; (load times deflection, N mm: lower is stiffer) at each step;
; ink mesh examples/layout.ink --target bracket --out bracket.stl to print it.
(import "std")

(def domain (move (box 60 2 24) 30 0 0))         ; x from 0 to 60, z from -12 to 12
(study hang (shape domain) (elastic (material petg))
  (fix wall (- x 1))
  (load tip (force 0 0 -100N) (region (max (- 58 x) (+ z 9))))
  (res 1))
(topology layout (study hang) (volume-fraction 0.4))

; The layout as a part: material where its density is over 1/2.
(body bracket (shape (max domain (study-field layout solid))))

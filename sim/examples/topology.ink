; Topology optimization: the stiffest cantilever that uses 35% of a thin
; plate, fixed at a wall and loaded at its free end. The result is a field,
; so it is checked, filleted and exported like any other shape. (In a thick
; block the same load gives a hollow box beam: webs on both faces.)
(import "std")

(def domain (move (box 120 4 50) 60 0 0))
(body envelope (shape domain))

(study cantilever (shape domain) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (- x 1))
  (load tip (force 0 0 -200N) (region (max (- 118 x) (- (abs z) 4))))
  (res 1))
(topology layout (study cantilever) (volume-fraction 0.35) (iterations 60) (filter 3))

; The optimized part: material where the density exceeds 1/2.
(def optimized (max domain (study-field layout solid)))
(body bracket (shape optimized) (attr color "#6C8EBF"))

; Verify the result with an ordinary study of the new shape.
(study verify (shape (shape-of bracket)) (elastic (E 2000MPa) (nu 0.3))
  (fix wall (- x 1))
  (load tip (force 0 0 -200N) (region (max (- 118 x) (- (abs z) 4))))
  (res 0.5))

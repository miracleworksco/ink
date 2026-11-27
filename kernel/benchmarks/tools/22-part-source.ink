; Source part for benchmark 22: a 60 x 40 x 10 plate with a counterbored
; hole (bore r 8 through, counterbore r 12, 4 deep from the top).
; Benchmark 22's input, benchmarks/inputs/22-part.stl, is its mesh:
;   ink mesh benchmarks/tools/22-part-source.ink --res 0.2 --out benchmarks/inputs/22-part.stl
(import "std")
(def plate (move (box 60 40 10) 0 0 5))
(def bore (cylinder-z 8 100))
(def counterbore (move (cylinder-z 12 8) 0 0 10))
(body part (shape (subtract plate (min bore counterbore))))

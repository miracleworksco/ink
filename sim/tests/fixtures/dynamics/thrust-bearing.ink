; A single-direction thrust bearing (51100: 24 x 10 x 9 mm): it pushes
; along its +z only.
(import "std")
(def spec (get (read-json "ratings.json") "thrust"))
(body bearing (shape (subtract (cylinder-z 12 9) (cylinder-z 5 10))) (frame load (at)))

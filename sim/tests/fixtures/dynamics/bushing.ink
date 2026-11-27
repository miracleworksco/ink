; Linear ball bushings, carrying across their bores: a plain one (LM8UU)
; at its middle, and a flanged one (LMF10LUU: 55 long, flange 6) at its
; middle, L/2 - H = 21.5 mm along its bore from its flange's face.
(import "std")
(def spec (get (read-json "ratings.json") "bushing"))
(body bushing (shape (subtract (cylinder-z 7.5 24) (cylinder-z 4 30))) (frame load (at)))
(body flanged (shape (subtract (min (move (cylinder-z 9.5 55) 0 0 27.5) (move (cylinder-z 19.5 6) 0 0 -3)) (cylinder-z 5 130)))
  (frame load (at (move 0 0 21.5))))

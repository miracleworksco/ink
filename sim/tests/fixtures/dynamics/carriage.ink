; A guide carriage (MGN12H: 54 g, C0 5880 N, MR 38.22 N m) and its rail:
; the carriage carries on its rail under its centre, x along it.
(import "std")
(def spec (get (read-json "ratings.json") "carriage"))
(body carriage (shape (move (box 45 27 10) 0 0 5)) (attr mass 54g)
  (frame load (at))
  (frame base (at (rot-x 180) (move 0 0 3))))
(body rail (shape (move (box 300 12 8) 150 0 -4)))

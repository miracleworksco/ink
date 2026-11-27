; A shorter block for the same rail (MGN12C): ratings of its own.
(import "std")
(def spec (get (read-json "ratings.json") "carriage-short"))
(body carriage (shape (move (box 35 27 10) 0 0 5)) (attr mass 34g) (frame load (at)))

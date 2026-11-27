; A part whose record rates no load.
(import "std")
(def spec (get (read-json "ratings.json") "screw"))
(body part (shape (cylinder-z 1.5 8)))

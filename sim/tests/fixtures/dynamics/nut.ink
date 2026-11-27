; A ball screw's nut (SFU1605): it carries along its axis, either way.
(import "std")
(def spec (get (read-json "ratings.json") "nut"))
(body nut (shape (subtract (cylinder-z 14 40) (cylinder-z 8 50))) (frame load (at)))

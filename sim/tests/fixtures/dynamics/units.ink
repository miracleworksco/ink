; Mounted bearing units, carrying at their inserts' centres: a flange unit
; (UCF201: 15 mm in from its mounting face), a pillow block (UCP201: its
; shaft along its y, 30.2 mm up), and a unit whose load frame is its
; assembly's (4461N18: its outer ring spans 1.8 to 6.8 mm, so 4.3).
(import "std")
(def spec (get (read-json "ratings.json") "unit"))
(body flange (shape (subtract (move (box 86 86 30) 0 0 15) (cylinder-z 6.5 70))) (frame load (at (move 0 0 15))))
(body pillow (shape (subtract (move (box 127 38 46) 0 0 23) (place (cylinder-z 6.5 60) (at (rot-x -90) (move 0 0 30.2)))))
  (frame load (at (rot-x -90) (move 0 0 30.2))))
(body ring (shape (move (subtract (cylinder-z 14 5) (cylinder-z 10 6)) 0 0 4.3)))
(body plate (shape (move (subtract (cylinder-z 30 1.5) (cylinder-z 10 2)) 0 0 0.75)))
(assembly mounted (instance plate plate) (instance ring ring) (frame load (at (move 0 0 4.3))))

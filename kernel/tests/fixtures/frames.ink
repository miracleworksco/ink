(import "std")
(body plate (shape (box 20 20 4)) (frame top (at (move 0 0 2))))
(body pin (shape (box 4 4 6)) (frame seat (at (rot-x 180) (move 0 0 -3))))
(assembly root (instance base plate) (instance peg pin (frame-of base top) seat))

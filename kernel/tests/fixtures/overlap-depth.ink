; Two blocks overlapping 2.5 mm: `ink interfere` names the overlap and gives its `depth`, 1.25 mm at
; its middle (the first point found lay 0.94 deep, and a camera slider's 2.5 mm overlap read 0.054).
(import "std")
(t-span 0 1s)
(body a (shape (box 40 20 10)))
(body b (shape (box 10 10 10)))
(assembly root (instance a a) (instance b b (at (move 22.5 0 0))))

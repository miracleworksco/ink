; A part with params an importer gives values: its length and two taps.
(import "std")
(param length 100)
(param tap-start 0)
(param tap-end 0)
(body extrusion (shape (move (box 20 20 length) 0 0 (/ length 2)))
  (frame end-a (at (rot-x 180)))
  (frame end-b (at (move 0 0 length))))

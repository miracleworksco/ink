; A radial ball bearing (608: 22 x 8 x 7 mm, C 3300 N, C0 1370 N) as the
; dynamics package reads a support: its record's load-ratings and where
; they act, its load frame (z along the bore), `centre` along the bore
; from its middle (an angular-contact bearing's pressure centre).
(import "std")
(param centre 0)
(def spec (get (read-json "ratings.json") "radial"))
(body (bearing (centre centre)) (shape (subtract (cylinder-z 11 7) (cylinder-z 4 8))) (attr material "steel")
  (frame load (at (move 0 0 centre)))
  (frame face-a (at (move 0 0 3.5)))
  (frame face-b (at (rot-x 180) (move 0 0 -3.5))))

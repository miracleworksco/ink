; A turn's angle wrapped by mod at a whole turn, or chosen by a select
; between branches a whole turn apart where it switches, moves the turn by
; none: its speed limit reads the turn (its sine and cosine), continuous
; there exactly, so it holds. Wrapped at 180 the turn jumps half round:
; unknown at the wrap. Written unwrapped, (* 360 t), it holds.
(import "std")
(t-span 0 1.5)
(body arm (shape (box 10 2 2)))
(assembly root
  ; Speeding up (400 t deg/s), chosen a whole turn back at t = 1.
  (instance chosen arm (at (move 0 0 0)) (turn z (select (- t 1) (* 200 t t) (- (* 200 t t) 360))) (speed-limits -700 700))
  ; Speeding up (400 t deg/s), wrapped at 360 at t = 1.34.
  (instance wrapped arm (at (move 0 20 0)) (turn z (mod (* 200 t t) 360)) (speed-limits -700 700))
  ; Wrapped by 180: a half turn at t = 0.5, on a step's end, moves the
  ; arm's end 10 mm.
  (instance flipped arm (at (move 0 40 0)) (turn z (mod (* 360 t) 180)) (speed-limits -400 400)))
; expect-status: speed-limits:root/flipped unknown

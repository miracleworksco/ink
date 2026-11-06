; A turn's angle wrapped by a select or by mod jumps by 360 (or 180): its
; speed limit reads the drive, and a jump is unbounded speed, so it is
; unknown at the wrap, though the turn itself is continuous there. Written
; unwrapped, (* 360 t), it holds.
(import "std")
(t-span 0 1.5)
(body arm (shape (box 10 2 2)))
(assembly root
  (instance chosen arm (at (move 0 0 0)) (turn z (select (- t 1) (* 360 t) (- (* 360 t) 360))) (speed-limits -400 400))
  (instance wrapped arm (at (move 0 20 0)) (turn z (mod (* 360 t) 360)) (speed-limits -400 400))
  ; Wrapped by 180: a half turn at t = 0.5, on a step's end, moves the
  ; arm's end 10 mm.
  (instance flipped arm (at (move 0 40 0)) (turn z (mod (* 360 t) 180)) (speed-limits -400 400)))
; expect-status: speed-limits:root/chosen unknown
; expect-status: speed-limits:root/wrapped unknown
; expect-status: speed-limits:root/flipped unknown

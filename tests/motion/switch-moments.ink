; Switches read by value. A trapezoid's speed written as clamped pieces
; ends one ramp's clamp (its condition 0.5 - max(0, t - 0.25)) as the
; next one's starts (t - 0.75): two conditions, one moment, so each case
; pairs their sides at once and the speed, 500 deg/s either side, holds
; to 600. Written as a select chain whose branches are one number at
; each switch (t1 = 0.25, exact), it is continuous, its speed decided;
; one branch a degree up at t = 0.75 jumps there, and is unknown.
(import "std")
(t-span 0 1)
(body arm (shape (box 10 2 2)))
(def acc 2000)
(def v 500)
(def t1 (/ v acc))
(def t2 (- 1 t1))
(def ta (min t t1))
(def tc (clamp (- t t1) 0 (- 1 (* 2 t1))))
(def td (clamp (- t t2) 0 t1))
(def clamped (+ (* 0.5 acc ta ta) (* v tc) (- (* v td) (* 0.5 acc td td))))
(def p1 (* 0.5 acc t1 t1))
(def p2 (+ p1 (* v (- t2 t1))))
(fn chain (lift) (select (- t t1) (* 0.5 acc t t)
                   (select (- t t2) (+ p1 (* v (- t t1)))
                     (- (+ p2 lift (* v (- t t2))) (* 0.5 acc (- t t2) (- t t2))))))
(assembly root
  (instance clamped arm (at (move 0 0 0)) (turn z clamped) (speed-limits -600 600))
  (instance chain arm (at (move 0 20 0)) (turn z (chain 0)) (speed-limits -510 510))
  (instance lifted arm (at (move 0 40 0)) (turn z (chain 1)) (speed-limits -510 510)))
; expect-status: speed-limits:root/lifted unknown

; A motor's rotor and a line a third as fast, driven as a roller
; floor's are: by a speed pulse written w tau (h + ...), tau = w / 4a a
; float. Their rates in the pulse's expression, -w tau and w tau / 3, are
; each a ratio of integers past 2^53 (tau's digits) but exactly -3 apart,
; so the pair repeats every 1080 / (w tau) of it, the value all its reads
; of t pass through. The arms sweep across each other but pass 4.7 mm
; apart: proven over one period of it, and so over the eleven seconds.
(import "std")
(t-span 0 11)
(fn lncosh (z) (ln (* 0.5 (+ (exp z) (exp (neg z))))))
(fn pulse (t0 t1 w a)
  (let ((tau (/ w (* 4 a))) (h (/ (- t1 t0) (* 2 tau))) (s (/ (- t (* 0.5 (+ t0 t1))) tau)))
    (* w tau (+ h (lncosh (* 0.5 (+ s h))) (neg (lncosh (* 0.5 (- s h))))))))
(def drive (pulse 3.0 6.0 5800rpm (* 200 360)))
(body arm (shape (box 30 2 2)))
(assembly line (instance arm arm (at (move 25 0 0))))
(assembly root
  (instance line line (at) (turn z (/ drive 3)))
  (instance rotor line (at (move 50 0 0)) (turn z (neg drive))))
(assert apart (disjoint (shape-of root/line/arm) (shape-of root/rotor/arm)))

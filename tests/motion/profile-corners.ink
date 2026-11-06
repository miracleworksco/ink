; Motion profiles written with min and clamp on t are continuous by
; construction: their speed limits are decided with t ranged, each case of
; a switch the speed reads several times bounded apart. A trapezoid
; written as a select chain may jump where it switches (its branches meet
; in real numbers but not in floats, t1 = 1/3): its speed limits are
; unknown.
(import "std")
(t-span 0 1)
(body arm (shape (box 10 2 2)))

; A trapezoidal speed: 1500 deg/s^2 to 500 deg/s, cruise, back to rest.
(def acc 1500)
(def v 500)
(def t1 (/ v acc))
(def t2 (- 1 t1))
(def p1 (* 0.5 acc t1 t1))
(def p2 (+ p1 (* v (- t2 t1))))
(def trapezoid (select (- t t1) (* 0.5 acc t t)
                 (select (- t t2) (+ p1 (* v (- t t1)))
                   (- (+ p2 (* v (- t t2))) (* 0.5 acc (- t t2) (- t t2))))))

; A cycloid of 90 degrees over 0.8 s, then a dwell.
(def u (min (/ t 0.8) 1))
(def cycloid (* 90 (- u (/ (sin (* 2 pi u)) (* 2 pi)))))

; Eased keyframes (an agent's robot arm, 2026-10-25): a smoothstep of each
; segment's clamped phase.
(fn ease (s) (* s s (- 3 (* 2 s))))
(fn seg (t0 t1) (ease (clamp (/ (- t t0) (- t1 t0)) 0 1)))
(def keys (+ (* -40 (seg 0 0.3)) (* 80 (seg 0.45 0.7)) (* -40 (seg 0.85 1))))

(assembly root
  (instance a arm (at (move 0 0 0)) (turn z trapezoid) (speed-limits -500 500 (accuracy 0.01)))
  (instance exact arm (at (move 0 60 0)) (turn z trapezoid) (speed-limits -510 510))
  (instance b arm (at (move 0 20 0)) (turn z cycloid) (speed-limits -10 230))
  (instance c arm (at (move 0 40 0)) (turn z keys) (speed-limits -500 500)))
; expect-status: speed-limits:root/a unknown
; expect-status: speed-limits:root/exact unknown

; A spin-up a read through two waves with no common period (the sine of a
; in radians, the cosine of a in degrees), which drive an outer turn b. The
; probe reads t only through b, repeating every turn of it: proven over
; one turn, whatever a does. The two waves summed alone repeat every turn
; of a (1e7 turns over the span): proven over one.
(t-span 0 1)
(def a (* 3600000000 (square t)))
(def b (* 1000000000 (+ (sin a) (cos (deg->rad a)))))
(def s (sin (deg->rad b)))
(def f (- (* (- x s) (+ x s)) 0.1))
(assert inside (subset (probe 0 0 0) f))
(def w (+ (sin (deg->rad a)) (cos (deg->rad a))))
(assert bounded (value w <= 1.5))

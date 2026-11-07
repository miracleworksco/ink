; A block slides at 100 mm/s, driven through a solve (q = 100 t), under a
; 50 mm/s speed limit: the limit reads q's rate through the solve (the
; implicit function's, 100), bounded over t, so it fails.
(import "std")
(t-span 0 1)
(solve position (unknowns (q 0)) (== q (* 100 t)))
(body block (shape (box 2 2 2)))
(assembly root (instance block block (at (move 0 0 0)) (slide x q) (speed-limits -50 50)))
(assert rate (value (diff q t) > 99.999))
; Three consistent equations for two unknowns (least squares): their rates
; are the least-squares solve's, whose root is certified nowhere: measured.
(solve pair (unknowns (a 0) (b 0)) (== a (* 10 t)) (== b (* 20 t)) (== (+ a b) (* 30 t)))
(assert pair-rates (value (min (- (diff a t) 9.999) (- (diff b t) 19.999)) > 0))
; expect-status: speed-limits:root/block fails
; expect-status: rate holds
; expect-status: pair-rates holds estimate

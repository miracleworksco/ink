; A value over the motion: 10 sin(2 pi t) peaks at 10 at t = 0.25.
(import "std")
(t-span 0 1)
(def w (* 10 (sin (* 2 pi t))))
(assert under-peak (value w < 10.01))
(assert at-peak (value w < 9.99))
; expect-status: at-peak fails

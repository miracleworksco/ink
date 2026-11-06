; Benchmark 32 companion: the loop closure solved numerically (solve, warm
; started at each --t), against the closed form the benchmark's placements use.
; A t-dependent solve blocks --t-span for the whole file, so it lives here:
;   for t in 0 0.11611 0.25 0.5 0.66387 0.75: ink assert tools/32-solved.ink --t $t --name solved-matches-closed-form
; Each holds as an estimate (accuracy null): solved values are numerical
; inputs.
; expect-warning: the file reads t but states no t-span
(import "../32-four-bar.ink")
(solve loop
  (unknowns (th3s 0.97) (th4s 1.70))
  (== (+ bx (* b (cos th3s))) (+ d (* c (cos th4s))))
  (== (+ by (* b (sin th3s))) (* c (sin th4s))))
(assert solved-matches-closed-form
  (value (max (abs (- (* r2d th4s) rocker-deg)) (abs (- (* r2d th3s) coupler-deg))) < 1e-6))
(assert solved-rocker-in-range (value (* r2d th4s) > 88.9))

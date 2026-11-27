; Jumps in a drive: a jump is unbounded speed, so a speed limit on a
; drive that may jump (a select or floor switching on t, an atan2's wrap)
; is unknown there, never held: however small the jump (a nudge, a blip),
; and whatever accuracy it is asked to. The same drive without one (a
; pulse or a kink written with clamp and max, a select whose branches meet
; where it switches) is decided: a pulse past the limit fails at a moment
; inside it.
(import "std")
(t-span 0 1s)
(body arm (shape (box 10 2 2)))
(assembly root
  ; 100 deg/s, 400 deg/s from t = 0.53 to 0.54, then 100 again.
  (instance pulse arm (at (move 0 0 0)) (turn z (+ (* 100deg/s t) (* 300deg/s (clamp (- t 0.53s) 0 0.01s)))) (speed-limits -300deg/s 300deg/s))
  ; Still, then 100 deg/s from t = 0.5.
  (instance kink arm (at (move 0 20 0)) (turn z (* 100deg/s (max (- t 0.5s) 0))) (speed-limits -300deg/s 300deg/s))
  ; At 0, then at 30 degrees from t = 0.5.
  (instance jump arm (at (move 0 40 0)) (turn z (select (- t 0.5s) 0 30)) (limits -90 90) (speed-limits -300deg/s 300deg/s))
  ; 10 deg/s, nudged 0.004 degrees at t = 0.53: 0.0004 mm at the arm's end.
  (instance nudge arm (at (move 0 60 0)) (turn z (select (- t 0.53s) (* 10deg/s t) (+ 0.004 (* 10deg/s t)))) (speed-limits -300deg/s 300deg/s (accuracy 0.01deg/s)))
  ; 50 deg/s, then 100 deg/s from t = 0.5, written as a select: 25
  ; degrees both sides of the switch.
  (instance switch arm (at (move 0 80 0)) (turn z (select (- t 0.5s) (* 50deg/s t) (+ 25 (* 100deg/s (- t 0.5s))))) (speed-limits -300deg/s 300deg/s))
  ; 200 steps of 0.008 degrees, 1.6 in all (0.14 mm at the arm's end).
  (instance stairs arm (at (move 0 100 0)) (turn z (* 0.008 (floor (* 200 (/ t 1s))))) (speed-limits -1deg/s 1deg/s (accuracy 0.01deg/s)))
  ; 1000 steps of 0.009 mm, 9 mm in all.
  (instance treads arm (at (move 0 120 0)) (slide x (* 0.009 (floor (* 1000 (/ t 1s))))) (speed-limits -1mm/s 1mm/s (accuracy 0.01mm/s)))
  ; 100 mm/s with a 4 mm step at t = 0.5, its speed asked to 5 mm/s.
  (instance hop arm (at (move 0 140 0)) (slide x (select (- t 0.5s) (* 100mm/s t) (+ 4 (* 100mm/s t)))) (speed-limits -500mm/s 500mm/s (accuracy 5mm/s)))
  ; 3000 rpm with a 90 degree step at t = 0.5, its speed asked to 100 deg/s.
  (instance flip arm (at (move 0 160 0)) (turn z (select (- t 0.5s) (* 18000deg/s t) (+ 90 (* 18000deg/s t)))) (speed-limits 0 6000rpm (accuracy 100deg/s)))
  ; 10 mm/s, out 0.03 mm and back within 0.1 us across t = 0.5.
  (instance blip arm (at (move 0 180 0)) (slide x (+ (* 10mm/s t) (select (- t 0.49999995s) 0 (select (- t 0.50000005s) 0.03 0)))) (speed-limits -500mm/s 500mm/s))
  ; 10 mm/s, then 20 from t = 0.5, stepping 0.06 mm there.
  (instance kick arm (at (move 0 220 0)) (slide x (select (- t 0.5s) (* 10mm/s t) (+ 5.06 (* 20mm/s (- t 0.5s))))) (speed-limits -25mm/s 25mm/s))
  ; Two 0.04 mm steps (none past the accuracy) and a 0.07 mm one, in one
  ; starting step: the search that stops at the pair still finds the 0.07.
  (instance trio arm (at (move 0 240 0)) (slide x (+ (* 10mm/s t) (select (- t 0.505s) 0 0.04) (select (- t 0.515s) 0 0.04) (select (- t 0.55s) 0 0.07))) (speed-limits -20mm/s 20mm/s))
  ; Never 1.3 um from 0, its formula speed 0.4 mm/s: 159 jumps back.
  (instance saw arm (at (move 0 200 0)) (slide x (* 0.0004 (atan2 (sin (* 1000rad/s t)) (cos (* 1000rad/s t))))) (speed-limits 0.3mm/s 0.5mm/s)))
; expect-status: speed-limits:root/pulse fails
; expect-status: speed-limits:root/jump unknown
; expect-status: speed-limits:root/nudge unknown
; expect-status: speed-limits:root/stairs unknown
; expect-status: speed-limits:root/treads unknown
; expect-status: speed-limits:root/hop unknown
; expect-status: speed-limits:root/flip unknown
; expect-status: speed-limits:root/blip unknown
; expect-status: speed-limits:root/kick unknown
; expect-status: speed-limits:root/trio unknown
; expect-status: speed-limits:root/saw unknown

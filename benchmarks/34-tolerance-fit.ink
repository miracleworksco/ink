; Benchmark 34: ISO 286 fits at 8 mm, proven over their tolerance bands.
;
; Task: a running fit and a press fit, each a shaft through a 10 mm plate,
; with the diameters' tolerances declared on their params, so every claim
; is checked for every pair of sizes in the bands, not at the nominal 8.
; ISO 286 at 6-10 mm: IT6 9 um, IT7 15 um, IT8 22 um; shaft deviations
; f es -13 um, s ei +23 um; holes H, EI 0. So:
;   - running fit H8/f7: bore 8.000..8.022, shaft 7.972..7.987;
;   - press fit H7/s6: bore 8.000..8.015, pin 8.023..8.032.
; And a control: the loose spec a first draft gave, shaft 8 +-0.05 in a
; bore 8.00..8.10, which clears at the middle sizes but not at all.
;
; Hand calculation (radial, half the diametral difference):
;   - running fit: clearance (D - d) / 2 in [(8.000 - 7.987) / 2,
;     (8.022 - 7.972) / 2] = [0.0065, 0.025]; at the middle sizes (8.011,
;     7.9795) 0.01575. Diametral: D - d in [0.013, 0.050].
;   - press fit: interference (d - D) / 2 in [(8.023 - 8.015) / 2,
;     (8.032 - 8.000) / 2] = [0.004, 0.016]. A ring 4.0075 < r < 4.0115
;     (the largest bore radius to the smallest pin radius) is inside both
;     at every pair of sizes, so they overlap throughout: the pin is larger
;     than the bore at every pair. ink proves no overlap over a band: a
;     point found at one pair says nothing of the others (docs/CORE.md).
;   - loose spec: (D - d) / 2 in [(8.00 - 8.05) / 2, (8.10 - 7.95) / 2] =
;     [-0.025, 0.075]: a shaft 8.05 in a bore 8.00 interferes by 0.025.
;
; Pass criteria (ink assert, bands applied by default; ink interfere):
;   - running-clear-5um, running-gap-12um, press-overlaps: proven;
;   - the controls disproven, each with a witness in the bands that fails
;     the claim: running-clear-10um (a pair of points closer than 0.010 at
;     sizes with (D - d) / 2 < 0.010), running-gap-20um (D - d <= 0.020),
;     running-overlaps (the shaft smaller than the bore at some sizes),
;     loose-apart (a point inside
;     both, at sizes with d > D);
;   - pinned at the band middles (--param bore-d=8.011:8.011,shaft-d=
;     7.9795:7.9795) running-clear-10um is proven, so the tolerances are
;     what fails it. The basic sizes 8 / 8 would not do: an f7 shaft's 8
;     lies outside its own band, so each value is written in its band; and
;     8 in 8 touches;
;   - interfere: the press fit is the only overlap, as declared.

; The press fit stands beside the running fit.

(import "std")

; ---------------------------------------------------------------- running fit H8/f7
(param bore-d 8 (tol 0 0.022))          ; H8: 8.000..8.022
(param shaft-d 7.9795 (tol -0.0075 0.0075))   ; 8 f7: 7.972..7.987
(body plate (shape (subtract (box 20 20 10) (cylinder-z (/ bore-d 2) 12))))
(body shaft (shape (cylinder-z (/ shaft-d 2) 20)))

; ---------------------------------------------------------------- press fit H7/s6
(param pin-bore-d 8 (tol 0 0.015))      ; H7: 8.000..8.015
(param pin-d 8.0275 (tol -0.0045 0.0045))     ; 8 s6: 8.023..8.032
(body pin-plate (shape (move (subtract (box 20 20 10) (cylinder-z (/ pin-bore-d 2) 12)) 30 0 0)))
(body pin (shape (move (cylinder-z (/ pin-d 2) 20) 30 0 0)))

; ---------------------------------------------------------------- loose spec (control)
(param loose-bore-d 8 (tol 0 0.1))
(param loose-shaft-d 8 (tol -0.05 0.05))
(body loose-plate (shape (move (subtract (box 20 20 10) (cylinder-z (/ loose-bore-d 2) 12)) 0 30 0)))
(body loose-shaft (shape (move (cylinder-z (/ loose-shaft-d 2) 20) 0 30 0)))

; The press fit's overlap is intended in a thin band at the pin's surface.
(assembly root
  (instance plate plate) (instance shaft shaft)
  (instance pin-plate pin-plate) (instance pin pin)
  (contact pin pin-plate (max (- (abs (- (len2 (- x 30) y) 4.01)) 0.03) (- (abs z) 10.1))))

(assert running-clear-5um (clearance (shape-of shaft) (shape-of plate) 0.005))
(assert running-gap-12um (value (- bore-d shaft-d) > 0.012))
(assert press-overlaps (value (- pin-d pin-bore-d) > 0))

; Controls: each fails inside the bands. The first holds at the middle
; sizes, and at the written values with --nominal.
; expect-status: running-clear-10um fails
(assert running-clear-10um (clearance (shape-of shaft) (shape-of plate) 0.010))
; expect-status: running-gap-20um fails
(assert running-gap-20um (value (- bore-d shaft-d) > 0.02))
; expect-status: running-overlaps fails
(assert running-overlaps (value (- shaft-d bore-d) > 0))
; expect-status: loose-apart fails
(assert loose-apart (disjoint (shape-of loose-shaft) (shape-of loose-plate)) (accuracy 0.01))

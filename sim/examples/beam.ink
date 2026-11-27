; A 10 mm square aluminium bar on two supports 200 mm apart, pressed down
; with 1000 N (about 100 kg) at the middle. Is it stiff enough, and strong
; enough? By hand (E 68900 MPa, yield 240 MPa, I = 10^4 / 12 = 833 mm^4):
;   deflection  P L^3 / 48 E I = 2.90 mm, and 0.02 of shear: 2.93 mm
;   stress      (P L / 4) c / I = 300 MPa: a safety factor of 240 / 300 = 0.80
; ink require examples/beam.ink: stiff-enough holds (2.94 mm), and
; strong-enough fails (0.79): the bar would bend for good.
; expect-status: strong-enough fails estimate
(import "std")

(body bar (shape (move (box 210 10 10) 100 0 0)))   ; x from -5 to 205, y and z from -5 to 5

; A support: the bar's bottom edge at x0, held within 1 mm of it, free to turn about it.
(fn edge (x0) (- (len2 (- x x0) (+ z 5)) 1))

(study press
  (shape (shape-of bar))
  (elastic (material aluminum-6061))
  (fix pin (edge 0) (pivot (at 0 0 -5) (axis 0 1 0)))
  (fix roller (edge 200) (dofs z) (pivot (at 200 0 -5) (axis 0 1 0)))   ; held up only: free to slide
  (load hand (force 0 0 -1000N) (region (- (abs (- x 100)) 1)))          ; spread over the middle 2 mm
  (res 0.5))                                                             ; cells of 0.5 mm

(require stiff-enough (max-displacement press < 4))       ; mm
(require strong-enough (min-safety-factor press > 1.5))

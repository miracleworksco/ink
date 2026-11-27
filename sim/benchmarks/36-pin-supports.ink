; Benchmark 36: pin, roller and bearing supports, and mechanisms.
;
; Task: support parts the way they are held, not clamped. (a) A steel beam
; 210 x 10 x 10 on its bottom edge at x = 0 (a pin) and x = 200 (a
; roller), 100 N down at midspan. (b) A steel shaft d20, z 0..200, on two
; self-aligning bearings at z = 10 and 190 (span 180), 1000 N radial at
; z = 100; the left bearing holds the shaft's spin and its axial place, the
; right one lets it slide. (c) An aluminium lever (x -10..100, 12 x 6) on a
; pin through a bore r3 at the origin, turned by a pure moment of 1000 N mm:
; free, then held by a stop at its tip, and on a frictionless seat.
; A plain fix holds all the material in its region, so a support written as
; a band would clamp the part; each support here pivots or slides where the
; real one does ((pivot ...), (slide ...)).
;
; Hand calculation (steel E 200 GPa, nu 0.3, G 76923):
;   (a) I = 10^4 / 12 = 833.33: P L^3 / 48 E I = 0.1000 mm, shear
;       P L / (4 kappa G A) = 0.00077 (kappa 0.8497): 0.10077 mm at
;       midspan; reactions 50 / 50 N, no axial reaction; the ends turn by
;       P L^2 / 16 E I = 1.500e-3 rad (the pin +, about +y; the roller -);
;       sigma = (P L / 4) c / I = 30.0 MPa.
;   (b) I = pi 20^4 / 64 = 7854: P L^3 / 48 E I = 0.07735, shear 0.00210
;       (kappa 0.8864): 0.07945 mm; reactions 500 / 500 N; the bearings
;       tilt by P L^2 / 16 E I = 1.289e-3 rad.
;   (c) on the pin alone the lever turns freely: there is no static
;       solution. A stop that turns about its line at x = 98.5 takes
;       1000 / 98.5 = 10.152 N, the pin the opposite.
;
; Pass criteria (the claims; ink solve, b36 in ink-ext/tests/cli/benchmarks.rs):
;   - (a) midspan deflection on the axis within 3% of 0.10077; reactions
;     within 0.1% of 50; end turns within 3% of 1.5e-3;
;   - (b) deflection at z = 100 within 3% of 0.07945, reactions within 0.1%
;     of 500;
;   - (c) the free lever and the frictionless seat fail as a mechanism,
;     naming the turn about (0, 0, 1) through the origin and the 1000 N mm;
;     with the stop, its reaction within 1e-4 of 10.152 N;
;   - the six deflection claims hold numerically.

(import "std")

; ---------------------------------------------------------------- (a) a beam on a pin and a roller
(body beam (shape (move (box 210 10 10) 100 0 0)) (attr material "steel"))
; The material within 1 mm of the beam's bottom edge line at x0 (along y):
; a pin of radius 1, rigid, turning about the line.
(fn edge (x0) (- (len2 (- x x0) (+ z 5)) 1))
(def midspan (- (abs (- x 100)) 1))
(study simply (shape (shape-of beam)) (elastic (material steel))
  (fix pin (edge 0) (pivot (at 0 0 -5) (axis 0 1 0)))
  (fix roller (edge 200) (dofs z) (pivot (at 200 0 -5) (axis 0 1 0)))
  (load mid (force 0 0 -100N) (region midspan))
  (res 1))
; Control: the same edges held rigid.
(study clamped-edges (shape (shape-of beam)) (elastic (material steel))
  (fix pin (edge 0))
  (fix roller (edge 200) (dofs z) (support))
  (load mid (force 0 0 -100N) (region midspan))
  (res 1))

; ---------------------------------------------------------------- (b) a shaft on self-aligning bearings
(body shaft (shape (move (cylinder-z 10 200) 0 0 100)) (attr material "steel"))
; A bearing holds the outer 2 mm of the seat over its 2 mm width.
(fn seat (zc) (max (- (abs (- z zc)) 1) (- 8 (len2 x y))))
(study bearings (shape (shape-of shaft)) (elastic (material steel))
  ; tilts about x and y through its centre; holds the spin and the axial place
  (fix left (seat 10) (pivot (at 0 0 10) (axis 1 0 0)) (pivot (at 0 0 10) (axis 0 1 0)))
  ; a ball joint (no axis: turns about every axis through the point), sliding axially
  (fix right (seat 190) (pivot (at 0 0 190)) (slide 0 0 1))
  (load radial (force -1000N 0 0) (region (- (abs (- z 100)) 1)))
  (res 1))
(study rigid-seats (shape (shape-of shaft)) (elastic (material steel))
  (fix left (seat 10))
  (fix right (seat 190) (dofs x y))
  (load radial (force -1000N 0 0) (region (- (abs (- z 100)) 1)))
  (res 1))

; ---------------------------------------------------------------- (c) a lever on a pin
(body lever (shape (subtract (move (box 110 12 6) 45 0 0) (cylinder-z 3 20))) (attr material "aluminum-6061"))
(def bore-band (- (len2 x y) 4.5))                    ; the bore's wall, r 3 to 4.5
(def middle (- (abs (- x 50)) 5))
; A mechanism: nothing holds the turn the pin leaves free. ink solve says so.
(study hinged (shape (shape-of lever)) (elastic (material aluminum-6061))
  (fix pin bore-band (pivot (at 0 0 0) (axis 0 0 1)) (rigid))   ; a steel pin holds the bore round
  (load turn (torque 0 0 1000N*mm) (region middle))
  (res 0.5))
; A stop at the tip, a roller turning about its line at x = 98.5. On cells of
; 0.5 its peak stress, where the load's region begins, rose 3.275, 3.340 and
; 3.419 MPa over three grids, which bound nothing; on 0.25 they converge.
(study hinged-stop (shape (shape-of lever)) (elastic (material aluminum-6061))
  (fix pin bore-band (pivot (at 0 0 0) (axis 0 0 1)) (rigid))   ; a steel pin holds the bore round
  (fix stop (- 97 x) (dofs y) (pivot (at 98.5 0 0) (axis 0 0 1)))
  (load turn (torque 0 0 1000N*mm) (region middle))
  (res 0.25))
; Frictionless in its bore: held radially at every point, so free to turn.
(study frictionless (shape (shape-of lever)) (elastic (material aluminum-6061))
  (fix seat (- (abs (- (len2 x y) 3.75)) 0.75) (dofs normal) (support))
  (load turn (torque 0 0 1000N*mm) (region middle))
  (res 0.5))

; ---------------------------------------------------------------- claims (5% of the hand values)
(require beam-stiff (max-displacement simply < 0.1058))
(require beam-not-too-stiff (max-displacement simply > 0.0957))
(require shaft-stiff (max-displacement bearings < 0.0834))
(require shaft-not-too-stiff (max-displacement bearings > 0.0755))
(require lever-held (max-displacement hinged-stop < 1))
(require lever-strong (max-stress hinged-stop < 30MPa))

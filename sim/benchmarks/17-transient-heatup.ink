; Benchmark 17: transient heat-up of a heated block vs lumped capacitance.
;
; Task: a 20 mm copper cube with a 5 W heater patch on its bottom face,
; cooled on all faces by air at 25 C, h = 25 W/(m^2 K), starting at 25 C.
; How hot is it after 10 minutes, and when does it pass 80 C?
;
; Hand calc (Biot = h a / k = 25 x 0.02 / 398 = 0.0013: lumped is exact
; to about 0.1%):
;   A = 6 x 400 mm^2 = 0.0024 m^2, hA = 0.06 W/K, steady rise 83.33 K
;   tau = rho c V / (h A) = 8960 x 385 x 8e-6 / 0.06 = 459.9 s
;   T(600 s) = 25 + 83.33 (1 - exp(-600 / 459.9)) = 85.72 C
;   T = 80 C at t = -tau ln(1 - 55 / 83.33) = 496.2 s
;
; Pass criteria (ink solve, transient block):
;   - final temperature (t = 600) within 2% of the 60.72 K rise;
;   - the history crosses 80 C within 3% of 496 s.

(import "std")

(def cube (box 20 20 20))
(body block (shape cube) (attr material "copper") (attr color "#C87533"))

(study heat-up (shape (shape-of block)) (heat (material copper))
  (source heater (power 5W) (region (max (- (len2 x y) 5) (- z -9))))
  (convection air (face cube) (coefficient 25W/m2/K) (ambient 25degC))
  (transient (duration 600s) (steps 60) (initial 25degC))
  (res 1))

(require below-limit (max-temperature heat-up < 90degC))
; T(600 s) within 2% of the lumped 60.72 K rise (85.72 C, above): the
; block only heats, so its peak over time is the last step's
(require t600-lo (max-temperature heat-up > 84.506degC))
(require t600-hi (max-temperature heat-up < 86.934degC))

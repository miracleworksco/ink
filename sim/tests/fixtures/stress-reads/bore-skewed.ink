; An 80 x 40 x 4 plate with a 1 mm bore, pulled by 160 N: CalculiX 2.956 MPa
; at the bore. The bore is written (1 - r)(1 + 0.95 sin 25x), the same sign as
; 1 - r and so the same part, its field far from a distance's near the bore.
; True: lt-3-0 and above, gt-2-8, sf-gt-80. False: lt-2-8, gt-3-1.
(import "std")
(def part (max (box 80 40 4) (* (- 1 (sqrt (+ (square x) (square y)))) (+ 1 (* 0.95 (sin (* 25 x)))))))
(study s (shape part) (elastic (E 200000MPa) (nu 0.3) (yield 250MPa))
 (fix wall (+ x 39.9)) (load pull (force 160N 0 0) (region (- 39.9 x))) (res 1))
(require lt-2-8 (max-stress s < 2.8MPa))
(require lt-3-0 (max-stress s < 3.0MPa))
(require lt-3-05 (max-stress s < 3.05MPa))
(require lt-3-1 (max-stress s < 3.1MPa))
(require lt-3-2 (max-stress s < 3.2MPa))
(require lt-3-4 (max-stress s < 3.4MPa))
(require gt-2-8 (max-stress s > 2.8MPa))
(require gt-3-1 (max-stress s > 3.1MPa))
(require sf-gt-80 (min-safety-factor s > 80))

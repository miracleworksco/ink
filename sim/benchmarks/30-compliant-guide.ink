; Benchmark 30: a printable PETG two-leaf parallel-guide stage.
;
; Task: translate a shuttle 0.5 mm in y, with a generous yield margin and a
; first frequency above 85% of the hand value. Two leaves run x = 0..40,
; centres y = +-9, 1 mm thick in y and 6 mm wide in z. Base and shuttle are
; 6 x 24 x 6 mm. Each leaf has r0.75 circular root fillets, written as planar
; arcs and extruded. Printed flat in xy; the material is the library's
; isotropic in-plane PETG (E 1500 MPa, nu 0.38, density 1270 kg/m^3, yield
; 47 MPa). An elastic one-cycle check, not a fatigue or creep qualification.
;
; Hand calculation (rigid shuttle, two ideal fixed-guided Euler-Bernoulli
; leaves; root fillets, shear and shuttle flexibility ignored):
;   - I = b h^3 / 12 = 0.5 mm^4; k = 2 x 12 E I / L^3 = 0.28125 N/mm.
;   - F = k x 0.5 = 0.140625 N; compliance F delta = 0.0703125 N mm.
;   - root stress 3 E h delta / L^2 = 1.40625 MPa (yield factor 33.4). Under
;     the force, each leaf's moment is (F/2)(L/2)(1 - 2x/L): at x = 4 the
;     surface stress is 0.8 x 1.40625 = 1.125 MPa, whatever the stiffness.
;   - shuttle 1.09728 g, each leaf 0.3048 g; the guided shape
;     u/delta = 3 s^2 - 2 s^3 gives integral (u/delta)^2 ds = 13/35, so
;     m_eff = 1.09728 + 2 (13/35) 0.3048 = 1.323703 g;
;     f1 = sqrt(1e6 k / m_eff) / 2 pi = 73.362 Hz.
;   - The leaves rotate at most 1.5 delta / L = 0.019 rad (1.1 deg): linear
;     kinematics apply.
;
; Pass criteria (the claims; ink solve, b30 in ink-ext/tests/cli/benchmarks.rs):
;   - check clean; assert: six claims hold numerically, the hand-stress
;     claim proven. Thresholds are expressions of params (hand values times
;     a margin), so they follow --param stroke=... and force-scale=...
;   - solve --study translate: F / uy within 15% of k; surface stress at
;     (4, 9.5, 0) within 10% of 1.125 MPa; reported peaks below yield;
;     equilibrium imbalance < 1e-5; no warnings (no large rotation).
;   - solve --study ring: first frequency within 10% of 73.362 Hz.

(import "std")
(param L 40)
(param h 1)
(param b 6)
(param separation 18)
(param stroke 0.5)
(param force-scale 1)
(param root-r 0.75)
(def block-len 6)
(def block-width 24)
(def base (move (box block-len block-width b) (* -0.5 block-len) 0 0))
(def shuttle (move (box block-len block-width b) (+ L (* 0.5 block-len)) 0 0))
; Each leaf has four circular concave root arcs in the xy section.
; For x in [0,r], extra half-width r-sqrt(r^2-(r-x)^2); mirror at L.
(fn radiused-leaf (yc)
  (let ((q (max 0 (- root-r x) (- x (- L root-r))))
        (extra (- root-r (sqrt (- (square root-r) (square q))))))
    (max (neg x) (- x L) (- (abs (- y yc)) (+ (* 0.5 h) extra))
         (- (abs z) (* 0.5 b)))))
(def upper (radiused-leaf (* 0.5 separation)))
(def lower (radiused-leaf (* -0.5 separation)))
(def guide (min base shuttle upper lower))
(body guide (shape guide) (attr material "petg") (attr density 1.27g/cm3))
(def E (get (material petg) "E"))
(def petg-yield (get (material petg) "yield"))
(def I (/ (* b (pow h 3)) 12))
(def k-hand (/ (* 24 E I) (pow L 3)))
(def force-hand (* k-hand stroke))
(def stress-hand (/ (* 3 E h stroke) (square L)))
(def m-eff (* (get (material petg) "density")
  (+ (* block-len block-width b) (* 2 (/ 13 35) L b h)))) ; grams
(def freq-hand (/ (sqrt (/ (* 1e6 k-hand) m-eff)) (* 2 pi)))
; The applied force is force-scale times the hand force: hand travel and
; compliance at that force.
(def travel-hand (* force-scale stroke))
(def compliance-hand (* force-scale force-hand travel-hand))
(study translate (shape (shape-of guide)) (elastic (material petg))
  (fix clamp base)
  (load actuator (force 0 (* 1N (* force-scale force-hand)) 0) (region shuttle))
  (res 0.25))
(study ring (shape (shape-of guide)) (elastic (material petg))
  (fix clamp base) (modes 1) (res 0.25))
(require below-yield (max-stress translate < (* 1MPa petg-yield)))
(require safe (min-safety-factor translate > 5))
(require travel-upper (max-displacement translate < (* 1.15 travel-hand)))
(require compliance-lower (compliance translate > (* 1N*mm (* 0.8 compliance-hand))))
(require compliance-upper (compliance translate < (* 1N*mm (* 1.15 compliance-hand))))
(require first-mode (min-frequency ring > (* 1Hz (* 0.85 freq-hand))))
; ... and within 10% of it (73.362 Hz, above)
(require first-mode-lo (min-frequency ring > (* 1Hz (* 0.9 freq-hand))))
(require first-mode-hi (min-frequency ring < (* 1Hz (* 1.1 freq-hand))))
(require hand-stress-safe (value stress-hand < (/ petg-yield 5)))

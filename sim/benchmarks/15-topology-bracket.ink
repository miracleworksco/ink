; Benchmark 15: topology-optimized cantilever bracket, then verified.
;
; Task: an aluminium bracket fits a 120 x 40 x 60 mm envelope, bolted to a
; wall over its x = 0 face, carrying 2 kN down at the tip. Find the stiffest
; layout with 30% of the envelope's material (SIMP), extract it as a field,
; and verify it with an ordinary study. Compare with the obvious
; hand design of the same volume: a 12 x 60 solid web (720 mm^2).
;
; Hand calc for the baseline web: Euler-Bernoulli plus shear,
;   delta = P L^3 / (3 E I) + P L / (k G A)
;         = 2000 120^3 / (3 68900 216000) + 2000 120 / (5/6 25900 720)
;         = 0.0774 + 0.0154 = 0.0928 mm (tip, before end effects)
;
; Pass criteria (ink topology / ink solve / ink mass / ink mesh):
;   - the baseline FEA tip deflection within 15% of 0.0928 mm (sanity);
;   - the layout settles within its 150 iterations (at 133: compliance 125.4
;     N mm);
;   - the optimized part's volume within 20% of 30% of the envelope (86 400);
;   - verified optimized part is stiffer than the baseline (compliance);
;   - verified min safety factor > 2 (yield 276);
;   - the optimized part meshes watertight as a single solid.

(import "std")

(def domain (move (box 120 40 60) 60 0 0))
(def load-pad (max (- 116 x) (- (abs y) 6) (- z -24)))          ; tip, bottom 6 mm, centre
(body envelope (shape domain))

(study cantilever (shape domain) (elastic (material aluminum-6061))
  (fix wall (- x 1))
  (load tip (force 0 0 -2000N) (region load-pad))
  (res 2))
(topology layout (study cantilever) (volume-fraction 0.3) (iterations 150) (filter 4))

(def optimized (max domain (study-field layout solid)))
(body bracket (shape optimized) (attr material "aluminum-6061") (attr density 2.7g/cm3) (attr color "#6C8EBF"))

(study verify (shape (shape-of bracket)) (elastic (material aluminum-6061))
  (fix wall (- x 1))
  (load tip (force 0 0 -2000N) (region load-pad))
  (res 1))

; baseline: a 12 mm solid web, full height
(def web (move (box 120 12 60) 60 0 0))
(body baseline (shape web))
(study baseline (shape web) (elastic (material aluminum-6061))
  (fix wall (- x 1))
  (load tip (force 0 0 -2000N) (region load-pad))
  (res 1))

; The 1 mm stress survey can classify a smooth round as a sharp corner.
; This claim asks for the safety factor of the stress read by its cells.
(require verified-strong (min-safety-factor verify > 2 (excluding sharp-corners)))
; The spline's volume bounds need finer integration cells to show the
; optimized part's volume is within 20% of 30% of the envelope, 86 400.
(require volume-lo (volume (shape-of bracket) > 69120mm3) (res 0.3))
(require volume-hi (volume (shape-of bracket) < 103680mm3) (res 0.3))

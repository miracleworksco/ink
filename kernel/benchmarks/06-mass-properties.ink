; Benchmark 06: mass properties (volume, centre of mass, inertia tensor).
;
; Task: the mass, centre of mass and inertia tensor of the polygon plate
; (benchmark 04) and the swept pipe (01). Each body carries the ten raw
; moments as field attributes (rho, rho x, ..., rho y z); the benchmark's test
; (b06 in crates/ink-cli/tests/cli/benchmarks.rs) weighs each by its m
; (ink mass --density m) for
;   m = int rho,  c = int rho p / m,
;   I = int rho (|p|^2 E - p p^T) dV  - m (|c|^2 E - c c^T)   (about the COM)
; then compares the plate with its exact (analytic polygon: shoelace) moments.
;
; Pass criteria (ink mass --res 0.25), for the plate: mass within
; 0.1% (relative), centre of mass within 0.001 mm, inertia within 0.1% of the
; largest diagonal term.

(import "std")
(import "04-polygon-plate.ink")
(import "01-swept-pipe.ink")

(def rho-plate 0.0027)          ; g/mm^3 (aluminium)
(def rho-pipe 0.00896)          ; g/mm^3 (copper)

(body plate-mp
  (shape (shape-of plate))
  (attr m rho-plate)
  (attr mx (* rho-plate x)) (attr my (* rho-plate y)) (attr mz (* rho-plate z))
  (attr mxx (* rho-plate x x)) (attr myy (* rho-plate y y)) (attr mzz (* rho-plate z z))
  (attr mxy (* rho-plate x y)) (attr mxz (* rho-plate x z)) (attr myz (* rho-plate y z)))

(body pipe-mp
  (shape (shape-of swept-pipe))
  (attr m rho-pipe)
  (attr mx (* rho-pipe x)) (attr my (* rho-pipe y)) (attr mz (* rho-pipe z))
  (attr mxx (* rho-pipe x x)) (attr myy (* rho-pipe y y)) (attr mzz (* rho-pipe z z))
  (attr mxy (* rho-pipe x y)) (attr mxz (* rho-pipe x z)) (attr myz (* rho-pipe y z)))

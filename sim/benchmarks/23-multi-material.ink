; Benchmark 23: multi-material assembly with by-body properties.
;
; Task (a): a bonded bar, 10 x 10 section: aluminium 6061 over x 0..50 and
; steel over x 50..100, clamped at x = 0 and pulled with 10 kN at x = 100.
; Mass properties of the pair come from each body's density attribute.
;
; Hand calc (a): springs in series,
;   delta = F / A (L1 / E1 + L2 / E2) = 100 (50/68900 + 50/200000) = 0.09757 mm
;   mass 5000 (0.0027 + 0.00785) = 52.75 g, COM x = (13.5 25 + 39.25 75) / 52.75
;   = 62.204 mm
;
; Pass criteria (ink solve, ink mass):
;   - (a) tip axial displacement within 3% of 0.09757 mm;
;   - (a) mass within 0.3% of 52.75 g, COM x within 0.1 mm of 62.204.

(import "std")

(body al-part (shape (move (box 50 10 10) 25 0 0)) (attr material "aluminum-6061") (attr density 2.7g/cm3) (attr color "#B8C2CC"))
(body st-part (shape (move (box 50 10 10) 75 0 0)) (attr material "steel") (attr density 7.85g/cm3) (attr color "#606060"))
(assembly bar (instance al al-part) (instance st st-part))
; the hand calc's mass 52.75 g within 0.3%, COM x 62.204 within 0.1 mm (proven)
(require mass-lo (integral bar (attr density) > 52.592) (res 0.1))
(require mass-hi (integral bar (attr density) < 52.908) (res 0.1))
(require com-x (mean bar (vec (- x 62.304) (- 62.104 x)) (weight (attr density)) < 0) (res 0.1))

(study pull (shape (shape-of bar))
  (elastic (E (* 1MPa (by-body bar E))) (nu (by-body bar nu)))
  (fix wall (- x 0.5))
  (load pull (force 10000N 0 0) (region (- 99.5 x)))
  (res 1))

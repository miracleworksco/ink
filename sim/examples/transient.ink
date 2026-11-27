; Heat over time: loads are expressions of t. A PETG bar heated at its
; tip by a 1 W source for 60 s, then left to conduct, its wall end held at
; 20 C. Its tip passes PETG's glass transition: the bar's 1-D conduction
; (uniform section, heated through the tip's 5 mm) peaks at 87.4 C as the
; heat stops; ink reads about 86.4 C, its steps' error about 0.9. On 60
; steps its three runs read 85.8, 84.7 and 83.4: not yet converging as
; their order says (steps of 40 s are long against the 60 s pulse), so a
; claim on them is not settled.
; expect-status: stays-below-glass fails estimate
(import "std")

(body bar (shape (move (box 100 10 10) 50 0 0)))
(study heat-pulse (shape (shape-of bar)) (heat (k 0.2W/m/K) (density 1270kg/m3) (capacity 1300J/kg/K))
  (temperature wall (- x 0.5) 20degC)
  (source torch (power (* 1W (select (- t 60s) 1 0))) (region (- 95 x)))
  (transient (duration 600s) (steps 120) (initial 20degC))     ; initial temperature, C
  (res 2.5))

(require stays-below-glass (max-temperature heat-pulse < 80degC))

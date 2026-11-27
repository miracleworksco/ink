; nafems-t3
; Source for inputs and printed targets: https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-t3.htm
; Temperature atx.08 m, t 32 s; target 36.60 C.
; Bar length.1 m, square.01 m, k 35, rho 7200, cp 440.5;initial 0 C,
; left 0 C, right 100*sin(pi*t/40)C; lateral adiabatic. Exact 1D thermal
; problem has identical finite 3D solution. The changing prescribed
; temperature must be represented exactly; parser refusal is a result.
; scoreboard: {"id":"nafems-t3","physics":"thermal transient","target":36.6,"target_kind":"converged reference","units":"C","study":"s","qualified":true,"sources":["https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-t3.htm"],"quantity":"Temperature atx.08m,t32s; target36.60C.","probe":"temperature","points":[[80,0,0]]}
(import "std")
(def part (move (box 100 10 10) 50 0 0))
(body part (shape part))
(study s (shape part)
  (heat (k 35W/m/K) (density 7200kg/m3) (capacity 440.5J/kg/K))
  (temperature left (- x 0.000001) 0)
  (temperature right (- 99.999999 x) (+ 0degC (* 100 (sin (* 4.5deg/s t)))))
  (transient (duration 32s) (initial 0degC)))

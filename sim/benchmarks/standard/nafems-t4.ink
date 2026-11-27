; nafems-t4
; Source for inputs and printed targets: https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-t4.htm
; Temperature atE=(.6, .2)m; target 18.3 C.
; Rectangle.6mx 1 m, thickness 1 m per public 3D reproduction; k 52, h 750.
; y 0 fixed 100 C, x 0 insulated, x.6 and y 1 convect to ambient 0 C. Extruded
; ends insulated make planar solution exact for this finite 3D solid.
; Exact separated solution is 18.253756537134212 C, consistent with target
; rounded to 0.1 C. Analytic qualification is retained with the report.
; scoreboard: {"id":"nafems-t4","physics":"thermal steady","target":18.3,"target_kind":"converged reference","units":"C","study":"s","qualified":true,"sources":["https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-t4.htm"],"quantity":"Temperature atE=(.6,.2)m; target18.3C.","qualification":"Separated exact solution18.253756537134212C rounds to published18.3C; same finite3D PDE","probe":"temperature","points":[[600,200,0]],"target_rounding":0.05}
(import "std")
(def part (move (box 600 1000 1000) 300 500 0))
(body part (shape part))
(study s (shape part)
  (heat (k 52W/m/K))
  (temperature bottom (- y 0.000001) 100degC)
  (convection right (- 599.999999 x) (coefficient 750W/m2/K) (ambient 0))
  (convection top (- 999.999999 y) (coefficient 750W/m2/K) (ambient 0)))

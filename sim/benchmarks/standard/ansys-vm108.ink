; ansys-vm108
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM108.html
; Temperature at x 10 ft, y 0; target 40 F. Radius 20 ft, length 5 ft, k 1 BTU/(h ft F).
; Lateral temperature 80 F*cos(theta), insulated ends. Exact 3D T(F)=80*x/R.
; Absolute conversion is (T(F)-32)*5/9; multiplying an absolute Celsius
; temperature by cos(theta) would be wrong. Field-valued temperature
; boundaries are intentionally submitted as specified; current parser may refuse.
; scoreboard: {"id":"ansys-vm108","physics":"thermal steady","target":4.444444444444445,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM108.html"],"quantity":"Temperature at x10ft,y0; target40F. Radius20ft,length5ft,k1BTU/(h ft F).","probe":"temperature","points":[[3048,0,0]]}
(import "std")
(def part (cylinder-z 20ft 5ft))
(body part (shape part))
(study s (shape part)
  (heat (k 1.7307346663713912W/m/K))
  (temperature side (- 6095.999999 (len2 x y)) (+ 0degC (* (/ 5 9) (- (* 80 (/ x (len2 x y))) 32)))))

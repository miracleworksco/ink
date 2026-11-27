; ansys-vm58
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM58.html
; Centerline maximum temperature; printed 419.9 F, rounded to 0.1 F.
; Radius .375 in, length 1 ft, insulated ends, k=13 BTU/(h ft F),
; h=5 BTU/(h ft^2 F), ambient 70 F, generation 111311.7 BTU/(h ft^3).
; Exact finite-solid solution T=Ta+q R/(2h)+q(R^2-r^2)/(4k).
; Conversions use IT BTU=1055.05585262 J, ft=.3048 m, F difference=5/9 K.
; scoreboard: {"id":"ansys-vm58","physics":"thermal steady","target":215.5,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM58.html"],"quantity":"Centerline maximum temperature; printed 419.9 F, rounded to 0.1 F.","published_target":419.9,"published_units":"F","target_rounding":0.02777777777777778,"metric":"max-temperature"}
(import "std")
(def part (cylinder-z 9.524999999999999 1ft))
(body part (shape part))
(study s (shape part)
  (heat (k 22.499550662828085W/m/K))
  (convection side (- 9.524999 (len2 x y)) (coefficient 28.39131670556744W/m2/K) (ambient 70degF))
  (source generated (power 100.08377589552867W) (region -1)))
; This wide measurement envelope requests the reported grid band; it is not a target test.
(require measurement (max-temperature s < 1000000degC))

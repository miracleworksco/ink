; ansys-vm92
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM92.html
; Inner surface maximum; printed rounded target 2957 F.
; Wall 9 in k=.8 and 5 in k=.1 BTU/(h ft F), area 1 ft^2; sides adiabatic.
; Inner h 12@3000 F, outer h 2@80 F, h in BTU/(h ft^2 F).
; Exact 3D series resistance: q/A=(3000-80)/(1/12+.75/.8+(5/12)/.1+1/2).
; Maximum Tf-q/(A hi). Same IT BTU/ft/F conversion as VM 58.
; scoreboard: {"id":"ansys-vm92","physics":"thermal steady","target":1625.0,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM92.html"],"quantity":"Inner surface maximum; printed rounded target2957 F.","target_rounding":0.2777777777777778,"metric":"max-temperature"}
(import "std")
(def part (move (box 14in 1ft 1ft) 7in 0 0))
(body part (shape part))
(study s (shape part)
  (heat (k (* 1W/m/K (select (- x 9in) 1.384587733097113 0.17307346663713913))))
  (convection inside (- x 0.000001) (coefficient 68.13916009336185W/m2/K) (ambient 3000degF))
  (convection outside (- 355.599999 x) (coefficient 11.356526682226976W/m2/K) (ambient 80degF)))
; This wide measurement envelope requests the reported grid band; it is not a target test.
(require measurement (max-temperature s < 1000000degC))

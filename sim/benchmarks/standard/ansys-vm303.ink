; ansys-vm303
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM303.html
; Maximum temperature on conducting bar; exact 128 C atx 6 mm.
; L 10 mm, area.0004 mm^2, k 1 W/(mm K)=1000 W/(m K), Q 6 W/mm^3, ends 20/80 C.
; The surrounding block has zero conductivity, so removing it gives exactly
; the same bar heat equation and adiabatic sides. Choose a square.02 mm
; section consistent with stated area; section shape cannot affect this 1D
; solution. T(x)=20+6x+3x(10-x), x in mm; total power=.024 W.
; scoreboard: {"id":"ansys-vm303","physics":"thermal steady","target":128,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM303.html"],"quantity":"Maximum temperature on conducting bar; exact128C atx6mm.","metric":"max-temperature"}
(import "std")
(def part (move (box 10 0.02 0.02) 5 0 0))
(body part (shape part))
(study s (shape part)
  (heat (k 1000W/m/K))
  (temperature left (- x 0.000001) 20degC)
  (temperature right (- 9.999999 x) 80degC)
  (source generated (power 0.024W) (region -1)))
; This wide measurement envelope requests the reported grid band; it is not a target test.
(require measurement (max-temperature s < 1000000degC))

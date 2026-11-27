; ansys-vm280
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM280.html
; Maximum temperature at inner radius.06 m; printed 222.300 C.
; Radii .06, .0675, .1175, .1675 m; k 49, .15, .48 W/(m K).
; Inner h 85@230 C, outer h 18@35 C. Axial length not specified: use 1 m,
; ends adiabatic, exact radial 3D solution independent of length.
; Thermal resistances are logarithmic layer sums plus both convection films.
; scoreboard: {"id":"ansys-vm280","physics":"thermal steady","target":222.3,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM280.html"],"quantity":"Maximum temperature at inner radius.06m; printed222.300C.","target_rounding":0.0005,"metric":"max-temperature"}
(import "std")
(def part (extrude-z (ring-2d 60 167.5) 0 1000))
(body part (shape part))
(study s (shape part)
  (heat (k (* 1W/m/K (select (- (len2 x y) 67.5) 49 (select (- (len2 x y) 117.5) 0.15 0.48)))))
  (convection inside (- (len2 x y) 60.000001) (coefficient 85W/m2/K) (ambient 230degC))
  (convection outside (- 167.499999 (len2 x y)) (coefficient 18W/m2/K) (ambient 35degC)))
; This wide measurement envelope requests the reported grid band; it is not a target test.
(require measurement (max-temperature s < 1000000degC))

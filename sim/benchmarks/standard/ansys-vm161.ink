; ansys-vm161
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM161.html
; Heat loss through outer convection per 1 ft; published 362.0 BTU/h.
; Radii 1.535, 1.75, 2.25 in; length 1 ft with insulated ends; k 25/.11 BTU/(h ft F).
; Inner h 40@300 F, outer h 4@80 F in BTU/(h ft^2 F).
; Exact 3D resistance 1/(hi*2pi*ri*L)+ln(r 2/r 1)/(2pi*k 1*L)
; +ln(ro/r 2)/(2pi*k 2*L)+1/(ho*2pi*ro*L). Report heat in W.
; scoreboard: {"id":"ansys-vm161","physics":"thermal steady","target":106.09172740234443,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM161.html"],"quantity":"Heat loss through outer convection per1ft; published362.0BTU/h.","units_override":"W","target_rounding":0.014653553508611113,"probe_scale":-1,"probe":"reaction:outside","points":[]}
(import "std")
(def part (extrude-z (ring-2d 1.535in 2.25in) 0 1ft))
(body part (shape part))
(study s (shape part)
  (heat (k (* 1W/m/K (select (- (len2 x y) 1.75in) 43.26836665928478 0.19038081330085302))))
  (convection inside (- (len2 x y) 38.989001) (coefficient 227.13053364453953W/m2/K) (ambient 300degF))
  (convection outside (- 57.149999 (len2 x y)) (coefficient 22.713053364453952W/m2/K) (ambient 80degF)))

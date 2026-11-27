; ansys-vm15-uniform-clamped
; Input/target source: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM15.html
; Signed center-bottom z displacement: target -.08736 in under downward load.
; Disk radius 40 in, thickness 1 in, E 30e6 psi, nu.3. Clamp cylindrical
; rim; pressure 6 psi on the top, bottom free. Asked point is bottom center.
; Kirchhoff plate target w=-p*a^4/(64D). Mindlin transverse shear adds
; magnitude p*a^2/(4*kappa*G*t), about+.286%, with kappa 5/6.
; Exact 3D rim and through-thickness corrections are not independently
; bounded; no pass or miss is counted against the plate target.
; scoreboard: {"id":"ansys-vm15-uniform-clamped","physics":"static","target":-2.2189439999999996,"target_kind":"idealisation","units":"mm","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM15.html"],"quantity":"Signed center-bottom z displacement: target -.08736 in under downward load.","target_signed":true,"probe":"displacement_mm.2","points":[[0,0,-12.7]]}
(import "std")
(def part (cylinder-z 40in 1in))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0.3))
  (fix rim (- 1015.999999 (len2 x y)))
  (load surface (pressure 0.041368543759008MPa) (region (- 12.699999 z))))

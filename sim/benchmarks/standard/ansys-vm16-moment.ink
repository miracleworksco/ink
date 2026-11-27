; ansys-vm16-moment
; Input/target source: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM16.html
; Signed tip-center y displacement: published +.005 in = +.127 mm.
; Beam length 10 in, height 2 in, unit thickness 1 in; E 30e6 psi, nu 0.
; Root clamped; end moment+2000 in lbf about z. A vanishingly thin end
; region distributes a pure linear axial traction, not a shear force.
; For nu 0 the exact 3D displacement ux=-k*x*y, uy=k*x^2/2, uz=0,
; k=M/(EI), satisfies the free sides and rigid root. Thus the center
; vertical displacement has no beam-to-solid correction in this case.
; The separate surface stress target at 1 in from the root is not substituted.
; scoreboard: {"id":"ansys-vm16-moment","physics":"static","target":0.127,"target_kind":"exact","units":"mm","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM16.html"],"quantity":"Signed tip-center y displacement: published +.005 in = +.127 mm.","probe":"displacement_mm.1","points":[[254,0,0]]}
(import "std")
(def part (move (box 10in 2in 1in) 5in 0 0))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0))
  (fix root (- x 0.000001))
  (load end (torque 0 0 225969.65805523336N*mm) (region (- 253.999999 x))))

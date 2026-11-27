; ansys-vm16-shear
; Input/target source: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM16.html
; Signed tip-center y displacement: published +.005 in = +.127 mm.
; Beam length 10 in, height 2 in, unit thickness 1 in; E 30e6 psi, nu 0.
; Root clamped; end shear 300 lbf in+y. Euler-Bernoulli target omits shear.
; Timoshenko adds F*L/(kappa*G*A), kappa 5/6, G=E/2: +.00012 in,
; +2.4% of the printed bending target. The rigid root and end traction
; introduce 3D boundary layers whose net correction is not bounded here.
; This idealisation is ungraded, including any apparent match.
; scoreboard: {"id":"ansys-vm16-shear","physics":"static","target":0.127,"target_kind":"idealisation","units":"mm","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM16.html"],"quantity":"Signed tip-center y displacement: published +.005 in = +.127 mm.","probe":"displacement_mm.1","points":[[254,0,0]]}
(import "std")
(def part (move (box 10in 2in 1in) 5in 0 0))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0))
  (fix root (- x 0.000001))
  (load end (force 0 1334.46648457815N 0) (region (- 253.999999 x))))

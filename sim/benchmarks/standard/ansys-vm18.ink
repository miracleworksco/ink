; ansys-vm18
; Input/target source: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM18.html
; Signed tip-center z displacement: published -2.648 in.
; Quarter-circle centerline radius 100 in, solid circular diameter 2 in;
; E 30e6 psi, nu.3. Clamp at(x 100 in, y 0), load 50 lbf downward at
; (x 0, y 100 in). Geometry is a true torus quadrant, not straight segments.
; The source beam target includes bending and torsion (Timoshenko,
; Strength of Materials I, equation 241). Transverse shear adds roughly
; -F*pi*R/(2*kappa*G*A)=-.0000733 in, kappa=6(1+nu)/(7+6nu).
; Curved-section and rigid-root corrections are not bounded; ungraded.
; scoreboard: {"id":"ansys-vm18","physics":"static","target":-67.25919999999999,"target_kind":"idealisation","units":"mm","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM18.html"],"quantity":"Signed tip-center z displacement: published -2.648 in.","qualification":"Curved-beam target; shear/section/clamp solid correction not bounded","probe":"displacement_mm.2","points":[[0,2540,0]]}
(import "std")
(def part (max (torus-z 100in 1in) (neg x) (neg y)))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0.3))
  (fix root (- y 0.000001))
  (load tip (force 0 0 -222.41108076302498N) (region (- x 0.000001))))

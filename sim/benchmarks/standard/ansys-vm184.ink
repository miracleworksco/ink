; ansys-vm184
; Input/target source: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM184.html
; Tip-center z displacement under 1 lbf weak-axis shear: target+.432 in.
; Cantilever 6 in long, section .2 in along y and .1 in along z;
; E 10e6 psi, nu.3. Root clamped; end load in+z.
; The source also has axial and other-direction load cases; they are not
; counted as extra problems. Shear adds F*L/(kappa*G*A)=.0000936 in,
; about+.022%, with kappa 5/6 and G=E/[2(1+nu)]. The clamp's 3D
; correction is not bounded; the beam comparison remains ungraded.
; scoreboard: {"id":"ansys-vm184","physics":"static","target":10.9728,"target_kind":"idealisation","units":"mm","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM184.html"],"quantity":"Tip-center z displacement under1 lbf weak-axis shear: target+.432 in.","qualification":"Beam shear correction given; finite-solid clamp correction unqualified","probe":"displacement_mm.2","points":[[152.4,0,0]]}
(import "std")
(def part (move (box 6in 0.2in 0.1in) 3in 0 0))
(body part (shape part))
(study s (shape part)
  (elastic (E 68947.57293168MPa) (nu 0.3))
  (fix root (- x 0.000001))
  (load tip (force 0 0 4.4482216152605N) (region (- 152.399999 x))))

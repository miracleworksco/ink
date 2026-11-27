; ansys-vm10
; Source for all dimensions, loads, material and printed target: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM10.html
; Compressive extreme-fiber axial stress at midspan top: target-700 psi.
; T beam: length 100 in (source arbitrary length), flange 9 in wide by 4 in
; high at bottom, stem 1.5 in wide by 16 in high above it; total height 20 in.
; E 30e6 psi, end bending moment 100000 in lbf; root clamped.
; Source omits nu: retain.3 and leave true 3D comparison ungraded.
; Section centroid is 6 in above bottom, I=2000 in^4; sigma=-M*(y-6)/I
; gives-700 psi at top and+300 psi at bottom. Pure bending theory away
; from ends is exact Saint-Venant stress; clamp corrections decay with
; distance/section width, but amplitude is not independently bounded.
; The end/clamp singular maxima are not substituted for beam-fiber stress.
; scoreboard: {"id":"ansys-vm10","physics":"static","target":-4.8263301052175995,"target_kind":"idealisation","units":"MPa","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM10.html"],"quantity":"Compressive extreme-fiber axial stress at midspan top: target-700psi.","units_override":"MPa","qualification":"Beam fiber target; nu unspecified and finite solid clamp decay amplitude unbounded","probe":"stress_mpa.xx","points":[[1270,508,0]]}
(import "std")
(def part (min (move (box 100in 4in 9in) 50in 2in 0) (move (box 100in 16in 1.5in) 50in 12in 0)))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0.3))
  (fix root (- x 0.000001))
  (load end (torque 0 0 11298482.902761668N*mm) (region (- 2539.999999 x))))

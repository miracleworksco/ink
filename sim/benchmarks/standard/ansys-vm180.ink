; ansys-vm180
; Source for all dimensions, loads, material and printed target: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM180.html
; Inner-fiber tangential stress at 45 degrees: tensile target+655 psi.
; Quarter annulus radii 3.5/4.5 in, axial thickness 1 in, E 30e6 psi, nu 0.
; Root at y 0 held; end at x 0 carries pure moment-100 in lbf about z.
; Source uses a load-spreading end beam and ignores its end nodes.
; Midarc is read here to avoid end singularities. Curved-beam theory
; r_n=(r_o-r_i)/ln(r_o/r_i), e=R-r_n gives
; sigma=M/(A*e)*(r_n/r-1), close to 655 psi. The exact 2D Airy
; solution and clamp/end layers can shift it; with nu 0 the corresponding
; plane-elastic field can extrude into 3D, but those boundary effects
; are not bounded here. This target remains ungraded.
; scoreboard: {"id":"ansys-vm180","physics":"static","target":4.51606602702504,"target_kind":"idealisation","units":"MPa","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM180.html"],"quantity":"Inner-fiber tangential stress at45degrees: tensile target+655psi.","tensor_direction":[-0.7071067811865475,0.7071067811865475,0],"units_override":"MPa","qualification":"Curved-beam target; exact Airy and finite end/clamp correction unbounded","probe":"stress_mpa","points":[[62.86179284748406,62.86179284748406,0]]}
(import "std")
(def part (max (extrude-z (ring-2d 3.5in 4.5in) -0.5in 0.5in) (neg x) (neg y)))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0))
  (fix root (- y 0.000001))
  (load end (torque 0 0 -11298.482902761669N*mm) (region (- x 0.000001))))

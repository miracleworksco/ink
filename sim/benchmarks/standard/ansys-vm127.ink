; ansys-vm127
; Source for all dimensions, loads, material and printed target: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM127.html
; First buckling load, target 38.553 lbf from hinged Euler column.
; Square .5 in, length 200 in, E 30e6 psi, 1 lbf compression. Source omits nu;
; use.3 and keep target ungraded. Rigid pin end sections can turn about
; both bending axes; loaded end slides axially. Whole solid retains both
; bending directions. Shear-flexible estimate P=P_E/(1+P_E/(kGA)),
; k 5/6, gives lower than Euler; finite solid pin/Poisson correction unknown.
; scoreboard: {"id":"ansys-vm127","physics":"buckling","target":38.553,"target_kind":"idealisation","units":"lbf (factor on 1lbf)","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM127.html"],"quantity":"First buckling load, target38.553lbf from hinged Euler column.","qualification":"Euler target; nu unspecified and finite solid pin correction unqualified","target_factor_to_load":4.4482216152605,"metric":"min-buckling-factor"}
(import "std")
(def part (move (box 200in 0.5in 0.5in) 100in 0 0))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0.3))
  (fix pin (- x 0.000001) (pivot (at 0 0 0) (axis 0 1 0)) (pivot (at 0 0 0) (axis 0 0 1)))
  (fix roller (- 5079.999999 x) (pivot (at 5080 0 0) (axis 0 1 0)) (pivot (at 5080 0 0) (axis 0 0 1)) (slide 1 0 0))
  (load compression (force -4.4482216152605N 0 0) (region (- 5079.999999 x))))
; Wide envelope requests the actual three-grid band, not agreement with target.
(require measurement (min-buckling-factor s > 0))

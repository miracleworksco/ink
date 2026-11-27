; ansys-vm96
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM96.html
; Temperature at centerline z=.25 ft; target 15.6 F from published table.
; Cylinder radius=length=.5 ft, k 1 BTU/(h ft F). Bottom and side 0 F, top 40 F.
; The target is a graphical reference, not exact to its last digit;
; no accuracy interval is published. Full 3D solid has the same axisymmetric PDE.
; Corners have discontinuous temperature boundary values; regions below
; assign the top face priority at its intersection with the side.
; Post-run exact Bessel series gives 15.35650235415488 F, versus
; the printed graph reading 15.6 F. Twenty and thirty modes agree.
; The original graph target remains ungraded.
; scoreboard: {"id":"ansys-vm96","physics":"thermal steady","target":-9.11111111111111,"target_kind":"converged reference","units":"C","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM96.html"],"quantity":"Temperature at centerline z=.25ft; target15.6F from published table.","qualification":"Graphical target has no independently accepted uncertainty","probe":"temperature","points":[[0,0,76.2]]}
(import "std")
(def part (move (cylinder-z 0.5ft 0.5ft) 0 0 0.25ft))
(body part (shape part))
(study s (shape part)
  (heat (k 1.7307346663713912W/m/K))
  (temperature bottom (- z 0.000001) 0degF)
  (temperature top (- 152.399999 z) 40degF)
  (temperature side (max (- 152.399999 (len2 x y)) (- z 152.399998)) 0degF))

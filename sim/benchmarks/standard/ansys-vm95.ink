; ansys-vm95
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM95.html
; Tip center temperature; square 1 in fin, length 8 in, base 100 F, ambient 0 F.
; k 25 BTU/(h ft F), lateralh 1 BTU/(h ft^2 F), insulated tip.
; Published 68.594 F is constant-section-temperature fin theory.
; Finite 3D section eigenvalues obey mu*tan(mu*b/2)=h/k.
; Their lowest eigenvalue is lower than fin sqrt(4h/(kb));
; that tends to increase the tip center temperature; higher modes/base layers
; remain unqualified. No pass/miss against the fin value is claimed.
; Post-run independent separated solid solution: tip 68.64347762641944 F,
; 0.04947762641944564 F above the printed fin target. Six and twelve
; transverse terms agree. Original idealisation remains ungraded in this score.
; scoreboard: {"id":"ansys-vm95","physics":"thermal steady","target":20.33,"target_kind":"idealisation","units":"C","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM95.html"],"quantity":"Tip center temperature; square1in fin, length8in, base100F, ambient0F.","qualification":"Finite-solid correction sign from separated cross-section eigenvalue; magnitude not yet accepted","probe":"temperature","points":[[203.2,0,0]]}
(import "std")
(def part (move (box 8in 1in 1in) 4in 0 0))
(body part (shape part))
(study s (shape part)
  (heat (k 43.26836665928478W/m/K))
  (temperature base (- x 0.000001) 100degF)
  (convection sides (max (- 0.000001 x) (- x 203.199999)) (coefficient 5.678263341113488W/m2/K) (ambient 0degF)))

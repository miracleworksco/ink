; ansys-vm28
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM28.html
; Temperature at y=.9141 m at 2000 s; published 98.9 C.
; Slab 0..1 m, lateral width.1 m (use.1 m adiabatic depth), k 54 W/mK,
; rho 7833 kg/m 3, cp 465J/kgK, initial 0 C, top convection 50 W/m^2K@1000 C.
; Bottom insulated: finite-domain approximation to semi-infinite solid.
; Exact halfspace formula erfc(eta)-exp(hd/k+(h/k)^2 alpha*t)
; *erfc(eta+h*sqrt(alpha*t)/k), multiplied 1000; d=1-y.
; Far-boundary difference at this time is exponentially small but not
; independently bounded here; this approximate target remains ungraded.
; scoreboard: {"id":"ansys-vm28","physics":"thermal transient","target":98.9,"target_kind":"idealisation","units":"C","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM28.html"],"quantity":"Temperature at y=.9141m at2000s; published98.9C.","qualification":"Semi-infinite target versus finite1m adiabatic slab not independently bounded","probe":"temperature","points":[[0,914.1,0]]}
(import "std")
(def part (move (box 100 1000 100) 0 500 0))
(body part (shape part))
(study s (shape part)
  (heat (k 54W/m/K) (density 7833kg/m3) (capacity 465J/kg/K))
  (convection top (- 999.999999 y) (coefficient 50W/m2/K) (ambient 1000degC))
  (transient (duration 2000s) (initial 0)))

; ansys-vm142
; Source for all dimensions, loads, material and printed target: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM142.html
; Axial stress at hole crown x 0, y.5 in, z 0: target+3018 psi.
; Square 12inx 12 in, hole diameter 1 in, thickness 1 in;E 30e6 psi, nu.3.
; Uniform tension 1000 psi on x ends. Quarter plane, half the thickness:
; x 0, y 0 and the midplane z 0 are mirrors; outer y and hole walls free.
; Target is a plane-stress stress-concentration formula. Thickness equals
; hole diameter, so no small thickness parameter justifies its 3D edge
; correction. A converged solid reference is required; sign and size
; remain unknown and no accuracy match or miss is counted.
; scoreboard: {"id":"ansys-vm142","physics":"static","target":20.808377510781025,"target_kind":"idealisation","units":"MPa","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM142.html"],"quantity":"Axial stress at hole crown x0,y.5in,z0: target+3018psi.","units_override":"MPa","qualification":"Plane-stress hole concentration; thickness/diameter1 needs accepted3D solid reference","probe":"stress_mpa.xx","points":[[0,12.7,0]]}
(import "std")
(def part (max (subtract (box 12in 12in 1in) (cylinder-z 0.5in 2in)) (neg x) (neg y) (neg z)))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0.3))
  (mirror sx x 0)
  (mirror sy y 0)
  (mirror sz z 0)
  (load end (pressure -6.894757293168MPa) (region (- 152.399999 x))))

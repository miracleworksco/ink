; ansys-vm25-pressure
; Input/target source: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM25.html
; Hoop stress sigma_yy at x 4 in, y 0, z 4 in: target+50000 psi (tension).
; Open-ended cylinder radii 4/8 in, E 30e6 psi, nu.3, internal pressure 30000 psi.
; Use an 8 in axial length. Quarter symmetry cuts are x=0 and y=0.
; Axial displacement is held at z=0 and free at the other end.
; Exact 3D Lame solution: A=p*a^2/(b^2-a^2), B=A*b^2,
; sigma_r=A-B/r^2, sigma_theta=A+B/r^2, uz=-2*nu*A*z/E.
; The radial/hoop fields satisfy this finite open-ended solid exactly;
; no shell or end-cap correction. Rotation subcase is excluded.
; scoreboard: {"id":"ansys-vm25-pressure","physics":"static","target":344.7378646584,"target_kind":"exact","units":"MPa","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM25.html"],"quantity":"Hoop stress sigma_yy at x4 in,y0,z4 in: target+50000 psi (tension).","target_rounding":0.003447378646584,"probe":"stress_mpa.yy","points":[[101.6,0,101.6]]}
(import "std")
(def part (max (extrude-z (ring-2d 4in 8in) 0 8in) (neg x) (neg y)))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0.3))
  (mirror sx x 0)
  (mirror sy y 0)
  (fix axial-datum (- z 0.000001) (dofs z) (support))
  (load bore (pressure 206.84271879503999MPa) (region (max (- (len2 x y) 101.600001) (- 0.000001 x) (- 0.000001 y) (- 0.000001 z) (- z 203.199999)))))

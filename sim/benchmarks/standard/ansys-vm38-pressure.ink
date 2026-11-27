; ansys-vm38-pressure
; Input/target source: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM38.html
; Hoop stress sigma_yy at x 4.4 in, y 0, z 4 in: target+18645 psi (rounded).
; Open-ended cylinder radii 4/8 in, E 30e6 psi, nu.3, internal pressure 12990 psi.
; Choose 8 in arbitrary axial length, with quarter symmetry and free axial
; strain. Hold axial displacement at one end of this 8 in slice.
; Exact 3D Lame stresses are A+B/r^2 in hoop, A-B/r^2 radial,
; A=p*a^2/(b^2-a^2), B=A*b^2. Axial displacement is -2*nu*A*z/E.
; No end-cap load. The pressure-only geometry is the same Lame family as
; VM 25, with different pressure and observation point; not independent
; physics coverage. Plastic/rotation subcases are excluded.
; scoreboard: {"id":"ansys-vm38-pressure","physics":"static","target":128.55274973111736,"target_kind":"exact","units":"MPa","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM38.html"],"quantity":"Hoop stress sigma_yy at x4.4 in,y0,z4 in: target+18645 psi (rounded).","target_rounding":0.003447378646584,"probe":"stress_mpa.yy","points":[[111.76,0,101.6]]}
(import "std")
(def part (max (extrude-z (ring-2d 4in 8in) 0 8in) (neg x) (neg y)))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0.3))
  (mirror sx x 0)
  (mirror sy y 0)
  (fix axial-datum (- z 0.000001) (dofs z) (support))
  (load bore (pressure 89.56289723825232MPa) (region (max (- (len2 x y) 101.600001) (- 0.000001 x) (- 0.000001 y) (- 0.000001 z) (- z 203.199999)))))

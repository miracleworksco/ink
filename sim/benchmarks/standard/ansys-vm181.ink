; ansys-vm181
; Input/target source: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM181.html
; First whole-plate natural frequency: first axisymmetric target 172.64 Hz.
; Disk radius 17 in, thickness.5 in, E 30e6 psi, nu.3;
; density .00073 lbf*s^2/in^4, converted with 1 lbf=4.4482216152605 N.
; Clamp the outer cylindrical rim; solve the whole plate, no symmetry.
; The higher targets 671.79 and 1505.07 Hz are higher axisymmetric modes,
; not the next whole-plate modes. Shear and rotary inertia lower thin-plate
; frequencies. Leading Mindlin estimate f/fK=[1+k_mode^2*(t^2/12+
; D/(kappa*G*t))]^-1/2, k_mode^2 about 10.2158/a^2, kappa 5/6.
; The 3D rim correction is unbounded here; no accuracy pass/miss assigned.
; scoreboard: {"id":"ansys-vm181","physics":"modal","target":172.64,"target_kind":"idealisation","units":"Hz","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM181.html"],"quantity":"First whole-plate natural frequency: first axisymmetric target172.64 Hz.","qualification":"Thinplate target; leading Mindlin shift downward, exact3D correction not accepted","metric":"min-frequency"}
(import "std")
(def part (cylinder-z 17in 0.5in))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0.3) (density 7801.433480086962kg/m3))
  (fix rim (- 431.799999 (len2 x y)))
  (modes 1))
; Wide envelope requests the actual three-grid band, not agreement with target.
(require measurement (min-frequency s > 0))

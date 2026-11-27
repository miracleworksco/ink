; ansys-vm78
; Source for all dimensions, loads, material and printed target: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM78.html
; Midspan center shear stress sigma_xz: target+7500 psi.
; Beam 10 in long, width 1 in along y, height 2 in along z, E 30e6 psi, nu 0.
; Root clamped;10000 lbf shear in+z on end. Source also asks Tsai-Wu
; failure index 225, which ink does not implement and is not replaced.
; Saint-Venant parabolic shear P*(h^2/4-z^2)/(2I) gives 7500 psi
; at midplane. The rigid-root 3D correction decays toward midspan,
; but no amplitude bound is accepted here; ungraded idealisation.
; scoreboard: {"id":"ansys-vm78","physics":"static","target":51.710679698759996,"target_kind":"idealisation","units":"MPa","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM78.html"],"quantity":"Midspan center shear stress sigma_xz: target+7500psi.","units_override":"MPa","qualification":"Parabolic shear target away from ends; clamp correction unbounded; failure-index subcase unsupported","probe":"stress_mpa.zx","points":[[127,0,0]]}
(import "std")
(def part (move (box 10in 1in 2in) 5in 0 0))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0))
  (fix root (- x 0.000001))
  (load end (force 0 0 44482.216152605N) (region (- 253.999999 x))))

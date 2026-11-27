; ansys-vm50
; Source for all dimensions, loads, material and printed target: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM50.html
; Lowest frequency of simply supported square 2 in beam, length 80 in.
; E 30e6 psi, density.000728 lbf*s^2/in^4, target 28.766 Hz.
; Nu unspecified; use.3 and keep solid comparison ungraded. Both end
; sections pivot in bending; one also slides along beam. No symmetry.
; Timoshenko leading estimate f/f_E=[1+(pi/L)^2*(I/A)*(1+E/(kG))]^-1/2
; with k 5/6 predicts lower frequency than Euler-Bernoulli;3D end effects
; not bounded. Density converts force units, not pound-mass.
; scoreboard: {"id":"ansys-vm50","physics":"modal","target":28.766,"target_kind":"idealisation","units":"Hz","study":"s","qualified":false,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM50.html"],"quantity":"Lowest frequency of simply supported square2in beam,length80in.","qualification":"Beam target; nu unspecified, shear/rotary shift predicted but finite-solid bound unavailable","metric":"min-frequency"}
(import "std")
(def part (move (box 80in 2in 2in) 40in 0 0))
(body part (shape part))
(study s (shape part)
  (elastic (E 206842.71879504MPa) (nu 0.3) (density 7780.0596897305595kg/m3))
  (fix pin (- x 0.000001) (pivot (at 0 0 0) (axis 0 1 0)) (pivot (at 0 0 0) (axis 0 0 1)))
  (fix roller (- 2031.999999 x) (pivot (at 2032 0 0) (axis 0 1 0)) (pivot (at 2032 0 0) (axis 0 0 1)) (slide 1 0 0))
  (modes 1))
; Wide envelope requests the actual three-grid band, not agreement with target.
(require measurement (min-frequency s > 0))

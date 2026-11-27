; solidworks-cantilever
; Input/target source: https://www.solidworks.com/sw/docs/Accurate_Design_Simulations_2018.pdf#page=11
; Signed tip-center y displacement: published -.004 in.
; Cantilever length 10 in, square 1 in, E 1e6 psi, nu.4; root clamped,
; transverse tip force 1 lbf downward. Source: whitepaper printed p 11.
; Euler-Bernoulli target omits shear. F*L/(kappa*G*A)=.0000336 in,
; +.84% in magnitude, with kappa 5/6 and G=E/[2(1+nu)]. The 3D rigid-root
; correction is not bounded, so the comparison is ungraded.
; SolidWorks prints -.004 in for both its solid and beam solutions.
; This rounded equality does not establish exactly zero error.
; scoreboard: {"id":"solidworks-cantilever","physics":"static","target":-0.1016,"target_kind":"idealisation","units":"mm","study":"s","qualified":false,"sources":["https://www.solidworks.com/sw/docs/Accurate_Design_Simulations_2018.pdf#page=11"],"quantity":"Signed tip-center y displacement: published -.004 in.","vendor_error":{"SolidWorks":"Rounded equality -.004in; exact error unavailable (hand-set mesh)"},"qualification":"Beam target plus shear estimate; clamp correction not bounded","probe":"displacement_mm.1","points":[[254,0,0]]}
(import "std")
(def part (move (box 10in 1in 1in) 5in 0 0))
(body part (shape part))
(study s (shape part)
  (elastic (E 6894.757293168MPa) (nu 0.4))
  (fix root (- x 0.000001))
  (load tip (force 0 -4.4482216152605N 0) (region (- 253.999999 x))))

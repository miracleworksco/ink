; afnor-sdll10
; Source for all dimensions, loads, material and printed target: https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v2/v2.02.010/index.html
; First bending frequency of clamped exponential-width beam, 143.303 Hz.
; Length .6 m, thickness .01 m, width .03*exp(-2*x_m)m;
; E 210 GPa, nu.3, rho 7800 kg/m 3; both ends clamped. Target from
; public SolidWorks AFNOR table p 66, supported by code_aster geometry
; and material definition. Thin beam target omits shear and rotary
; inertia, which lower bending frequency at leading order. Finite solid
; clamp response not bounded, so no pass/miss against beam number.
; scoreboard: {"id":"afnor-sdll10","physics":"modal","target":143.303,"target_kind":"idealisation","units":"Hz","study":"s","qualified":false,"sources":["https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v2/v2.02.010/index.html"],"quantity":"First bending frequency of clamped exponential-width beam,143.303Hz.","qualification":"Beam target, finite solid correction not bounded","vendor_error":{"SolidWorks":"1.97% solid;1.78% beam (hand-set meshes, AFNOR2019 p66)"},"metric":"min-frequency"}
(import "std")
(def part (max (neg x) (- x 600) (- (abs z) 5) (- (abs y) (* 15 (exp (* -0.002 x))))))
(body part (shape part))
(study s (shape part)
  (elastic (E 210000MPa) (nu 0.3) (density 7800kg/m3))
  (fix left (- x 0.000001))
  (fix right (- 599.999999 x))
  (modes 1))
; Wide envelope requests the actual three-grid band, not agreement with target.
(require measurement (min-frequency s > 0))

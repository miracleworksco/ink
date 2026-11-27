; afnor-sdll09
; Input/target source: https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v2/v2.02.009/index.html
; First bending frequency of tapered square cantilever, case 1.
; Length 1 m, root square.04 m, tip square.01 m, linear size taper;
; E 210 GPa and density 7800 kg/m^3 from code_aster reference problem.
; Source omits nu: retain ink's .3 assumption and leave comparison ungraded.
; Reference-solution formula lambda1=23.289 gives 55.519135565590965 Hz:
; f=lambda1*h_tip/(2*pi*L^2)*sqrt(E/(12*rho)).
; Source inconsistency: code_aster model A and SolidWorks
; print 54.18 Hz instead, matching a different E scale. Both targets remain
; quarantined; neither is chosen to score the measured value.
; Shear and rotary inertia lower beam frequencies at leading order;
; the tapered-mode correction and 3D root layer have not been bounded.
; scoreboard: {"id":"afnor-sdll09","physics":"modal","target":55.519135565590965,"target_kind":"idealisation","units":"Hz","study":"s","qualified":false,"sources":["https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v2/v2.02.009/index.html"],"quantity":"First bending frequency of tapered square cantilever, case1.","qualification":"Quarantined: reference formula and material imply55.519Hz, modelA prints54.18Hz; finite-solid correction also unqualified","metric":"min-frequency","published_reference_candidate":54.18,"quarantined":true}
(import "std")
(def part (max (- (abs y) (- 20 (* 0.015 x))) (- (abs z) (- 20 (* 0.015 x))) (neg x) (- x 1000)))
(body part (shape part))
(study s (shape part)
  (elastic (E 210000MPa) (nu 0.3) (density 7800kg/m3))
  (fix root (- x 0.000001))
  (modes 1))
; Wide envelope requests the actual three-grid band, not agreement with target.
(require measurement (min-frequency s > 0))

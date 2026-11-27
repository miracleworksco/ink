; afnor-tpla06
; Source for inputs and printed targets: https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.01.006/index.html
; Temperature at axisz.5 m; published fin formula 500*sinh(a*z)/sinh(a*L).
; Radius.01 m, length 1 m, k 33.33, h 10ambient 0, ends 0/500 C.
; Target assumes uniform temperature across section; source states<1%
; reference uncertainty. Lowest transverse root mu*J 1(muR)/J 0(muR)=h/k
; is slightly below fin a=sqrt(2h/(kR)); true center temperature is
; expected slightly higher. Boundary-layer contribution not bounded here.
; No pass or miss is assigned until that 3D correction is accepted.
; Post-run separated solid solution: center 10.414500988273232 C,
; 0.022889780117125724 C above fin theory. First omitted radial mode
; decays by exp(-190) at the midpoint. Frozen comparison remains ungraded.
; scoreboard: {"id":"afnor-tpla06","physics":"thermal steady","target":10.391611208156107,"target_kind":"idealisation","units":"C","study":"s","qualified":false,"sources":["https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.01.006/index.html"],"quantity":"Temperature at axisz.5m; published fin formula500*sinh(a*z)/sinh(a*L).","qualification":"Fin approximation source<1%; finite3D center correction not bounded","probe":"temperature","points":[[0,0,500]]}
(import "std")
(def part (move (cylinder-z 10 1000) 0 0 500))
(body part (shape part))
(study s (shape part)
  (heat (k 33.33W/m/K))
  (temperature bottom (- z 0.000001) 0)
  (temperature top (- 999.999999 z) 500degC)
  (convection side (max (- 9.999999 (len2 x y)) (- 0.000001 z) (- z 999.999999)) (coefficient 10W/m2/K) (ambient 0)))

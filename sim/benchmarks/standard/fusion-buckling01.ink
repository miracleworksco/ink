; fusion-buckling01
; All dimensions, properties, force, support, target and vendor error:
; https://help.autodesk.com/cloudhelp/ENU/Fusion-Simulate/files/GUID-B3839814-8098-42E5-90C4-C139252107EE.htm
; Square 25 mm, length 1000 mm, E=210000 MPa, nu=.3; fixed z=0,
; 1 N compression on z=1000. Asked: lowest linear buckling load, 16867 N.
; Target is beam theory: pi^2 E (25^4/12)/(4 L^2).
; A shear-flexible beam predicts P_E/(1+P_E/(k G A)), k=5/6;
; this is a lower correction, not a qualified converged 3D reference.
; No match/miss is scored until a solid reference qualifies the idealisation.
; scoreboard: {"id":"fusion-buckling01","physics":"buckling","target":16867,"target_kind":"idealisation","units":"N","metric":"min-buckling-factor","study":"s","qualified":false,"qualification":"Euler column; finite-solid correction not independently qualified","vendor_error":{"Fusion":"0.083% (hand-set 8 mm quadratic tet mesh)"}}
(import "std")
(def column (move (box 25 25 1000) 0 0 500))
(body column (shape column))
(study s (shape column) (elastic (E 210000MPa) (nu 0.3))
  (fix base (- z 0.000001))
  (load compression (force 0 0 -1N) (region (- 999.999999 z))))
(require measurement (min-buckling-factor s > 0))

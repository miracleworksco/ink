; nafems-fv52
; Source for all dimensions, loads, material and printed target: https://help.autodesk.com/cloudhelp/ENU/Fusion-Simulate/files/GUID-1ECFF5CF-9936-4E2F-BFEA-FED0D51E8140.htm
; First elastic out-of-plane frequency, fourth mode including 3 rigid modes.
; Plate 10mx 10mx 1 m, E 200 GPa, nu.3, rho 8000 kg/m 3. All bottom
; edges hold z only, as in the public Fusion solid reproduction.
; No artificial stabilisation is added: current ink may refuse the free
; planar rigid motions. Its min-frequency interface cannot select the
; fourth mode; if it solves, the primary comparison stays ungraded.
; Published solid target 45.897 Hz; Fusion 44.61 Hz, -2.80% with 1.67 m mesh.
; scoreboard: {"id":"nafems-fv52","physics":"modal","target":45.897,"target_kind":"converged reference","units":"Hz","study":"s","qualified":false,"sources":["https://help.autodesk.com/cloudhelp/ENU/Fusion-Simulate/files/GUID-1ECFF5CF-9936-4E2F-BFEA-FED0D51E8140.htm"],"quantity":"First elastic out-of-plane frequency, fourth mode including 3 rigid modes.","qualification":"Asked first elastic mode after three rigid motions; current min-frequency metric cannot select it","vendor_error":{"Fusion":"-2.80% (hand-set1.67m quadratic tets)"},"metric":"min-frequency"}
(import "std")
(def part (box 10000 10000 1000))
(body part (shape part))
(study s (shape part)
  (elastic (E 200000MPa) (nu 0.3) (density 8000kg/m3))
  (fix edges (max (+ z 499.999999) (- 4999.999999 (max (abs x) (abs y)))) (dofs z) (support))
  (modes 4))
; Wide envelope requests the actual three-grid band, not agreement with target.
(require measurement (min-frequency s > 0))

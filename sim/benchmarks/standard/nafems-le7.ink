; nafems-le7
; Source for all dimensions, loads, material and printed target: https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-le7.htm
; Outer axial stress sigma_zz at r 1.0125 m, z 1.4034 m; target 25.9 MPa.
; Mean radius 1 m, cylindrical length 1.5 m, hemisphere cap radius 1 m,
; wall .025 m, E 210 GPa, nu.3, internal pressure 1 MPa. Base holds z;
; axis regularity at apex follows from quarter solid symmetry.
; Build true constant-thickness offset cylinder plus hemisphere surfaces.
; Shell target omits radial stress and uses midsurface equilibrium;
; finite pressure-vessel correction scales t/R=.025, with an extra
; cylinder/cap bending layer at the asked position. Sign/magnitude of
; that local shift not independently qualified, so no pass/miss claimed.
; scoreboard: {"id":"nafems-le7","physics":"static","target":25.9,"target_kind":"idealisation","units":"MPa","study":"s","qualified":false,"sources":["https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-le7.htm"],"quantity":"Outer axial stress sigma_zz at r1.0125 m,z1.4034 m; target25.9 MPa.","units_override":"MPa","qualification":"Shell target, local finite-thickness correction not bounded","probe":"stress_mpa.zz","points":[[1012.5,0,1403.4]]}
(import "std")
(def part (max (neg x) (neg y) (subtract (min (move (cylinder-z 1012.5 1500) 0 0 750) (move (sphere 1012.5) 0 0 1500)) (min (move (cylinder-z 987.5 1500.000002) 0 0 749.999999) (move (sphere 987.5) 0 0 1500)))))
(body part (shape part))
(study s (shape part)
  (elastic (E 210000MPa) (nu 0.3))
  (fix base (- z 0.000001) (dofs z) (support))
  (mirror sx x 0)
  (mirror sy y 0)
  (load inside (pressure 1MPa) (region (max (- (select (- z 1500) (len2 x y) (len3 x y (- z 1500))) 987.500001) (- 0.000001 z) (- 0.000001 x) (- 0.000001 y)))))

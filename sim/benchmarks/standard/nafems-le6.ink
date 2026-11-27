; nafems-le6
; Source for all dimensions, loads, material and printed target: https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-le6.htm
; Maximum principal stress on bottom at rhombus center E; +.802 MPa.
; Side length 1 m, acute angle 30 deg, thickness .01 m; E 210 GPa, nu .3.
; All four boundary midsurface lines held in z, A held x/y, B held y.
; Pressure .7 kPa downward on top. Finite solid uses exactly the
; midsurface support lines, leaving edge rotations free as specified.
; Thin plate target omits transverse shear and normal stress; leading
; interior correction scales (t/span)^2 but the skew-corner layers
; are not independently bounded. No pass/miss against the shell number.
; scoreboard: {"id":"nafems-le6","physics":"static","target":0.802,"target_kind":"idealisation","units":"MPa","study":"s","qualified":false,"sources":["https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-le6.htm"],"quantity":"Maximum principal stress on bottom at rhombus center E; +.802 MPa.","principal_max":true,"units_override":"MPa","qualification":"Thin plate stress, finite skew solid correction not bounded","probe":"stress_mpa","points":[[933.0127018922194,250,-5]]}
(import "std")
(def part (max (neg y) (- y 500) (neg (- (* 0.5 x) (* 0.8660254037844386 y))) (- (- (* 0.5 x) (* 0.8660254037844386 y)) 500) (- (abs z) 5)))
(body part (shape part))
(study s (shape part)
  (elastic (E 210000MPa) (nu 0.3))
  (fix edges (max (- (abs z) 0.000001) (min (- y 0.000001) (- 499.999999 y) (- (- (* 0.5 x) (* 0.8660254037844386 y)) 0.000001) (- 499.999999 (- (* 0.5 x) (* 0.8660254037844386 y))))) (dofs z))
  (fix A (- (len3 x y z) 0.000001) (dofs x y))
  (fix B (- (len3 (- x 1000) y z) 0.000001) (dofs y) (support))
  (load pressure (pressure 0.0007MPa) (region (- 4.999999 z))))

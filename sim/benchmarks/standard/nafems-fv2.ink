; nafems-fv2
; Source for all dimensions, loads, material and printed target: https://help.autodesk.com/cloudhelp/ENU/Fusion-Simulate/files/GUID-1B36CACE-A2BB-4202-9042-2D44E42C6284.htm
; Lowest in-plane mode of pin-ended eight-spoke double cross;11.336 Hz.
; Each spoke 5 m long, square .125 m; overall axis span 10 m.
; E 200 GPa, rho 8000 kg/m 3, nu.3. Solid union of four 10 m bars,
; rotated 0, 45, 90, -45 degrees, shared central junction. End sections
; pivot about z, center translations held; bottom face holds z as in
; Fusion public solid reproduction. Grooves are replaced by ideal pins
; at exactly the section centroid, matching the stated beam supports.
; Shear/rotary inertia lower beam frequencies at order(width/length)^2;
; finite solid junction and restrained back face Poisson response can
; offset this. No accepted finite-solid bound is available.
; scoreboard: {"id":"nafems-fv2","physics":"modal","target":11.336,"target_kind":"idealisation","units":"Hz","study":"s","qualified":false,"sources":["https://help.autodesk.com/cloudhelp/ENU/Fusion-Simulate/files/GUID-1B36CACE-A2BB-4202-9042-2D44E42C6284.htm"],"quantity":"Lowest in-plane mode of pin-ended eight-spoke double cross;11.336Hz.","qualification":"Beam target, pin substitution and finite solid junction/backface correction unqualified","vendor_error":{"Fusion":"0.053% (hand-set31.25mm quadratic tets; first mode only)"},"metric":"min-frequency"}
(import "std")
(def part (min (box 10000 125 125) (box 125 10000 125) (rotate-z (box 10000 125 125) 45) (rotate-z (box 10000 125 125) -45)))
(body part (shape part))
(study s (shape part)
  (elastic (E 200000MPa) (nu 0.3) (density 8000kg/m3))
  (fix plane (+ z 62.499999) (dofs z) (support))
  (fix pin0 (- 4999.999999 (+ (* 1.0 x) (* 0.0 y))) (dofs x y) (pivot (at 5000.0 0.0 0) (axis 0 0 1)))
  (fix pin1 (- 4999.999999 (+ (* 0.7071067811865476 x) (* 0.7071067811865475 y))) (dofs x y) (pivot (at 3535.533905932738 3535.5339059327375 0) (axis 0 0 1)))
  (fix pin2 (- 4999.999999 (+ (* 6.123233995736766e-17 x) (* 1.0 y))) (dofs x y) (pivot (at 3.061616997868383e-13 5000.0 0) (axis 0 0 1)))
  (fix pin3 (- 4999.999999 (+ (* -0.7071067811865475 x) (* 0.7071067811865476 y))) (dofs x y) (pivot (at -3535.5339059327375 3535.533905932738 0) (axis 0 0 1)))
  (fix pin4 (- 4999.999999 (+ (* -1.0 x) (* 1.2246467991473532e-16 y))) (dofs x y) (pivot (at -5000.0 6.123233995736766e-13 0) (axis 0 0 1)))
  (fix pin5 (- 4999.999999 (+ (* -0.7071067811865477 x) (* -0.7071067811865475 y))) (dofs x y) (pivot (at -3535.5339059327384 -3535.5339059327375 0) (axis 0 0 1)))
  (fix pin6 (- 4999.999999 (+ (* -1.8369701987210297e-16 x) (* -1.0 y))) (dofs x y) (pivot (at -9.18485099360515e-13 -5000.0 0) (axis 0 0 1)))
  (fix pin7 (- 4999.999999 (+ (* 0.7071067811865474 x) (* -0.7071067811865477 y))) (dofs x y) (pivot (at 3535.5339059327366 -3535.5339059327384 0) (axis 0 0 1)))
  (modes 1))
; Wide envelope requests the actual three-grid band, not agreement with target.
(require measurement (min-frequency s > 0))

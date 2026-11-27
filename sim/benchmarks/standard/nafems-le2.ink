; nafems-le2
; Source for all dimensions, loads, material and printed target: https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-le2.htm
; Outer surface tangential stress at E, target +60 MPa.
; Cylinder mean radius 1 m, wall .01 m, axial width .5 m, arc angle 30 deg.
; E=210 GPa, nu=.3. Clamp the theta=0 end; hold axial displacement
; on both axial edges. At theta=30 deg apply uniform normal edge moment
; 1000 N m per m of edge, total 500 N m. Positive couple chosen so
; outer surface stress is tensile. E lies at axial z=.3 m, theta=20 deg.
; This is shell pure bending: 6M/t^2=60 MPa. Finite curved solid stress
; is asymmetric across wall (curved-beam correction of order t/R=.01);
; no accepted bound on that shift or the held-edge layer is available.
; scoreboard: {"id":"nafems-le2","physics":"static","target":60,"target_kind":"idealisation","units":"MPa","study":"s","qualified":false,"sources":["https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-le2.htm"],"quantity":"Outer surface tangential stress at E, target +60 MPa.","tensor_direction":[-0.3420201433256687,0.9396926207859084,0],"units_override":"MPa","qualification":"Shell stress, finite curvature and axial edge correction unbounded","probe":"stress_mpa","points":[[944.391083889838,343.73024404229704,300]]}
(import "std")
(def part (max (extrude-z (ring-2d 995 1005) 0 500) (neg y) (- (* 0.8660254037844386 y) (* 0.5 x))))
(body part (shape part))
(study s (shape part)
  (elastic (E 210000MPa) (nu 0.3))
  (fix clamp (- y 0.000001))
  (fix left (- z 0.000001) (dofs z) (support))
  (fix right (- 499.999999 z) (dofs z) (support))
  (load moment (torque 0 0 500000N*mm) (region (- (- (* 0.5 x) (* 0.8660254037844386 y)) 0.000001))))

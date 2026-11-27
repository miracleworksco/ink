; nafems-le1
; Input/target source: https://help.autodesk.com/cloudhelp/ENU/Fusion-Simulate/files/GUID-D258911F-A85E-44C3-A510-74585BA992A2.htm
; Tangential stress sigma_yy at D=(2000, 0, 0) mm: target +92.7 MPa.
; Quarter elliptical ring: inner semiaxes 2000/1000 mm, outer 3250/2750 mm,
; thickness 100 mm, E 210000 MPa, nu.3. Outer pressure 10 MPa outward.
; Symmetry: x 0, y 0 and the midsurface z 0 are mirrors (half the thickness).
; Target is plane stress. Finite 3D stress varies through the thickness;
; interior corrections scale as(thickness/curvature radius)^2, but the
; curved free-edge layer is not bounded here. No match or miss is scored.
; The outer-wall region's rims cross the cuts and end faces at an angle, not
; along them (a rim lying in a face, its integral missed most of its area):
; its field is z, x, y and 50 - z against 1500 mm (1 - E), E the outer
; ellipse's (x/3250)^2 + (y/2750)^2, which grows from the wall at 2/3250 to
; 2/2750 per mm: 1500 mm makes it about the distance in from the wall.
; scoreboard: {"id":"nafems-le1","physics":"static","target":92.7,"target_kind":"idealisation","units":"MPa","study":"s","qualified":false,"sources":["https://help.autodesk.com/cloudhelp/ENU/Fusion-Simulate/files/GUID-D258911F-A85E-44C3-A510-74585BA992A2.htm"],"quantity":"Tangential stress sigma_yy at D=(2000,0,0) mm: target +92.7 MPa.","vendor_error":{"Fusion":"0.14% (hand-set150mm quadratic tets)"},"qualification":"Plane-stress target, finite100mm3D edge correction unqualified","probe":"stress_mpa.yy","points":[[2000,0,0]]}
(import "std")
(def part
  (max (- (+ (square (/ x 3250)) (square (/ y 2750))) 1)
       (- 1 (+ (square (/ x 2000)) (square (/ y 1000))))
       (neg x) (neg y) (neg z) (- (abs z) 50)))
(body part (shape part))
(study s (shape part)
  (elastic (E 210000MPa) (nu 0.3))
  (mirror sx x 0)
  (mirror sy y 0)
  (mirror sz z 0)
  (load outer (pressure -10MPa) (region (max (+ (- z 50) (* 1500mm (- 1 (+ (square (/ x 3250)) (square (/ y 2750)))))) (- (* 1500mm (- 1 (+ (square (/ x 3250)) (square (/ y 2750))))) z) (- (* 1500mm (- 1 (+ (square (/ x 3250)) (square (/ y 2750))))) x) (- (* 1500mm (- 1 (+ (square (/ x 3250)) (square (/ y 2750))))) y)))))

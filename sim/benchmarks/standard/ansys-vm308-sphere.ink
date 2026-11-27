; ansys-vm308-sphere
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM308.html
; Temperature at radius.110 m; target 70.267 C.
; Inner/outer radii.1/.12 m; inner 100 C, outer 25 C; linear conductivity 7 to 2 W/(m K).
; Exact radial solution 100-75*I(r)/I(ro), I(r)=integral_ri^r ds/(s^n*k(s));
; n=1 for cylinder, n=2 for sphere. Cylinder ends insulated: finite length 240 mm
; is an exact extrusion of the infinite-cylinder equation, not an end correction.
; scoreboard: {"id":"ansys-vm308-sphere","physics":"thermal steady","target":70.267,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM308.html"],"quantity":"Temperature at radius.110m; target70.267C.","target_rounding":0.0005,"probe":"temperature","points":[[110,0,0]]}
(import "std")
(def part (subtract (sphere 120) (sphere 100)))
(body part (shape part))
(study s (shape part)
  (heat (k (* 1W/m/K (- 7 (* 0.25 (- (len3 x y z) 100))))))
  (temperature inside (- (len3 x y z) 100.000001) 100degC)
  (temperature outside (- 119.999999 (len3 x y z)) 25degC))

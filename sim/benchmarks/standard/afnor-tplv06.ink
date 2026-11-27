; afnor-tplv06
; Source for inputs and printed targets: https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.04.006/index.html
; Temperature at radius 1.5 m; target 32.5 C.
; Spherical annulus ri 1 m, ro 2 m, k 1, Q 100 W/m^3, both walls 20 C.
; Exact 3D T=20+Q/(6k)*[(b^2-a^2)*(1/a-1/r)/(1/a-1/b)-(r^2-a^2)].
; Use the whole sphere; total source 100*(4pi/3)*(2^3-1^3) W.
; scoreboard: {"id":"afnor-tplv06","physics":"thermal steady","target":32.5,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.04.006/index.html"],"quantity":"Temperature at radius1.5m; target32.5C.","probe":"temperature","points":[[1500,0,0]]}
(import "std")
(def part (subtract (sphere 2000) (sphere 1000)))
(body part (shape part))
(study s (shape part)
  (heat (k 1W/m/K))
  (temperature inside (- (len3 x y z) 1000.000001) 20degC)
  (temperature outside (- 1999.999999 (len3 x y z)) 20degC)
  (source generated (power 2932.153143350474W) (region -1)))

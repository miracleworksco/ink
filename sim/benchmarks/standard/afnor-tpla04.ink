; afnor-tpla04
; Source for inputs and printed targets: https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.01.004/index.html
; Temperature at r 1.5 m; target 32.62 C rounded.
; Annulus 1..2 m, k 1 W/mK, Q 100 W/m^3, both surfaces 20 C, insulated ends.
; Length unspecified: use 1 m; invariant radial solution is exact at any length.
; T=20+Q/(4k)*[(b^2-a^2)ln(r/a)/ln(b/a)-(r^2-a^2)].
; scoreboard: {"id":"afnor-tpla04","physics":"thermal steady","target":32.62,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.01.004/index.html"],"quantity":"Temperature at r1.5m; target32.62C rounded.","target_rounding":0.005,"probe":"temperature","points":[[1500,0,500]]}
(import "std")
(def part (extrude-z (ring-2d 1000 2000) 0 1000))
(body part (shape part))
(study s (shape part)
  (heat (k 1W/m/K))
  (temperature inside (- (len2 x y) 1000.000001) 20degC)
  (temperature outside (- 1999.999999 (len2 x y)) 20degC)
  (source generated (power 942.4777960769379W) (region -1)))

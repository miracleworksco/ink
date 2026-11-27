; afnor-tpla01
; Source for inputs and printed targets: https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.01.001/index.html
; Temperature at r.32 m; target 66.506 C; ri.30 m, ro.35 m, length.1 m, k 1.
; Inner 100 C, outer h 500@17.03444 C, ends insulated. Original lower-half inner
; boundary 100 C and upper-half incoming flux 1729.9091 W/m^2 produce the same
; radial logarithmic solution; here impose 100 C on whole bore. This is an
; analytically equivalent boundary substitution, not the original mixed BC test.
; The removed flux differs only at the printed input rounding. T(r)=Ti-C ln(r/ri),
; C=(Ti-Ta)/(ln(ro/ri)+1/(h ro)).
; scoreboard: {"id":"afnor-tpla01","physics":"thermal steady","target":66.506,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.01.001/index.html"],"quantity":"Temperature at r.32m; target66.506C; ri.30m,ro.35m,length.1m,k1.","adaptation":"Equivalent inner-temperature boundary, original mixed-flux coverage not claimed","target_rounding":0.0005,"probe":"temperature","points":[[320,0,50]]}
(import "std")
(def part (extrude-z (ring-2d 300 350) 0 100))
(body part (shape part))
(study s (shape part)
  (heat (k 1W/m/K))
  (temperature inside (- (len2 x y) 300.000001) 100degC)
  (convection outside (- 349.999999 (len2 x y)) (coefficient 500W/m2/K) (ambient 17.03444degC)))

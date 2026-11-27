; afnor-ttll01
; Source for inputs and printed targets: https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.21.001/index.html
; Temperature x.8 m at.2 s, target evaluated from cited odd Fourier series.
; Slab 0..2 m, k 1, rho*cp 1J/m^3K, initial 0, ends 100 C.
; Use adiabatic 2 m x 2 m section: dimensions not specified for transverse
; invariant direction, no effect on exact 3D solution. Density 1000 kg/m^3
; and capacity.001J/kgK give exactly the stated product.
; T=100*[1-4/pi*sum_odd sin(n*pi*x/2)/n*exp(-(n*pi/2)^2*t)].
; scoreboard: {"id":"afnor-ttll01","physics":"thermal transient","target":26.3672814896455,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.21.001/index.html"],"quantity":"Temperature x.8m at.2s, target evaluated from cited odd Fourier series.","probe":"temperature","points":[[800,0,0]]}
(import "std")
(def part (move (box 2000 2000 2000) 1000 0 0))
(body part (shape part))
(study s (shape part)
  (heat (k 1W/m/K) (density 1000kg/m3) (capacity 0.001J/kg/K))
  (temperature left (- x 0.000001) 100degC)
  (temperature right (- 1999.999999 x) 100degC)
  (transient (duration 0.2s) (initial 0)))

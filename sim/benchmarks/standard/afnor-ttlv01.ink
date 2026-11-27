; afnor-ttlv01
; Source for inputs and printed targets: https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.25.001/index.html
; Center temperature at 2400 s; source Gurney-Lurie chart 956 C with<2% uncertainty.
; Sphere diameter.2 m, k 48.822, rho 7200, cp 669, initial 20 C,
; outer h 232.5 W/m^2K with ambient 1000 C. Exact finite sphere geometry.
; Target is a chart reading, not a converged reference despite this
; schema category; independent series qualification remains outstanding.
; Post-run independent sphere series gives 952.8107813281413 C,
; within the chart uncertainty; ten and twenty modes agree. Frozen chart
; target remains ungraded; a chart reading is not a converged numerical reference.
; scoreboard: {"id":"afnor-ttlv01","physics":"thermal transient","target":956,"target_kind":"converged reference","units":"C","study":"s","qualified":false,"sources":["https://docaster-en-codeaster-doc-4113bb35b8213385c636b94f1536856d7f617.gitlab.io/manuals/man_v/v4/v4.25.001/index.html"],"quantity":"Center temperature at2400s; source Gurney-Lurie chart956C with<2% uncertainty.","qualification":"Published chart target<2% uncertainty; no accepted independent exact series yet","probe":"temperature","points":[[0,0,0]]}
(import "std")
(def part (sphere 100))
(body part (shape part))
(study s (shape part)
  (heat (k 48.822W/m/K) (density 7200kg/m3) (capacity 669J/kg/K))
  (convection surface -1 (coefficient 232.5W/m2/K) (ambient 1000degC))
  (transient (duration 2400s) (initial 20degC)))

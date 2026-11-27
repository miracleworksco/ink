; ansys-vm115
; Source for inputs and printed targets: https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM115.html
; Final center maximum temperature at 720 s; printed 132.10 F (three Fourier terms).
; Slab 1 ft thick, area 1 ft^2, adiabatic sides, initial 60 F, ends 32 F for t>0.
; k 20 BTU/(h ft F), density 500lb/ft^3, capacity.2 BTU/(lb F), Q 40000 BTU/(h ft^3).
; Exact 3D invariant solution uses odd cosine modes; initial excess 28 F,
; steady excessQ*(L^2-x^2)/(2k), half-thicknessL=.5 ft.
; Default 100 time steps: no steps/res override. Peak over history equals
; final center since it rises from 60 F toward the heated steady state.
; scoreboard: {"id":"ansys-vm115","physics":"thermal transient","target":55.611111111111114,"target_kind":"exact","units":"C","study":"s","qualified":true,"sources":["https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/ans_vm/Hlp_V_VM115.html"],"quantity":"Final center maximum temperature at720s; printed132.10F (three Fourier terms).","target_rounding":0.002777777777777778,"metric":"max-temperature"}
(import "std")
(def part (box 1ft 1ft 1ft))
(body part (shape part))
(study s (shape part)
  (heat (k 34.61469332742782W/m/K) (density 8009.231686980069kg/m3) (capacity 837.3600000000001J/kg/K))
  (temperature left (+ x 152.399999) 32degF)
  (temperature right (- 152.399999 x) 32degF)
  (source generated (power 11722.84280688889W) (region -1))
  (transient (duration 720s) (initial 60degF)))
; This wide measurement envelope requests the reported grid band; it is not a target test.
(require measurement (max-temperature s < 1000000degC))

; Materials and fluids: plain IR data, engineering reference inputs.
; Every scalar cites a source ID below; full title/maker/revision/URL/hash and
; per-property conditions are in std/material-sources.json, following catalog
; source manifests. UNSOURCED is an explicit open gap, not validated data.
; Units: E/strength MPa, density g/mm3 (written kg/m3), k/conductivity W/m K,
; capacity J/kg K, viscosity Pa s. Temper, product form and temperature matter.
; Strength is a specification minimum only where labelled so. Datasheet typical
; strength is NOT a guaranteed allowable: a claim that reads one by get is an
; estimate, naming it (typical_strengths). Evidence labels are returned by
; (material NAME).
; Printed values apply only to the named process/direction.
; Compatibility aliases resolve silently; lookup retains the requested name and
; adds canonical-name. Other spelling hints remain refused. See std/README.md.

; pla: Prusament PLA, FFF horizontal ISO 527-1 specimens; MK3, 0.20 mm layers, 100% rectilinear infill, 2 perimeters. See source p3 for full process; not across-layer Z strength.
(def pla-prusament-fff-horizontal (data
  (condition "Prusament PLA, FFF horizontal ISO 527-1 specimens; MK3, 0.20 mm layers, 100% rectilinear infill, 2 perimeters. See source p3 for full process; not across-layer Z strength.")
  (evidence "std/material-sources.json#pla-prusament-fff-horizontal")
  (strength-kind "yield: datasheet typical")
  (E 2300) ; prusa-pla; typical; p2, tensile modulus; E ±100 MPa; yield ±3 MPa
  (nu 0.36) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1240kg/m3) ; prusa-pla; typical; p2, ISO 1183; filament/material density, not infill-adjusted.
  (yield 51) ; prusa-pla; datasheet typical; p2, tensile yield; E ±100 MPa; yield ±3 MPa
  (k 0.13) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 1800) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def pla pla-prusament-fff-horizontal) ; accepted compatibility alias

; petg: Prusament PETG, FFF horizontal ISO 527-1 specimens; MK3, 0.20 mm layers, 100% rectilinear infill, 2 perimeters. See source p3 for full process; not across-layer Z strength.
(def petg-prusament-fff-horizontal (data
  (condition "Prusament PETG, FFF horizontal ISO 527-1 specimens; MK3, 0.20 mm layers, 100% rectilinear infill, 2 perimeters. See source p3 for full process; not across-layer Z strength.")
  (evidence "std/material-sources.json#petg-prusament-fff-horizontal")
  (strength-kind "yield: datasheet typical")
  (E 1500) ; prusa-petg; typical; p2, tensile modulus; E ±100 MPa; yield ±2 MPa
  (nu 0.38) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1270kg/m3) ; prusa-petg; typical; p2, ISO 1183; filament/material density, not infill-adjusted.
  (yield 47) ; prusa-petg; datasheet typical; p2, tensile yield; E ±100 MPa; yield ±2 MPa
  (k 0.2) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 1200) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def petg petg-prusament-fff-horizontal) ; accepted compatibility alias

; abs: Ensinger TECARAN ABS grey, machined stock shapes; specimens directly after machining under German standard climate. Datasheet values, not valid for prints.
(def abs-tecaran-stock (data
  (condition "Ensinger TECARAN ABS grey, machined stock shapes; specimens directly after machining under German standard climate. Datasheet values, not valid for prints.")
  (evidence "std/material-sources.json#abs-tecaran-stock")
  (strength-kind "yield: datasheet typical")
  (E 1700) ; ensinger-manual; typical; p88; tensile modulus ISO 527-2, density ISO 1183.
  (nu 0.35) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1040kg/m3) ; ensinger-manual; typical; p88; tensile modulus ISO 527-2, density ISO 1183.
  (yield 32) ; ensinger-manual; datasheet typical; p88, ISO 527-2; maker explicitly says these are not minimum values.
  (k 0.17) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 1400) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def abs abs-tecaran-stock) ; accepted compatibility alias

; nylon-pa12: Ensinger TECAMID 12 natural, machined stock shapes; specimens directly after machining under German standard climate. Datasheet values, not valid for prints.
(def pa12-tecamid12-stock (data
  (condition "Ensinger TECAMID 12 natural, machined stock shapes; specimens directly after machining under German standard climate. Datasheet values, not valid for prints.")
  (evidence "std/material-sources.json#pa12-tecamid12-stock")
  (strength-kind "yield: datasheet typical")
  (E 1800) ; ensinger-manual; typical; p91; tensile modulus ISO 527-2, density ISO 1183.
  (nu 0.39) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1020kg/m3) ; ensinger-manual; typical; p91; tensile modulus ISO 527-2, density ISO 1183.
  (yield 54) ; ensinger-manual; datasheet typical; p91, ISO 527-2; maker explicitly says these are not minimum values.
  (k 0.3) ; ensinger-manual; typical; p91, ISO 22007-4:2008; specific heat converted J/g K to J/kg K.
  (capacity 1800) ; ensinger-manual; typical; p91, ISO 22007-4:2008; specific heat converted J/g K to J/kg K.
))
(def nylon-pa12 pa12-tecamid12-stock) ; accepted compatibility alias

; polycarbonate: Ensinger TECANAT natural, machined stock shapes; specimens directly after machining under German standard climate. Datasheet values, not valid for prints.
(def pc-tecanat-stock (data
  (condition "Ensinger TECANAT natural, machined stock shapes; specimens directly after machining under German standard climate. Datasheet values, not valid for prints.")
  (evidence "std/material-sources.json#pc-tecanat-stock")
  (strength-kind "yield: datasheet typical")
  (E 2200) ; ensinger-manual; typical; p92; tensile modulus ISO 527-2, density ISO 1183.
  (nu 0.37) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1190kg/m3) ; ensinger-manual; typical; p92; tensile modulus ISO 527-2, density ISO 1183.
  (yield 69) ; ensinger-manual; datasheet typical; p92, ISO 527-2; maker explicitly says these are not minimum values.
  (k 0.25) ; ensinger-manual; typical; p92, ISO 22007-4:2008; specific heat converted J/g K to J/kg K.
  (capacity 1300) ; ensinger-manual; typical; p92, ISO 22007-4:2008; specific heat converted J/g K to J/kg K.
))
(def polycarbonate pc-tecanat-stock) ; accepted compatibility alias

; polypropylene: Ensinger tecafine-pp, unfilled machined stock shapes; datasheet values, not valid for prints.
(def polypropylene-tecafine-pp-stock (data
  (condition "Ensinger tecafine-pp, unfilled machined stock shapes; datasheet values, not valid for prints.")
  (evidence "std/material-sources.json#polypropylene-tecafine-pp-stock")
  (strength-kind "yield: datasheet typical")
  (E 1400) ; ensinger-pp; typical; ISO 527-2 tensile modulus
  (nu 0.42) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 900kg/m3) ; ensinger-pp; typical; ISO 1183
  (yield 32) ; ensinger-pp; datasheet typical; ISO 527-2 tensile strength at yield; no guaranteed minimum.
  (k 0.2) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 1920) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def polypropylene polypropylene-tecafine-pp-stock) ; accepted compatibility alias

; hdpe: Ensinger tecafine-pe300, unfilled machined stock shapes; datasheet values, not valid for prints.
(def hdpe-tecafine-pe300-stock (data
  (condition "Ensinger tecafine-pe300, unfilled machined stock shapes; datasheet values, not valid for prints.")
  (evidence "std/material-sources.json#hdpe-tecafine-pe300-stock")
  (strength-kind "yield: datasheet typical")
  (E 1100) ; ensinger-hdpe; typical; ISO 527-2 tensile modulus
  (nu 0.42) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 960kg/m3) ; ensinger-hdpe; typical; ISO 1183
  (yield 23) ; ensinger-hdpe; datasheet typical; ISO 527-2 tensile strength at yield; no guaranteed minimum.
  (k 0.46) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 1900) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def hdpe hdpe-tecafine-pe300-stock) ; accepted compatibility alias

; pet: Ensinger tecapet-white, unfilled machined stock shapes; datasheet values, not valid for prints.
(def pet-tecapet-white-stock (data
  (condition "Ensinger tecapet-white, unfilled machined stock shapes; datasheet values, not valid for prints.")
  (evidence "std/material-sources.json#pet-tecapet-white-stock")
  (strength-kind "yield: datasheet typical")
  (E 3100) ; ensinger-pet; typical; ISO 527-2 tensile modulus
  (nu 0.37) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1360kg/m3) ; ensinger-pet; typical; ISO 1183
  (yield 79) ; ensinger-pet; datasheet typical; ISO 527-2 tensile strength at yield; no guaranteed minimum.
  (k 0.31) ; ensinger-pet; typical; ISO 22007-4:2008
  (capacity 1200) ; ensinger-pet; typical; ISO 22007-4:2008; J/g K converted to J/kg K.
))
(def pet pet-tecapet-white-stock) ; accepted compatibility alias

; nylon-pa66: Ensinger TECAMID 66 natural, machined stock shapes; specimens directly after machining under German standard climate. Datasheet values, not valid for prints.
(def pa66-tecamid66-stock (data
  (condition "Ensinger TECAMID 66 natural, machined stock shapes; specimens directly after machining under German standard climate. Datasheet values, not valid for prints.")
  (evidence "std/material-sources.json#pa66-tecamid66-stock")
  (strength-kind "yield: datasheet typical")
  (E 3500) ; ensinger-manual; typical; p90; tensile modulus ISO 527-2, density ISO 1183.
  (nu 0.39) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1150kg/m3) ; ensinger-manual; typical; p90; tensile modulus ISO 527-2, density ISO 1183.
  (yield 84) ; ensinger-manual; datasheet typical; p90, ISO 527-2; maker explicitly says these are not minimum values.
  (k 0.36) ; ensinger-manual; typical; p90, ISO 22007-4:2008; specific heat converted J/g K to J/kg K.
  (capacity 1500) ; ensinger-manual; typical; p90, ISO 22007-4:2008; specific heat converted J/g K to J/kg K.
))
(def nylon-pa66 pa66-tecamid66-stock) ; accepted compatibility alias

; pom: Ensinger TECAFORM AH natural, machined stock shapes; specimens directly after machining under German standard climate. Datasheet values, not valid for prints.
(def pom-c-tecaform-ah-stock (data
  (condition "Ensinger TECAFORM AH natural, machined stock shapes; specimens directly after machining under German standard climate. Datasheet values, not valid for prints.")
  (evidence "std/material-sources.json#pom-c-tecaform-ah-stock")
  (strength-kind "yield: datasheet typical")
  (E 2800) ; ensinger-manual; typical; p88; tensile modulus ISO 527-2, density ISO 1183.
  (nu 0.35) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1410kg/m3) ; ensinger-manual; typical; p88; tensile modulus ISO 527-2, density ISO 1183.
  (yield 67) ; ensinger-manual; datasheet typical; p88, ISO 527-2; maker explicitly says these are not minimum values.
  (k 0.39) ; ensinger-manual; typical; p88, ISO 22007-4:2008; specific heat converted J/g K to J/kg K.
  (capacity 1400) ; ensinger-manual; typical; p88, ISO 22007-4:2008; specific heat converted J/g K to J/kg K.
))
(def pom pom-c-tecaform-ah-stock) ; accepted compatibility alias

; ptfe: Ensinger tecaflon-ptfe, unfilled machined stock shapes; datasheet values, not valid for prints.
(def ptfe-tecaflon-ptfe-stock (data
  (condition "Ensinger tecaflon-ptfe, unfilled machined stock shapes; datasheet values, not valid for prints.")
  (evidence "std/material-sources.json#ptfe-tecaflon-ptfe-stock")
  (E 500) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (nu 0.46) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 2150kg/m3) ; ensinger-ptfe; typical; ISO 1183
  ; removed yield: UNSOURCED; unsourced; Removed: maker reports tensile strength, not yield; legacy yield 10 MPa had no source.
  (k 0.2) ; ensinger-ptfe; typical; ISO 22007-4:2008
  (capacity 1000) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def ptfe ptfe-tecaflon-ptfe-stock) ; accepted compatibility alias

; peek: Ensinger tecapeek-natural, unfilled machined stock shapes; datasheet values, not valid for prints.
(def peek-tecapeek-natural-stock (data
  (condition "Ensinger tecapeek-natural, unfilled machined stock shapes; datasheet values, not valid for prints.")
  (evidence "std/material-sources.json#peek-tecapeek-natural-stock")
  (strength-kind "yield: datasheet typical")
  (E 4200) ; ensinger-peek; typical; ISO 527-2 tensile modulus
  (nu 0.4) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1310kg/m3) ; ensinger-peek; typical; ISO 1183
  (yield 116) ; ensinger-peek; datasheet typical; ISO 527-2 tensile strength at yield; no guaranteed minimum.
  (k 0.27) ; ensinger-peek; typical; ISO 22007-4:2008
  (capacity 1100) ; ensinger-peek; typical; ISO 22007-4:2008; J/g K converted to J/kg K.
))
(def peek peek-tecapeek-natural-stock) ; accepted compatibility alias

; nylon-pa6-mc901: MC901 cast PA6, absolutely dry reference values; KHK Tables 11.1–11.2. Not printed.
(def nylon-pa6-mc901 (data
  (condition "MC901 cast PA6, absolutely dry reference values; KHK Tables 11.1–11.2. Not printed.")
  (evidence "std/material-sources.json#nylon-pa6-mc901")
  (E 3432) ; khk-plastic; typical; Table 11.1: ASTM D638 E, D792 density, nu no standard; Table 11.2: ASTM C177 k, capacity 1.67 kJ/kg K converted to J/kg K.
  (nu 0.4) ; khk-plastic; typical; Table 11.1: ASTM D638 E, D792 density, nu no standard; Table 11.2: ASTM C177 k, capacity 1.67 kJ/kg K converted to J/kg K.
  (density 1160kg/m3) ; khk-plastic; typical; Table 11.1: ASTM D638 E, D792 density, nu no standard; Table 11.2: ASTM C177 k, capacity 1.67 kJ/kg K converted to J/kg K.
  (k 0.23) ; khk-plastic; typical; Table 11.1: ASTM D638 E, D792 density, nu no standard; Table 11.2: ASTM C177 k, capacity 1.67 kJ/kg K converted to J/kg K.
  (capacity 1670) ; khk-plastic; typical; Table 11.1: ASTM D638 E, D792 density, nu no standard; Table 11.2: ASTM C177 k, capacity 1.67 kJ/kg K converted to J/kg K.
))

; tpu: Polymaker PolyFlex TPU95 filament reference density. Small-strain E and other legacy properties remain unsourced; no as-printed strength is supplied.
(def tpu-polyflex95-reference (data
  (condition "Polymaker PolyFlex TPU95 filament reference density. Small-strain E and other legacy properties remain unsourced; no as-printed strength is supplied.")
  (evidence "std/material-sources.json#tpu-polyflex95-reference")
  (E 26) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (nu 0.48) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1200kg/m3) ; polymaker-tds-tpu95; typical; ISO 1183
  (k 0.19) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 1800) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def tpu tpu-polyflex95-reference) ; accepted compatibility alias

; polystyrene: Unfilled polystyrene, unspecified grade/process; legacy physical values are unsourced.
(def polystyrene-unfilled-unspecified (data
  (condition "Unfilled polystyrene, unspecified grade/process; legacy physical values are unsourced.")
  (evidence "std/material-sources.json#polystyrene-unfilled-unspecified")
  (E 3200) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (nu 0.34) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1050kg/m3) ; nist-thermal; typical; p19, solid polystyrene at 20 C; unspecified grade. Not foam.
  ; removed tensile: UNSOURCED; unsourced; Removed unsourced legacy tensile strength 45 MPa.
  (k 0.15) ; nist-thermal; typical; p19, solid polystyrene at 20 C; unspecified grade. Not foam.
  (capacity 1160) ; nist-thermal; typical; p19, solid polystyrene at 20 C; unspecified grade. Not foam.
))
(def polystyrene polystyrene-unfilled-unspecified) ; accepted compatibility alias

; aluminum-6061: 6061-T6 extrusion, Hydro specified mechanical limits, wall <=25 mm. Physical properties typical; capacity at 100 C, other physical values at 20–25 C.
(def aluminum-6061-t6-extrusion (data
  (condition "6061-T6 extrusion, Hydro specified mechanical limits, wall <=25 mm. Physical properties typical; capacity at 100 C, other physical values at 20–25 C.")
  (evidence "std/material-sources.json#aluminum-6061-t6-extrusion")
  (strength-kind "yield: specification minimum")
  (E 68900) ; nasa-6061; typical; Summary density 2.70 g/cm3; §7.21 nu; §7.222 modulus 10 million psi converted and rounded to 68.9 GPa.
  (nu 0.33) ; nasa-6061; typical; Summary density 2.70 g/cm3; §7.21 nu; §7.222 modulus 10 million psi converted and rounded to 68.9 GPa.
  (density 2700kg/m3) ; nasa-6061; typical; Summary density 2.70 g/cm3; §7.21 nu; §7.222 modulus 10 million psi converted and rounded to 68.9 GPa.
  (yield 240) ; hydro-6061; specification minimum; Extruded mechanical property limits, T6/T6511; minimum Rp0.2.
  (k 167) ; hydro-6061; typical; Typical thermal conductivity at 25 C.
  (capacity 896) ; kaiser-6061; typical; 6061 physical table, specific heat at 100 C; not a 20 C measurement.
))
(def aluminum-6061 aluminum-6061-t6-extrusion) ; accepted compatibility alias

; aluminum-6063: 6063-T5 extruded profiles, wall <=25 mm; lowest EN 755-2:2016 minimum across this wall range. Physical constants typical 20 C; T6 physical table used for E/density/capacity.
(def aluminum-6063-t5-extrusion (data
  (condition "6063-T5 extruded profiles, wall <=25 mm; lowest EN 755-2:2016 minimum across this wall range. Physical constants typical 20 C; T6 physical table used for E/density/capacity.")
  (evidence "std/material-sources.json#aluminum-6063-t5-extrusion")
  (strength-kind "yield: specification minimum")
  (E 69000) ; hydro-6063-en; typical; p2 typical physical constants at 20 C, T6 row; E/density/capacity used also for T5 as stated in condition.
  (nu 0.33) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 2700kg/m3) ; hydro-6063-en; typical; p2 typical physical constants at 20 C, T6 row; E/density/capacity used also for T5 as stated in condition.
  (yield 110) ; hydro-6063-en; specification minimum; p1 EN 755-2:2016 minimum Rp0.2, lowest of <=10 and >10–25 mm bands.
  (k 209) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 901) ; hydro-6063-en; typical; p2 typical physical constants at 20 C, T6 row; E/density/capacity used also for T5 as stated in condition.
))
(def aluminum-6063 aluminum-6063-t5-extrusion) ; accepted compatibility alias

; aluminum-6063-t6: 6063-T6 extruded profiles, wall <=25 mm; lowest EN 755-2:2016 minimum across this wall range. Physical constants typical 20 C; T6 physical table used for E/density/capacity.
(def aluminum-6063-t6-extrusion (data
  (condition "6063-T6 extruded profiles, wall <=25 mm; lowest EN 755-2:2016 minimum across this wall range. Physical constants typical 20 C; T6 physical table used for E/density/capacity.")
  (evidence "std/material-sources.json#aluminum-6063-t6-extrusion")
  (strength-kind "yield: specification minimum")
  (E 69000) ; hydro-6063-en; typical; p2 typical physical constants at 20 C, T6 row; E/density/capacity used also for T5 as stated in condition.
  (nu 0.33) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 2700kg/m3) ; hydro-6063-en; typical; p2 typical physical constants at 20 C, T6 row; E/density/capacity used also for T5 as stated in condition.
  (yield 160) ; hydro-6063-en; specification minimum; p1 EN 755-2:2016 minimum Rp0.2, lowest of <=10 and >10–25 mm bands.
  (k 201) ; hydro-6063-en; typical; p2 T6, typical at 20 C.
  (capacity 901) ; hydro-6063-en; typical; p2 typical physical constants at 20 C, T6 row; E/density/capacity used also for T5 as stated in condition.
))
(def aluminum-6063-t6 aluminum-6063-t6-extrusion) ; accepted compatibility alias

; aluminum-6005a-t5: 6005A-T5 extrusion, wall <=25 mm; Hydro North American limits; physical properties typical at 20–25 C.
(def aluminum-6005a-t5-extrusion (data
  (condition "6005A-T5 extrusion, wall <=25 mm; Hydro North American limits; physical properties typical at 20–25 C.")
  (evidence "std/material-sources.json#aluminum-6005a-t5-extrusion")
  (strength-kind "yield: specification minimum")
  (E 70000) ; hydro-6005a-en; typical; p2 typical physical constants at 20 C; T6 physical row used as reference for T5/T61.
  (nu 0.33) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 2710kg/m3) ; hydro-6005a-en; typical; p2 typical physical constants at 20 C; T6 physical row used as reference for T5/T61.
  (yield 215) ; hydro-6005a-na; specification minimum; Extruded mechanical property limits, <=25 mm.
  (k 193) ; hydro-6005a-na; typical; Extruded physical table, T5, 25 C.
  (capacity 898) ; hydro-6005a-en; typical; p2 typical physical constants at 20 C; T6 physical row used as reference for T5/T61.
))
(def aluminum-6005a-t5 aluminum-6005a-t5-extrusion) ; accepted compatibility alias

; aluminum-6005a-t61: 6005A-T61 extrusion, wall <=25 mm; Hydro North American limits; physical properties typical at 20–25 C.
(def aluminum-6005a-t61-extrusion (data
  (condition "6005A-T61 extrusion, wall <=25 mm; Hydro North American limits; physical properties typical at 20–25 C.")
  (evidence "std/material-sources.json#aluminum-6005a-t61-extrusion")
  (strength-kind "yield: specification minimum")
  (E 70000) ; hydro-6005a-en; typical; p2 typical physical constants at 20 C; T6 physical row used as reference for T5/T61.
  (nu 0.33) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 2710kg/m3) ; hydro-6005a-en; typical; p2 typical physical constants at 20 C; T6 physical row used as reference for T5/T61.
  (yield 240) ; hydro-6005a-na; specification minimum; Extruded mechanical property limits, <=25 mm.
  (k 188) ; hydro-6005a-na; typical; Extruded physical table, T61, 25 C.
  (capacity 898) ; hydro-6005a-en; typical; p2 typical physical constants at 20 C; T6 physical row used as reference for T5/T61.
))
(def aluminum-6005a-t61 aluminum-6005a-t61-extrusion) ; accepted compatibility alias

; aluminum-5052: 5052-H32 sheet, thickness >0.63 through 50 mm; BlueScope mechanical limits. Physical legacy values except where sourced remain unsourced.
(def aluminum-5052-h32-sheet (data
  (condition "5052-H32 sheet, thickness >0.63 through 50 mm; BlueScope mechanical limits. Physical legacy values except where sourced remain unsourced.")
  (evidence "std/material-sources.json#aluminum-5052-h32-sheet")
  (strength-kind "yield: specification minimum")
  (E 70300) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (nu 0.33) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 2680kg/m3) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (yield 160) ; bluescope-5052; specification minimum; p3 mechanical property specification, minimum yield.
  (k 138) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 880) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def aluminum-5052 aluminum-5052-h32-sheet) ; accepted compatibility alias

; aluminum-7075: Bare 7075-T6 sheet, 0.008–0.249 in thick. Yield is lowest minimum across that range; density/capacity typical alloy data; k is T651 plate reference, not measured T6 sheet.
(def aluminum-7075-t6-sheet (data
  (condition "Bare 7075-T6 sheet, 0.008–0.249 in thick. Yield is lowest minimum across that range; density/capacity typical alloy data; k is T651 plate reference, not measured T6 sheet.")
  (evidence "std/material-sources.json#aluminum-7075-t6-sheet")
  (strength-kind "yield: specification minimum")
  (E 71700) ; nasa-7075; typical; §7.21 nu; §7.222 typical modulus10.4 million psi converted and rounded to71.7 GPa.
  (nu 0.33) ; nasa-7075; typical; §7.21 nu; §7.222 typical modulus10.4 million psi converted and rounded to71.7 GPa.
  (density 2800kg/m3) ; kaiser-7075; typical; p2, density at 20 C, specific heat at 100 C, k for T651 at 20 C (explicit temper proxy).
  (yield 434) ; kaiser-sheet-plate; specification minimum; Table 7.2, 63 ksi minimum in thinnest band; 63*6.894757293=434.3697 MPa, conservatively rounded down.
  (k 130) ; kaiser-7075; typical; p2, density at 20 C, specific heat at 100 C, k for T651 at 20 C (explicit temper proxy).
  (capacity 960) ; kaiser-7075; typical; p2, density at 20 C, specific heat at 100 C, k for T651 at 20 C (explicit temper proxy).
))
(def aluminum-7075 aluminum-7075-t6-sheet) ; accepted compatibility alias

; steel: SSAB ASTM A36 hot rolled plate, thickness 2.54–152 mm. Typical physical references: maker E/nu/density and generic-steel handbook thermal values.
(def steel-a36-hot-rolled-plate (data
  (condition "SSAB ASTM A36 hot rolled plate, thickness 2.54–152 mm. Typical physical references: maker E/nu/density and generic-steel handbook thermal values.")
  (evidence "std/material-sources.json#steel-a36-hot-rolled-plate")
  (strength-kind "yield: specification minimum")
  (E 200000) ; ssab-questions; typical; Questions 3,10,17; approximate room-temperature values for SSAB steels, not a certificate for a particular A36 heat.
  (nu 0.3) ; ssab-questions; typical; Questions 3,10,17; approximate room-temperature values for SSAB steels, not a certificate for a particular A36 heat.
  (density 7850kg/m3) ; ssab-questions; typical; Questions 3,10,17; approximate room-temperature values for SSAB steels, not a certificate for a particular A36 heat.
  (yield 248) ; ssab-a36; specification minimum; 36 ksi required minimum; 36*6.894757293=248.2113 MPa, conservatively rounded down.
  (k 62) ; nist-thermal; typical; p16, generic steel thermal reference at 20 C; not a measured A36 heat. Density7850kg/m3.
  (capacity 480) ; nist-thermal; typical; p16, generic steel thermal reference at 20 C; not a measured A36 heat. Density7850kg/m3.
))
(def steel steel-a36-hot-rolled-plate) ; accepted compatibility alias

; stainless-304: Outokumpu Core 304/4301 hot rolled plate/sheet, EN 10088-2 limits in source Table 5. Not cold worked.
(def stainless-304-4301-hot-rolled (data
  (condition "Outokumpu Core 304/4301 hot rolled plate/sheet, EN 10088-2 limits in source Table 5. Not cold worked.")
  (evidence "std/material-sources.json#stainless-304-4301-hot-rolled")
  (strength-kind "yield: specification minimum")
  (E 200000) ; outokumpu-core; typical; Physical properties table, typical 20 C except k/capacity as tabulated.
  (nu 0.265) ; nasa-snp; typical; p229, SS304 room-temperature reference; not a certificate for Outokumpu product.
  (density 7900kg/m3) ; outokumpu-core; typical; Physical properties table, typical 20 C except k/capacity as tabulated.
  (yield 210) ; outokumpu-core; specification minimum; Table 5 hot rolled strip/plate minimum Rp0.2.
  (k 15) ; outokumpu-core; typical; Physical properties table, typical 20 C except k/capacity as tabulated.
  (capacity 500) ; outokumpu-core; typical; Physical properties table, typical 20 C except k/capacity as tabulated.
))
(def stainless-304 stainless-304-4301-hot-rolled) ; accepted compatibility alias

; copper: KME Cu-ETP annealed physical reference at 20 C. No tensile yield minimum supplied: removed legacy 70 MPa, and no universal safety factor is available.
(def copper-cu-etp-annealed (data
  (condition "KME Cu-ETP annealed physical reference at 20 C. No tensile yield minimum supplied: removed legacy 70 MPa, and no universal safety factor is available.")
  (evidence "std/material-sources.json#copper-cu-etp-annealed")
  (E 130000) ; kme-c110; typical; Physical properties at 20 C, annealed.
  (nu 0.34) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 8920kg/m3) ; kme-c110; typical; Physical properties at 20 C, annealed.
  ; removed yield: UNSOURCED; unsourced; Removed unsourced70 MPa. KME R220 gives <=140 MPa informational maximum, not a minimum.
  (k 394) ; kme-c110; typical; Physical properties at 20 C, annealed.
  (capacity 394) ; kme-c110; typical; Physical properties at 20 C, annealed.
))
(def copper copper-cu-etp-annealed) ; accepted compatibility alias

; brass: Wieland C26000 H01 quarter-hard rod, diameter <=2.5 in stocked range. Yield is minimum at 0.5% extension, not0.2% proof stress. Physical values are maker references; nu remains unsourced.
(def brass-c26000-h01-rod (data
  (condition "Wieland C26000 H01 quarter-hard rod, diameter <=2.5 in stocked range. Yield is minimum at 0.5% extension, not0.2% proof stress. Physical values are maker references; nu remains unsourced.")
  (evidence "std/material-sources.json#brass-c26000-h01-rod")
  (strength-kind "yield: specification minimum")
  (E 110317) ; wieland-c260; typical; Physical properties table,20 C, maker republishes CDA physical facts.
  (nu 0.34) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 8530kg/m3) ; wieland-c260; typical; Physical properties table,20 C, maker republishes CDA physical facts.
  (yield 140) ; wieland-c260; specification minimum; H01 rod >1 inch, lowest minimum across stocked diameter range <=2.5 in; ASTM B927/B927M-23.
  (k 121.2) ; wieland-c260; typical; Physical properties table,20 C, maker republishes CDA physical facts.
  (capacity 377.1) ; wieland-c260; typical; Physical properties table,20 C, maker republishes CDA physical facts.
))
(def brass brass-c26000-h01-rod) ; accepted compatibility alias

; cast-iron: KIPP EN-GJL-200 grey cast iron reference. Not an equivalence claim for FC200 or ASTM class30. Listed strength ranges are typical, not certified design allowables.
(def cast-iron-en-gjl-200 (data
  (condition "KIPP EN-GJL-200 grey cast iron reference. Not an equivalence claim for FC200 or ASTM class30. Listed strength ranges are typical, not certified design allowables.")
  (evidence "std/material-sources.json#cast-iron-en-gjl-200")
  (strength-kind "tensile: datasheet typical; compressive: datasheet typical")
  (E 88000) ; kipp-iron; typical; Density7.15g/cm3; E lower endpoint of 88–113 GPa range.
  (nu 0.26) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 7150kg/m3) ; kipp-iron; typical; Density7.15g/cm3; E lower endpoint of 88–113 GPa range.
  (tensile 200) ; kipp-iron; datasheet typical; Lower endpoint of 200–300 MPa maker range; no product/thickness-specific minimum.
  (compressive 720) ; kipp-iron; datasheet typical; Maker compressive strength; no guaranteed minimum.
  (k 50) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 460) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def cast-iron cast-iron-en-gjl-200) ; accepted compatibility alias

; tungsten-carbide: Goodfellow WC94-Co6 rod. Physical ranges are reference values; no guaranteed strength supplied.
(def wc94-co6-goodfellow-rod (data
  (condition "Goodfellow WC94-Co6 rod. Physical ranges are reference values; no guaranteed strength supplied.")
  (evidence "std/material-sources.json#wc94-co6-goodfellow-rod")
  (E 600000) ; goodfellow-wc6co; typical; Maker physical properties: nu upper endpoint0.21–0.22; k lower endpoint60–80 at 20 C; capacity lower endpoint200–480 at 25 C. Scalar endpoints are not uncertainty bounds propagated by ink.
  (nu 0.22) ; goodfellow-wc6co; typical; Maker physical properties: nu upper endpoint0.21–0.22; k lower endpoint60–80 at 20 C; capacity lower endpoint200–480 at 25 C. Scalar endpoints are not uncertainty bounds propagated by ink.
  (density 14950kg/m3) ; goodfellow-wc6co; typical; Maker physical properties: nu upper endpoint0.21–0.22; k lower endpoint60–80 at 20 C; capacity lower endpoint200–480 at 25 C. Scalar endpoints are not uncertainty bounds propagated by ink.
  ; removed compressive: UNSOURCED; unsourced; Removed unsourced legacy strength; no guaranteed value for this named rod.
  ; removed bending: UNSOURCED; unsourced; Removed unsourced legacy strength; no guaranteed value for this named rod.
  (k 60) ; goodfellow-wc6co; typical; Maker physical properties: nu upper endpoint0.21–0.22; k lower endpoint60–80 at 20 C; capacity lower endpoint200–480 at 25 C. Scalar endpoints are not uncertainty bounds propagated by ink.
  (capacity 200) ; goodfellow-wc6co; typical; Maker physical properties: nu upper endpoint0.21–0.22; k lower endpoint60–80 at 20 C; capacity lower endpoint200–480 at 25 C. Scalar endpoints are not uncertainty bounds propagated by ink.
))
(def tungsten-carbide wc94-co6-goodfellow-rod) ; accepted compatibility alias

; glass: Annealed soda-lime silicate glass, AGC SUNMAX reference per EN572. Bending is characteristic45 MPa, not tensile and not a design allowable.
(def glass-soda-lime-annealed-en572 (data
  (condition "Annealed soda-lime silicate glass, AGC SUNMAX reference per EN572. Bending is characteristic45 MPa, not tensile and not a design allowable.")
  (evidence "std/material-sources.json#glass-soda-lime-annealed-en572")
  (strength-kind "bending: datasheet")
  (E 70000) ; agc-soda-lime; typical; Physical and mechanical properties per EN572.
  (nu 0.2) ; agc-soda-lime; typical; Physical and mechanical properties per EN572.
  (density 2500kg/m3) ; agc-soda-lime; typical; Physical and mechanical properties per EN572.
  ; removed tensile: agc-soda-lime; datasheet; Misnamed: source45MPa is characteristic bending strength.
  (k 1) ; agc-soda-lime; typical; Physical and mechanical properties per EN572.
  (capacity 720) ; agc-soda-lime; typical; Physical and mechanical properties per EN572.
  (bending 45) ; agc-soda-lime; datasheet; EN572 characteristic bending strength, before design/surface/duration factors.
))
(def glass glass-soda-lime-annealed-en572) ; accepted compatibility alias

; borosilicate-glass: SCHOTT BOROFLOAT33 flat borosilicate glass. Bending25 MPa per DIN52292T1, no tensile limit.
(def glass-borofloat33-schott (data
  (condition "SCHOTT BOROFLOAT33 flat borosilicate glass. Bending25 MPa per DIN52292T1, no tensile limit.")
  (evidence "std/material-sources.json#glass-borofloat33-schott")
  (strength-kind "bending: datasheet")
  (E 64000) ; schott-borofloat-pdf; typical; E/nu at 20 C; k at 90 C; capacity averaged20–100 C.
  (nu 0.2) ; schott-borofloat-pdf; typical; E/nu at 20 C; k at 90 C; capacity averaged20–100 C.
  (density 2230kg/m3) ; schott-borofloat-pdf; typical; E/nu at 20 C; k at 90 C; capacity averaged20–100 C.
  ; removed tensile: UNSOURCED; unsourced; Removed45MPa borrowed from soda-lime bending; not borosilicate tensile.
  (k 1.2) ; schott-borofloat-pdf; typical; E/nu at 20 C; k at 90 C; capacity averaged20–100 C.
  (capacity 830) ; schott-borofloat-pdf; typical; E/nu at 20 C; k at 90 C; capacity averaged20–100 C.
  (bending 25) ; schott-borofloat-pdf; datasheet; DIN52292T1 bending strength; characteristic/reference, not a design allowable.
))
(def borosilicate-glass glass-borofloat33-schott) ; accepted compatibility alias

; water: NIST water, liquid293.15K at 1bar; reference thermophysical values.
(def water-liquid-20c-1bar (data
  (condition "NIST water, liquid293.15K at 1bar; reference thermophysical values.")
  (evidence "std/material-sources.json#water-liquid-20c-1bar")
  (viscosity 0.0010016) ; nist-water; typical; Query293.15K,1bar. Cp4.18406J/gK converted to4184.06J/kgK.
  (density 998.207kg/m3) ; nist-water; typical; Query293.15K,1bar. Cp4.18406J/gK converted to4184.06J/kgK.
  (conductivity 0.598012) ; nist-water; typical; Query293.15K,1bar. Cp4.18406J/gK converted to4184.06J/kgK.
  (capacity 4184.06) ; nist-water; typical; Query293.15K,1bar. Cp4.18406J/gK converted to4184.06J/kgK.
))
(def water water-liquid-20c-1bar) ; accepted compatibility alias

; air: Legacy air reference, unspecified humidity/pressure; all four values unsourced. Supply measured/validated data for flow.
(def air-unspecified-reference (data
  (condition "Legacy air reference, unspecified humidity/pressure; all four values unsourced. Supply measured/validated data for flow.")
  (evidence "std/material-sources.json#air-unspecified-reference")
  (viscosity 1.81e-5) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 1.204kg/m3) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (conductivity 0.026) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 1005) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def air air-unspecified-reference) ; accepted compatibility alias

; oil: Legacy unspecified oil reference; no viscosity grade, maker or temperature. All four values unsourced; not a named real oil.
(def oil-unspecified-reference (data
  (condition "Legacy unspecified oil reference; no viscosity grade, maker or temperature. All four values unsourced; not a named real oil.")
  (evidence "std/material-sources.json#oil-unspecified-reference")
  (viscosity 0.1) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (density 880kg/m3) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (conductivity 0.13) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
  (capacity 1900) ; UNSOURCED; unsourced; Legacy physical reference value; no qualifying source located.
))
(def oil oil-unspecified-reference) ; accepted compatibility alias

(def other-names (data
  (aluminium-6061 aluminum-6061)
  (al6061 aluminum-6061)
  (aluminium-7075 aluminum-7075)
  (al7075 aluminum-7075)
  (ss304 stainless-304)
  (pa12 nylon-pa12)
  (nylon nylon-pa12)
  (pc polycarbonate)
  (aluminium aluminum-6061)
  (aluminum aluminum-6061)
  (aluminium-6063 aluminum-6063)
  (aluminium-5052 aluminum-5052)
  (stainless stainless-304)
  (stainless-steel stainless-304)
  (pp polypropylene)
  (pa66 nylon-pa66)
  (nylon-66 nylon-pa66)
  (acetal pom)
  (delrin pom)
  (ps polystyrene)
  (borosilicate borosilicate-glass)
  (soda-lime-glass glass)
  (carbide tungsten-carbide)))

; Source bibliography (individual facts; source documents are not vendored).
; prusa-pla: Prusa Polymers, Technical datasheet: Prusament PLA by Prusa Polymers. v1.1, 16-02-2022.
; https://prusament.com/wp-content/uploads/2022/10/PLA_Prusament_TDS_2021_10_EN.pdf
; prusa-petg: Prusa Polymers, Technical datasheet: Prusament PETG by Prusa Polymers. v1.1, 16-02-2022.
; https://storage.googleapis.com/prusa3d-content-prod-14e8-wordpress-prusament-prod/2023/10/9f8d2165-tds_prusament-petg_n_en.pdf
; ensinger-manual: Ensinger, Engineering plastics — The Manual. 10/12, E9911075A011GB.
; https://www.ensingerplastics.com/en-us/-/media/ensinger/files/document-teaser-files/brochures/shapes/engineering-plastics-manual-english.ashx
; ensinger-pet: Ensinger, TECAPET white material data. Undated live maker page; retrieved 2026-11-12.
; https://www.ensingerplastics.com/en/shapes/tecapet-white
; ensinger-peek: Ensinger, TECAPEEK natural material data. Undated live maker page; retrieved 2026-11-12.
; https://www.ensingerplastics.com/en/shapes/peek-tecapeek-natural
; ensinger-pp: Ensinger, TECAFINE PP natural material data. Undated live maker page; retrieved 2026-11-12.
; https://www.ensingerplastics.com/en-gb/shapes/tecafine-pp-natural
; ensinger-hdpe: Ensinger, TECAFINE PE 300 natural material data. Undated live maker page; retrieved 2026-11-12.
; https://www.ensingerplastics.com/en-sg/shapes/polyethylene-tecafine-pe-300-natural
; ensinger-ptfe: Ensinger, TECAFLON PTFE natural material data. Undated live maker page; retrieved 2026-11-12.
; https://www.ensingerplastics.com/en-gb/shapes/tecaflon-ptfe-natural
; khk-plastic: KHK, Gear Technical Reference: Design of Plastic Gears. Undated live page; retrieved 2026-11-12.
; https://khkgears.net/new/gear_knowledge/gear_technical_reference/design-of-plastic-gears.html
; polymaker-tds-tpu95: Polymaker, PolyFlex TPU95 Technical Data Sheet. v5.5, 2025-12-29.
; https://polymaker.com/wp-content/uploads/TDS_Polymaker_PolyFlex-TPU95_V5.5_2025_12_29_EN.pdf
; hydro-6061: Hydro, Alloy 6061. 2019/01 (2/2019).
; https://www.hydro.com/Document/Index?id=560718&name=Alloy+6061.pdf
; nasa-6061: NASA / R. F. Muraca and J. S. Whittick, Materials data handbook: Aluminum alloy 6061. NASA CR-123772, 2nd edition, May 1972.
; https://ntrs.nasa.gov/api/citations/19720022808/downloads/19720022808.pdf
; hydro-6005a-en: Hydro, Werkstoffdatenblatt: EN AW-6005A, Stranggepresste Profile. Rev.4, 03/2023.
; https://www.hydro.com/globalassets/08-about-hydro/hydro-worldwide/austria/nenzing/datenblatter-legierungen/hydro-en-aw-6005a/
; hydro-6005a-na: Hydro, Alloy 6005A. 2019/02 (3/2019).
; https://www.hydro.com/globalassets/01-products--services/extruded-profiles/americas/ena-resources/alloy-data-sheets/hydro_2019_data_sheet_6005a_revise_1.pdf
; hydro-6063-en: Hydro, Werkstoffdatenblatt: EN AW-6063, Stranggepresste Profile. Rev.4, 03/2023.
; https://www.hydro.com/globalassets/08-about-hydro/hydro-worldwide/austria/nenzing/datenblatter-legierungen/hydro-en-aw-6063/
; bluescope-5052: BlueScope Distribution, Aluminium 5052 Data Sheet. Undated; October 2022 URL edition.
; https://www.bluescopedistribution.com.au/wp-content/uploads/2022/10/Aluminium-5052-Data-Sheet-BlueScope-Distribution.pdf
; kaiser-sheet-plate: Kaiser Aluminum, Aluminum Sheet & Plate. KA-SP-AA5/6-1.09, January 2009.
; https://online.kaiseraluminum.com/depot/PublicProductInformation/Document/1010/Kaiser_Aluminum_Sheet___Plate.pdf
; kaiser-7075: Kaiser Aluminum, Sheet Coil & Plate Alloy 7075 — Technical Data. Rev.05/06.
; https://online.kaiseraluminum.com/depot/PublicProductInformation/Document/1017/Kaiser_Aluminum_7075_Sheet_Coil_and_Plate.pdf
; ssab-a36: SSAB, ASTM A36 — Mechanical Properties. Undated live maker page; retrieved 2026-11-12.
; https://www.ssab.com/en-us/brands-and-products/commercial-steel/structural-steel/astm-a36
; outokumpu-core: Outokumpu, Core range datasheet. 1560EN:4, November2022.
; https://www.outokumpu.com/-/media/files/products/core/outokumpu-core-range-datasheet.pdf?modified=20251117111909&revision=025e9931-a1d5-4c8f-8ff5-f881d38916da
; kme-c110: KME, Cu-ETP: Material data sheet. 01/2021.
; https://www.kme.com/fileadmin/DOWNLOADCENTER/COPPER%20DIVISON/4%20Industrial%20Rolled/3_Datasheets/Datasheets_NEW_2021/Cu-ETP_01_2021_e.pdf
; wieland-c260: Wieland Concast, Alloy C26000: 70/30 Brass. Mechanical properties per ASTM B927/B927M-23.
; https://www.concast.com/files/C26000.pdf
; kipp-iron: KIPP, EN-GJL-200 material properties. Undated live maker page; retrieved 2026-11-12.
; https://www.kippusa.com/en-us/en-gjl-200
; goodfellow-wc6co: Goodfellow, Tungsten Carbide-Cobalt Rod | WC94-Co6. Undated live maker page; retrieved 2026-11-12.
; https://www.goodfellow.com/usa/tungsten-carbide-cobalt-rod-group
; agc-soda-lime: AGC Glass Europe, SUNMAX Premium: Product factsheet. 08/2020.
; https://www.agc-yourglass.com/sites/default/files/agc_docs/Factsheet_Sunmax_0820_EN.pdf
; schott-borofloat-pdf: SCHOTT, BOROFLOAT 33 — Technical Data. Undated PDF mirrored by Edmund Optics; retrieved 2026-11-12.
; https://www.edmundoptics.com/ViewDocument/schott-borofloat-technical-data-sheet-english.pdf
; nist-water: NIST, NIST Chemistry WebBook SRD69: water, thermophysical properties. Live query293.15K,1bar; retrieved2026-11-12.
; https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C7732185&Type=IsoBar&Digits=6&P=1&TLow=293.15&THigh=293.15&TInc=1&RefState=DEF&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm
; ssab-questions: SSAB, 20 questions about steel. Undated live maker reference; retrieved 2026-11-12.
; https://www.ssab.com/en/support/product-material-data/steel/20-questions
; nasa-7075: NASA / R. F. Muraca and J. S. Whittick, Materials data handbook: Aluminum alloy7075. NASA CR-123773, 2nd edition, May1972.
; https://ntrs.nasa.gov/api/citations/19720022809/downloads/19720022809.pdf
; nasa-snp: NASA, Space Nuclear Propulsion Material Property Handbook. SNP-HDBK-0008, revision3.1, 08/05/2024.
; https://ntrs.nasa.gov/api/citations/20240004217/downloads/SNP-HDBK-0008_SNP-Material-Handbook.pdf?attachment=true
; nist-thermal: NBS / Daniel Gross, Data Sources for Parameters Used in Predictive Modeling of Fire Growth and Smoke Spread. NBSIR85-3223, September1985.
; https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=101567
; kaiser-6061: Kaiser Aluminum, Sheet Coil & Plate Alloy 6061 — Technical Data. Rev. 05/06.
; https://online.kaiseraluminum.com/depot/PublicProductInformation/Document/1015/Kaiser_Aluminum_5052_Sheet_Coil_and_Plate.pdf

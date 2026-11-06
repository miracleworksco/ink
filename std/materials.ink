; Materials and fluids, as plain IR data. A body's (attr material "steel")
; names one, or is data with their keys; ink mass reads its density and
; packages the rest. Free text goes in (attr material-note "..."). Values
; are typical at room temperature, in the IR's units (MPa, g/mm^3, W/(m K),
; J/(kg K)), densities written in kg/m3. Datasheets vary by grade, process
; and print settings: where it matters, define your own, (def my-alloy
; (data (E ...) ...)).
;
; Keys: E (Young's modulus), nu (Poisson's ratio), density, yield, k
; (conductivity), capacity (specific heat).
;
; yield is a ductile yield stress. Brittle and grained materials have none
; (no one limit fits concrete in tension or wood across the grain); they
; carry named strengths (MPa) from a named grade.


; Printed polymers, as printed, in-plane. Parts are weaker across layers.
(def pla (data (E 3500) (nu 0.36) (density 1240kg/m3) (yield 50) (k 0.13) (capacity 1800)))
(def petg (data (E 2100) (nu 0.38) (density 1270kg/m3) (yield 45) (k 0.2) (capacity 1200)))
(def abs (data (E 2200) (nu 0.35) (density 1050kg/m3) (yield 40) (k 0.17) (capacity 1400)))
(def nylon-pa12 (data (E 1700) (nu 0.39) (density 1010kg/m3) (yield 45) (k 0.25) (capacity 1700)))
(def polycarbonate (data (E 2300) (nu 0.37) (density 1200kg/m3) (yield 60) (k 0.2) (capacity 1250)))

; Moulded and extruded polymers: unfilled, room temperature, dry.
(def polypropylene (data (E 1500) (nu 0.42) (density 905kg/m3) (yield 33) (k 0.2) (capacity 1920)))
(def hdpe (data (E 1000) (nu 0.42) (density 955kg/m3) (yield 26) (k 0.46) (capacity 1900)))
(def pet (data (E 2800) (nu 0.37) (density 1380kg/m3) (yield 55) (k 0.2) (capacity 1200)))
(def nylon-pa66 (data (E 2900) (nu 0.39) (density 1140kg/m3) (yield 80) (k 0.25) (capacity 1670)))
(def pom (data (E 3000) (nu 0.35) (density 1410kg/m3) (yield 65) (k 0.31) (capacity 1470)))
(def ptfe (data (E 500) (nu 0.46) (density 2170kg/m3) (yield 10) (k 0.25) (capacity 1000)))
(def peek (data (E 3600) (nu 0.4) (density 1300kg/m3) (yield 100) (k 0.25) (capacity 1340)))
; Cast PA6, MC901 (KHK's nylon gears): KHK, Design of Plastic Gears, Tables
; 11.1-11.2, dry. E is the tensile modulus. KHK gives no tensile yield
; (tensile strength 96, compressive yield 103), so none here.
(def nylon-pa6-mc901 (data (E 3432) (nu 0.4) (density 1160kg/m3) (k 0.23) (capacity 1670)))
; Thermoplastic polyurethane, Shore 95A: an elastomer, E at small strain,
; no yield.
(def tpu (data (E 26) (nu 0.48) (density 1210kg/m3) (k 0.19) (capacity 1800)))
; General-purpose polystyrene: brittle, so a tensile strength, no yield.
(def polystyrene (data (E 3200) (nu 0.34) (density 1050kg/m3) (tensile 45) (k 0.14) (capacity 1300)))

; Metals.
(def aluminum-6061 (data (E 68900) (nu 0.33) (density 2700kg/m3) (yield 276) (k 167) (capacity 896)))
; 6063 as T5 (extrusions), 5052 as H32 (sheet).
(def aluminum-6063 (data (E 68900) (nu 0.33) (density 2700kg/m3) (yield 145) (k 209) (capacity 900)))
; 6063-T6 (80/20's T-slot profiles): ASM MatWeb "Aluminum 6063-T6", the
; Aluminum Association's typical values.
(def aluminum-6063-t6 (data (E 68900) (nu 0.33) (density 2700kg/m3) (yield 214) (k 200) (capacity 900)))
; 6005A extrusions: T5 (MISUMI's HFS frames), T61 (AndyMark's hex stock).
; E, density, capacity: Hydro's EN AW-6005A data sheet, typical at 20 C.
; k by temper, and yield: its North American 6005A sheet; yield is the
; minimum up to 6.3 mm (no typical value found). nu: MakeItFrom 6005A-T5.
(def aluminum-6005a-t5 (data (E 70000) (nu 0.33) (density 2710kg/m3) (yield 215) (k 193) (capacity 898)))
(def aluminum-6005a-t61 (data (E 70000) (nu 0.33) (density 2710kg/m3) (yield 240) (k 188) (capacity 898)))
(def aluminum-5052 (data (E 70300) (nu 0.33) (density 2680kg/m3) (yield 193) (k 138) (capacity 880)))
(def aluminum-7075 (data (E 71700) (nu 0.33) (density 2810kg/m3) (yield 503) (k 130) (capacity 960)))
(def steel (data (E 200000) (nu 0.3) (density 7850kg/m3) (yield 250) (k 50) (capacity 490)))
(def stainless-304 (data (E 193000) (nu 0.29) (density 8000kg/m3) (yield 215) (k 16.2) (capacity 500)))
(def copper (data (E 117000) (nu 0.34) (density 8960kg/m3) (yield 70) (k 398) (capacity 385)))
(def brass (data (E 100000) (nu 0.34) (density 8500kg/m3) (yield 200) (k 120) (capacity 380)))

; Others: brittle or grained, so no yield.
; Grey cast iron, FC200 / ASTM class 30: tensile is the minimum tensile
; strength, compressive typical.
(def cast-iron (data (E 100000) (nu 0.26) (density 7200kg/m3) (tensile 200) (compressive 750) (k 50) (capacity 460)))
; Cemented tungsten carbide, WC-6%Co (cutting tools): bending is the
; transverse rupture strength.
(def tungsten-carbide (data (E 600000) (nu 0.22) (density 14900kg/m3) (compressive 5000) (bending 3000) (k 85) (capacity 210)))
; Annealed soda-lime float glass (EN 572-1): tensile is the characteristic
; bending strength, short-term, before any design factor.
(def glass (data (E 70000) (nu 0.22) (density 2500kg/m3) (tensile 45) (k 1) (capacity 840)))
; Borosilicate 3.3 (ISO 3585; labware): tensile as for soda-lime glass.
(def borosilicate-glass (data (E 64000) (nu 0.2) (density 2230kg/m3) (tensile 45) (k 1.2) (capacity 830)))

; Each material has one name, above. A model writing another is refused
; with the name to write: the name in another case or with spaces for
; hyphens, or an other-names entry, (other-name material). A generic name
; is not a material: its hint is its commonest grade, to write where the
; source names none. other-names is data, not a material.
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

; Fluids, for flow packages. Keys: viscosity (Pa s), density, conductivity,
; capacity.
(def water (data (viscosity 1.0e-3) (density 998kg/m3) (conductivity 0.6) (capacity 4186)))
(def air (data (viscosity 1.81e-5) (density 1.204kg/m3) (conductivity 0.026) (capacity 1005)))
(def oil (data (viscosity 0.1) (density 880kg/m3) (conductivity 0.13) (capacity 1900)))

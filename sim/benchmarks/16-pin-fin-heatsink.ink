; Benchmark 16: pin-fin heat sink vs fin theory.
;
; Task: an aluminium (k 167) pin-fin heat sink, base 40 x 40 x 4, 6 x 6 pins
; d = 3, L = 20 on a 6.5 mm pitch. A 10 W chip heats a r = 8 patch under
; the base; every surface loses heat to 25 C air at h = 25 W/(m^2 K).
;
; Hand calc (isothermal base, fins with the corrected-length tip):
;   m = sqrt(4 h / (k d)) = sqrt(4 25 / (167 0.003)) = 14.13 1/m
;   Lc = L + d/4 = 20.75 mm, m Lc = 0.2932, eta = tanh(m Lc) / (m Lc) = 0.9724
;   fins: 36 pi d Lc eta = 36 x 195.56 x 0.9724 = 6846 mm^2 effective
;   bare base: top 1600 - 36 pi 1.5^2 = 1345.5, bottom 1600, sides 640
;   h A_eff = 25 W/(m^2 K) x (6846 + 3585.5) mm^2 = 0.2608 W/K
;   base temperature rise 10 / 0.2608 = 38.35 K -> 63.35 C
; The chip, the hottest point, sits above the base by the conduction across
; it and the spreading from the r = 8 patch (Lee, Song, Au and Moran's
; constriction model, 1995: a source of radius 8 on a disc of the base's
; area, b = 22.57, 4 thick, h_eff = 0.2608 W/K / 1600 mm^2 = 163 W/(m^2 K)):
;   0.15 K + 2.24 K, so about 40.7 K at the chip (6% over 38.35).
;
; Pass criteria (ink solve):
;   - max temperature rise within 10% of 38.35 K;
;   - heat balance: convection removes 10 W (within 0.5%).

(import "std")

(param pins 6)          ; per side (load-time)
(param pitch 6.5)
(param d 3)
(param L 20)
(def base (move (box 40 40 4) 0 0 2))
(def pin (move (cylinder-z (* 0.5 d) L) 0 0 (+ 4 (* 0.5 L))))
(def row0 (* -0.5 (* pitch (- pins 1))))
(def pin-array (pattern-linear (pattern-linear (move pin row0 row0 0) pins pitch 0 0) pins 0 pitch 0))
(body heatsink (shape (min base pin-array)) (attr material "aluminum-6061") (attr color "#B8C2CC"))

(study cooling (shape (shape-of heatsink)) (heat (k 167W/m/K))
  (source chip (power 10W) (region (max (- (len2 x y) 8) (- z 1))))
  (convection air (face heatsink) (coefficient 25W/m2/K) (ambient 25degC))
  (res 0.5))

(require chip-cool (max-temperature cooling < 70degC))
; the rise within 10% of fin theory's 38.35 K over the 25 C air (above)
(require rise-lo (max-temperature cooling > 59.515degC))
(require rise-hi (max-temperature cooling < 67.185degC))

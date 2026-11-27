; A 2 W aluminium heater plate on a printed PETG block: each body conducts
; by its own material, and the block's temperature is claimed alone.
(import "std")
(body block (shape (move (box 40 40 8) 0 0 4)) (attr material "petg"))
(body plate (shape (move (box 32 32 8) 0 0 12)) (attr material "aluminum-6061"))
(assembly root (instance block block) (instance plate plate))

(study warm (shape (shape-of root))
  (heat (k (* 1W/m/K (by-body root k 0.2))))    ; each body's k, else 0.2, given its unit whole
  (source heater (power 2W) (region (shape-of root/plate)))
  (convection air (- 0.5 z) (coefficient 10W/m2/K) (ambient 25degC))   ; all but the bottom
  (res 1))

(require hottest-under-120 (max-temperature warm < 120degC))
; One body's temperature: the block lies where the field is under 100 (C).
(require block-under-100 (subset (shape-of root/block) (- (study-field warm temperature) 100)))

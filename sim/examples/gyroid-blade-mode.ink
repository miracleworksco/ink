; A tapered aluminium blade made of an open gyroid lattice on a solid root,
; vibrating in one of its natural modes: the modes are solved on the lattice
; itself (no mesh), the mode's shape and stress come back as fields, and the
; blade drawn is the lattice warped by its shape times a scale that swings
; with t (Space plays it; one period a second, the swing exaggerated to be
; seen: a viewer states the scale). Coloured by how far each point swings: the
; still band partway along is the mode's node.
(import "std")

(param mode 3)         ; which mode is played (from 1)
(param swing 6)        ; the shape's peak shown, mm (the mode is 1 mm at its peak)

; The blade: 120 mm long, tapering from 24 to 12 mm thick, on a solid root.
(def taper (+ 1 (* -0.5 (/ (+ x 60) 120))))      ; 1 at the root, 0.5 at the tip
(def envelope (max (- (abs x) 60) (- (abs y) (* 16 taper)) (- (abs z) (* 12 taper))))
(def lattice (lattice-fill envelope (gyroid 12) 1.8))   ; walls at least 1.8 mm thick
(def root (move (box 14 40 32) -67 0 0))
(def blade (fillet-union root lattice 2 90))

(study vib (shape blade) (elastic (E 69000MPa) (nu 0.33) (density 2.7g/cm3))
  (fix base (+ x 72))
  (modes 4)
  (res 0.8))

(t-span 0 1s)
(def k (* swing (sin (* 360deg/s t))))
(def ux (study-field vib displacement-x (mode mode)))
(def uy (study-field vib displacement-y (mode mode)))
(def uz (study-field vib displacement-z (mode mode)))
(body blade-mode
  (shape (warp blade (- x (* k ux)) (- y (* k uy)) (- z (* k uz))))
  ; How far each point swings (the mode's shape, 0 at the root to 1 mm at
  ; its peak), the root bare metal.
  (attr color (mix (vec 0.6 0.62 0.66) (color-map (study-field vib displacement (mode mode)) 0 1)
                   (clamp (* 8 (study-field vib displacement (mode mode))) 0 1)))
  (attr metallic 0.5)
  (attr roughness 0.35))

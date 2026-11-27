; A domed speaker grille, 172 mm across, perforated by 2,713 holes
; on a hex grid, their size a halftone: 1.0 mm in the open, 1.9 mm inside
; the letters of a logo (or, with --param art=1, the ripples of two
; sources), shrinking toward the rim. The holes are one field,
; a hex lattice whose hole radius reads a field of position, so the logo is
; drawn by the perforation itself, and the count costs nothing: the grid
; is never expanded into features.
(import "std")

(param pitch 2.6)      ; hole spacing, mm (hex)
(param skin 1.5)       ; dome thickness, mm
(param dome-r 76)      ; perforated dome's radius at its foot, mm
(param sphere 300)     ; the dome's spherical radius, mm
(param art 0)          ; the halftone: 0 the logo, 1 two sources' ripples

(def rho (len2 x y))

; The dome: a spherical shell meeting z = 0 at dome-r, its sag 9.8 mm.
(def foot (sqrt (- (square sphere) (square dome-r))))
(def dome (max (- (abs (- (len3 x y (+ z foot)) (- sphere (* 0.5 skin)))) (* 0.5 skin))
               (- rho dome-r) (neg (+ z 2))))

; The mounting flange: a ring 74 to 86 mm in radius, 2 mm thick, with four
; 4.5 mm screw holes on an 81 mm radius, blended smoothly into the dome.
(def flange (max (- (abs (- rho 80)) 6) (- (abs (+ z 1)) 1)))
(def screws (pattern-polar (move (cylinder-z 2.25 20) (* 81 (cos (/ pi 4))) (* 81 (sin (/ pi 4))) 0) 4))
(def body (subtract (smooth-union dome flange 3) screws))

; The halftone: k is 1 inside the logo's letters and 0 outside, graded
; over 3 mm; near the rim the holes shrink to 60 %. The logo is "ink" in
; std's sans at its boldest, 44 mm to the caps, centred on the dome.
(def logo (lettering-2d (lettering-styled "ink" (font "sans" (wght 900)) 44 3 1 0.5)))
; Or the ripples of two sources 50 mm apart, 16 mm in wavelength: k is
; their summed wave, 0 to 1 (a stereo pair's interference pattern).
(fn ripple (sx) (cos (* (/ (* 2 pi) 16) (len2 (- x sx) y))))
(def waves (clamp (+ 0.5 (* 0.4 (+ (ripple -25) (ripple 25)))) 0 1))
(def logo-k (clamp (- 0.5 (/ logo 3)) 0 1))

; The grille whose halftone is k (0 to 1), its holes from radius lo (k 0)
; to 0.95 mm (k 1).
(def py (* pitch (sqrt 3)))
(fn grille-with (k lo)
  (let ((hole-r (* (mix lo 0.95 k) (clamp (/ (- 72 rho) 10) 0.6 1)))
        ; A hex lattice is two rectangular ones, the second offset by half a
        ; cell: the nearest hole is the nearer of the two lattices' nearest.
        ; Holes are whole: one is drilled only where its centre is within
        ; 71 mm of the axis.
        (lattice-hole (lambda (ox oy)
          (let ((cx (- (* pitch (floor (+ 0.5 (/ (+ x ox) pitch)))) ox))
                (cy (- (* py (floor (+ 0.5 (/ (+ y oy) py)))) oy)))
            (select (- (len2 cx cy) 71) (- (len2 (- x cx) (- y cy)) hole-r) 10))))
        (holes (max (min (lattice-hole 0 0) (lattice-hole (* 0.5 pitch) (* 0.5 py)))
                    (- (abs z) 20))))
    (max (subtract body holes) (cylinder-z 90 30))))

(body grille
  (shape (if (== art 0) (grille-with logo-k 0.5) (grille-with waves 0.3)))
  (attr material "aluminum-5052")
  (attr color (vec 0.09 0.09 0.10))
  (attr metallic 0.4)
  (attr roughness 0.5))

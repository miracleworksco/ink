; Benchmark 03: lofted, twisted blade with an elastic study.
;
; Task: root NACA 0018-like section, chord 60; tip NACA 0010-like, chord 30;
; span 150 mm, 30 deg of linear twist about the quarter-chord line (the z
; axis), chord linear from root to tip. Aluminium; root clamped, 50 N tip load
; normal to the root chord (+y).
;
; Reference: Euler-Bernoulli on the same field (tools/blade_beam.py): tip
; deflection 0.945 mm, peak stress 29.0 MPa.
;
; expect-status: strong-enough unknown
; Pass criteria (ink require):
;   - tip deflection under 5 mm (holds numerically);
;   - the blade within its envelope; safety factor above 1.5.

(import "std")

(param span 150)
(param c-root 60)           ; chord (mm)
(param t-root 0.18)         ; thickness / chord
(param c-tip 30)
(param t-tip 0.10)
(param twist 30)            ; deg, root to tip
(param t-te 0.01)           ; trailing-edge thickness / chord: 0.6 mm at the root, 0.3 at the tip
; tip load 50 N: study forces must be literal numbers ("expected a constant
; number" for a param), so the load cannot be a design parameter.

; Symmetric NACA 4-digit section in (x, y): chord c along +x, quarter chord at
; the origin, thickness ratio th. The closed section (-0.1036 coefficient)
; with its half-thickness raised linearly from the leading edge to t-te / 2
; at the trailing edge: a finite trailing edge 1% of the chord thick, as
; machined and cast blades have (about 0.5 to 1%).
(fn naca-2d (c th)
  (let ((u (+ (/ x c) 0.25))                        ; 0 at LE, 1 at TE
        (uc (clamp u 0 1))
        (yt (+ (* 5 th c (+ (* 0.2969 (sqrt uc)) (* -0.1260 uc) (* -0.3516 uc uc)
                            (* 0.2843 uc uc uc) (* -0.1036 uc uc uc uc)))
               (* 0.5 t-te c uc))))
    (max (- (abs y) yt) (* c (neg u)) (* c (- u 1)))))

(def s (clamp (/ z span) 0 1))                      ; 0 at root, 1 at tip
(def chord (mix c-root c-tip s))

; Rotate a 2D-section field by the local twist angle (z-dependent rotation).
(fn twisted (g)
  (let ((a (deg->rad (* twist s))))
    (warp g (+ (* (cos a) x) (* (sin a) y)) (- (* (cos a) y) (* (sin a) x)) z)))

; Naive morph of the real-size sections with std `loft-z` (kept for
; comparison; not used): mid-span chord 37.0 instead of 45.0.
(def blade-naive (twisted (loft-z (naca-2d c-root t-root) (naca-2d c-tip t-tip) 0 span)))

; Normalised morph: unit-chord sections, morphed, then scaled by chord(s).
(def blade-field
  (max (twisted (* chord (warp (blend (naca-2d 1 t-root) (naca-2d 1 t-tip) s)
                               (/ x chord) (/ y chord) z)))
       (slab-z 0 span)))

(body blade
  (shape blade-field)
  (attr color "#9DA7B3")
  (attr material "aluminum-6061") (attr material-note "6061-T6")
  (attr density 2.7g/cm3))

(study tip-load
  (shape blade-field)
  (elastic (E 69000MPa) (nu 0.33) (density 2700kg/m3) (yield 240MPa))
  (fix root (- z 2))                                 ; bottom 2 mm clamped
  (load tip (force 0 50N 0) (region (- (- span 2) z)))   ; top 2 mm
  (res 0.9))

; The blade's faces lie on its envelope's: no search separates the faces at
; the model's floor within its budget, so it is decided at cells of 0.05 (to
; 0.1 mm).
(require in-envelope (within (shape-of blade) (-20 -20 0) (45 20 150)) (res 0.05))
; With its trailing edge 1% of the chord thick (0.6 mm at the root, 0.3 at
; the tip: two or more of the solve's finest cells), every grid resolves
; it. The tip, 0.877 mm, is decided within its grid error of the coarser
; grids'. The factor, 8.751, rose 8.655 and 8.691 on the coarser two: grids
; moving one way bound nothing, so it is not settled here; on cells of 0.45
; (7 minutes) the three converge and it holds, 8.760. (Closed, the edge was a
; wedge thinner than those cells for 0.5 mm back from it, which they may
; hold too stiff: neither claim was decided. The beam reference's 0.945 mm
; is a closed section's, a little softer at its tip.)
(require strong-enough (min-safety-factor tip-load > 1.5))
(require tip-deflection (max-displacement tip-load < 5))

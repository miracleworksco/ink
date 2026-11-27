; Benchmark 08: helical compression spring, stiffness vs hand calc.
;
; Task: a steel compression spring (wire d = 3, mean coil diameter D = 20,
; n = 4 active turns, pitch 8) between two 5 mm end plates. Fix the bottom
; plate, guide the top plate (x, y held, z free) and push it down with 50 N.
; The spring is a tube swept along a helix polyline (24 segments per turn).
;
; Hand calc (steel, E = 200 GPa, nu = 0.3, G = E / 2(1+nu) = 76.9 GPa):
;   k = G d^4 / (8 D^3 n) = 76923 * 81 / (8 * 8000 * 4) = 24.34 N/mm
;   deflection under 50 N: 2.054 mm
;   Wahl shear stress: C = D/d = 6.67, Kw = (4C-1)/(4C-4) + 0.615/C = 1.225,
;   tau = Kw 8 F D / (pi d^3) = 115.5 MPa; von Mises = sqrt(3) tau = 200 MPa.
; End correction: the wire is fused to a plate until its centre is d/2
; above it, i.e. for 1.5/8 = 0.19 turn at each end, so n_active = 3.625 and
; k = 26.86 N/mm, deflection 1.861 mm. (Real closed-and-ground springs lose
; about one turn per end the same way.) Where the wire leaves a plate it
; rises at the pitch's 7 degrees, so wire and plate would meet in a wedge
; whose tip is thinner than any cell; a real end sits flat on its plate (or
; is welded to it), so the wire is joined to the plates by a smooth union
; of width 1 mm, a fillet of about 0.85 mm where they meet square, which
; fills the wedge's tip and leaves the coil as it is.
;
; Pass criteria (ink require, then the cached ink solve):
;   - top plate deflection within 10% of 1.861 mm (1.675 .. 2.047);
;   - the 99th-percentile von Mises stress within 25% of 200 MPa (the peak
;     sits at the wire/plate junction, a sharp re-entrant corner).

(import "std")

(param d-wire 3)
(param D-coil 20)
(param turns 4)      ; load-time: it sets the helix point count
(param pitch 8)
(param plate-t 5)
(param force 50)

(def coil-pts (helix-points (* 0.5 D-coil) pitch turns (* 24 turns)))
(def wire (tube coil-pts (* 0.5 d-wire)))
(def free-length (* pitch turns))
(def bottom-plate (move (cylinder-z 13 plate-t) 0 0 (* -0.5 plate-t)))
(def top-plate (move (cylinder-z 13 plate-t) 0 0 (+ free-length (* 0.5 plate-t))))
(body spring (shape (smooth-union wire (min bottom-plate top-plate) 1)) (attr material "steel") (attr color "#707880"))

(study squeeze (shape (shape-of spring))
  (elastic (E 200000MPa) (nu 0.3) (yield 1600MPa))          ; spring steel
  (fix base (- z (- 0 (- plate-t 1))))                     ; bottom 1 mm of the base plate
  (fix guide (- (- (+ free-length plate-t) 1) z) (dofs x y))   ; top plate: no sideways motion
  (load push (force 0 0 (* 1N (neg force))) (region (- (- (+ free-length plate-t) 1) z)))
  (res 0.6))

(require deflection-hi (max-displacement squeeze < 2.047))
(require deflection-lo (max-displacement squeeze > 1.675))

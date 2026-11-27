; Benchmark 25: natural frequencies of a cantilever plate and a bracket.
;
; Task (a): an aluminium 6061 plate 100 x 20 x 2, clamped over its first 5 mm
; (x < 0) and free elsewhere. First six natural frequencies.
; Task (b): an L-bracket: a 6 mm aluminium wall flange (30 x 40), held by its
; back face, with an arm 48 x 15 x 4 carrying a 12 mm steel sensor cube on
; its end (x 48..60, centred on the arm). Densities and stiffness per body
; from their materials. Is the first mode above 350 Hz?
;
; Hand calc (a), Euler-Bernoulli, f = beta^2 / (2 pi L^2) sqrt(E I / (rho A)),
; E = 68900 MPa, rho = 2700 kg/m^3, A = 40, I = 13.33 (flat), 1333 (on edge):
;   bending 1 (beta 1.8751): 163.2 Hz  (plate stiffness E/(1-nu^2): 172.9)
;   bending 2 (4.6941): 1022.8 Hz; bending 3 (7.8548): 2863.9 Hz
;   lateral (in-plane) bending 1: 1632.1 Hz
;   torsion n (St-Venant, free warping): (2n-1)/(4L) sqrt(G J / (rho Ip)),
;   J = 0.312 b h^3 = 49.97, Ip = b h (b^2 + h^2) / 12: 1491.6, 4474.9 Hz
; Hand calc (b): the arm as an Euler-Bernoulli beam (L 48) carrying a rigid
; body (13.565 g, centroid 6 mm past the arm's end, Jc = m (12^2+12^2)/12),
; exact frequency equation solved by bisection:
;   vertical (z) bending: 420.6 Hz, second 4243.6 Hz; sideways (y): 1577.1 Hz
;   (Rayleigh with a tip point load at the cube's centroid: 415.9 Hz; the
;   naive "point mass at the beam end" formula gives 496.0, 18% high)
;   total mass 40.78 g.
;
; Pass criteria (ink solve):
;   - (a) bending modes 1-3 within 3% of Euler-Bernoulli; lateral bending and
;     torsion 1 within 8% (shear, clamp restraint of the section);
;   - (b) first mode within 5% of 420.6 Hz; (min-frequency shake > 350)
;     holds; study mass within 0.5% of 40.78 g.

(import "std")

; (a) cantilever plate
(param L 100)
(param b 20)
(param h 2)
(def plate (move (box (+ L 5) b h) (* 0.5 (- L 5)) 0 0))
(body plate (shape plate) (attr material "aluminum-6061") (attr color "#B8C2CC"))

(study ring (shape (shape-of plate))
  (elastic (material aluminum-6061))
  (fix clamp x)
  (modes 6)
  (res 0.5))

; (b) wall bracket with a steel sensor cube
(param arm-len 48)
(param blk 12)
(def flange (move (box 6 30 40) -3 0 0))
(def arm (move (box arm-len 15 4) (* 0.5 arm-len) 0 0))
(def cube (move (box blk blk blk) (+ arm-len (* 0.5 blk)) 0 0))
(body bracket-al (shape (min flange arm)) (attr material "aluminum-6061") (attr color "#B8C2CC"))
(body sensor (shape cube) (attr material "steel") (attr color "#606060"))
(assembly bracket (instance frame bracket-al) (instance sensor sensor))

(study shake (shape (shape-of bracket))
  (elastic (E (* 1MPa (by-body bracket E))) (nu (by-body bracket nu)) (density (* 1g/mm3 (by-body bracket density))))
  (fix wall (face flange (normal -1 0 0)))
  (modes 4)
  (res 0.5))

(require plate-above-150 (min-frequency ring > 150Hz))
; bending 1 within 3% of Euler-Bernoulli's 163.2 Hz (above)
(require bend1-lo (min-frequency ring > 158.3Hz))
(require bend1-hi (min-frequency ring < 168.1Hz))
(require bracket-above-350 (min-frequency shake > 350Hz))

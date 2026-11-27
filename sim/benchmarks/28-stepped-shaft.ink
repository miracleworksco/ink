; Benchmark 28: stepped shaft with shoulder fillets, a keyway and a pulley hub.
;
; Task: a steel drive shaft on two bearings (centres z = 10 and 190, span
; 180): d20 seats (z -20..20 incl. a coupling stub, and 180..200), a d30 body
; (20..180), r3 shoulder fillets (std's fillet-union), and an end-milled
; 8 x 4 keyway at the pulley (z 85..115, on +y). The pulley puts a belt pull
; P = 3000 N (-x) and a torque T = 100 N m on the shaft at z = 100; the torque
; is taken out at the coupling end (z = -20). A separate aluminium hub (bore
; 15.01) and a 7.98 x 7 key form the assembly.
;
; Hand calculation (steel E 200 GPa, nu 0.3):
;   - d30: I = pi 30^4/64 = 39761 mm^4, J = 79522; d20: I = 7854, J = 15708.
;   - simply supported, R = P/2 = 1500 N; M(z) = 1500 (z - 10) up to z = 100.
;   - z = 60, r = 14: sigma = M r / I = 75000 x 14 / 39761 = 26.41 MPa;
;     tau = T r / J = 17.61 MPa; von Mises sqrt(s^2 + 3 t^2) = 40.34 MPa.
;   - coupling stub z = 0, r = 9 (no bending): tau = 57.30, von Mises 99.2.
;   - left shoulder z = 20 (d/D = 20/30, r/d = 0.15): nominal sigma =
;     32 M / pi d^3 = 19.10 (M = 15000), tau = 16 T / pi d^3 = 63.66.
;     Peterson (Shigley's power-law fits of the charts, D/d = 1.5):
;     Kt bending = 0.938 (r/d)^-0.258 = 1.53, Kts torsion ~ 0.855 (r/d)^-0.235
;     = 1.34 (interpolated between the D/d 1.33 and 2 fits);
;     peak von Mises ~ sqrt((1.53 x 19.10)^2 + 3 (1.34 x 63.66)^2) = 150 MPa.
;   - right shoulder z = 180: pure bending, 1.53 x 19.10 = 29.2 MPa.
;   - deflection at z = 100 (stepped EI, unit-load integral, python):
;     bending 0.0461 + shear (kappa 0.886) 0.0032 = 0.0493 mm.
;   - fits: hub bore 15.01 on a 15.000 shaft (0.010 radial clearance); key
;     7.98 in an 8.00 keyway (0.01 per side) and an 8.04 hub slot (0.03), key
;     floor gap 0.005; key reaches 3 mm into the hub slot.
;
; Pass criteria:
;   - von Mises at (14, 0, 60) within 5% of 40.34, bending and torsion parts
;     separately within 5%; stub torsion within 5% of 57.3;
;   - peak von Mises at the left shoulder within 25% of Peterson (150);
;     right shoulder (bending only) within 25% of 29.2;
;   - deflection at z = 100 within 10% of 0.0493;
;   - hub and key clear the shaft (proven), the key engages shaft and hub
;     (proven), a hub bored to 14.99 is disproven; interfere: no overlaps.
; expect-status: control-tight-hub-disproven fails

(import "std")

(param fillet-r 3)           ; shoulder fillet radius
(param hub-bore-r 15.01)     ; pulley hub bore (shaft 15.000: 0.01 radial clearance)

; ---------------------------------------------------------------- shaft (axis z)
; Coupling stub + left bearing seat d20 over z -20..20, body d30 over 20..180,
; right bearing seat d20 over 180..200. Bearings centred at z = 10 and 190.
; The seats overlap the body by 10 mm so the shoulder is one concave edge.
(def seat-l (move (cylinder-z 10 50) 0 0 5))          ; z -20 .. 30
(def body-mid (move (cylinder-z 15 160) 0 0 100))     ; z 20 .. 180
(def seat-r (move (cylinder-z 10 30) 0 0 185))        ; z 170 .. 200
(def filleted (fillet-union (fillet-union seat-l body-mid fillet-r 90) seat-r fillet-r 90))
; End-milled keyway on +y: 8 wide, 4 deep (floor y = 11), z 85..115.
(def keyway (max (- (len2 x (- z (clamp z 89 111))) 4) (- 11 y)))
(def shaft-field (subtract filleted keyway))
(body shaft (shape shaft-field) (attr material "steel") (attr density 7.85g/cm3) (attr color "#8A8F96"))

; ---------------------------------------------------------------- pulley hub + key
; Hub OD 70, z 88..112, bore hub-bore-r, through key slot 8.04 wide up to y 18.3.
(def hub-slot (max (- (abs x) 4.02) (neg y) (- y 18.3)))
(def hub-bore (- (len2 x y) hub-bore-r))
(def hub-field (subtract (move (cylinder-z 35 24) 0 0 100) (min hub-bore hub-slot)))
(body hub (shape hub-field) (attr material "aluminum-6061") (attr color "#C8A06A"))
; Key 7.98 wide, 7 tall (y 11.005 .. 18.005), round ends, z 86.01 .. 113.99.
(def key-field (max (- (len2 x (- z (clamp z 90 110))) 3.99) (- 11.005 y) (- y 18.005)))
(body key (shape key-field) (attr material "steel") (attr color "#404040"))

(assembly root
  (instance shaft shaft)
  (instance hub hub)
  (instance key key))

; ---------------------------------------------------------------- fit claims
; running clearance hub/shaft (0.010 radial by design). Written as disjoint
; from the grown shaft: the (clearance a b d) form of the same claims stays
; unknown down to res 0.05 (see Status).
(require hub-clears-shaft (disjoint (- (shape-of root/shaft) 0.005) (shape-of root/hub)))
; ... and not sloppy: hub material within r 15.02 all round, away from the slot
(require hub-bore-tight (subset (max (extrude-z (ring-2d 15.02 15.03) 88.5 111.5) (- y -5)) (shape-of root/hub)))
; key: side clearance 0.01 in the shaft keyway, 0.03 in the hub slot
(require key-clears-shaft (disjoint (- (shape-of root/shaft) 0.004) (shape-of root/key)))
(require key-clears-hub (disjoint (- (shape-of root/hub) 0.02) (shape-of root/key)))
; the claim form the language offers for this, kept for the record
(require key-clears-hub-clearance-form (clearance (shape-of root/hub) (shape-of root/key) 0.02))
; key engages both: it fills y 12..14 (inside the shaft keyway) and y 16..17.5
; (inside the hub slot, above the bore) over the hub length
(def key-low (move (box 6 2 20) 0 13 100))
(def key-high (move (box 6 1.5 20) 0 16.75 100))
(require key-in-keyway (subset key-low (shape-of root/key)))
(require key-in-hub-slot (subset key-high (shape-of root/key)))
(require key-high-in-slot (disjoint key-high (shape-of root/hub)))
(require key-high-outside-shaft (disjoint key-high (shape-of root/shaft)))
; control: a hub bored to 14.99 (interference) must be disproven
(def tight-hub (subtract (move (cylinder-z 35 24) 0 0 100) (min (- (len2 x y) 14.99) hub-slot)))
(require control-tight-hub-disproven (disjoint shaft-field tight-hub) (accuracy 0.002))

; ---------------------------------------------------------------- load case
; Belt pull P = 3000 N (-x) at the hub, torque T = 100 N m delivered at the hub
; and taken out at the coupling end (z = -20). Each bearing holds the outer
; 2 mm of its seat over its 2 mm width as a ball joint, every turn about its
; centre free (a self-aligning bearing: no torque, no bending moment); the
; right one also slides axially. The coupling face holds y (reacts torque).
; Torque is a couple of y forces on the material at |x| > 12 (a circular
; segment whose centroid is at 13.2109 mm): Ft = T / (2 x 13.2109).
(def P 3000)
(def T 100000)
(def seg-c 12)
(def xbar 13.210919583)
(def Ft (/ T (* 2 xbar)))
(fn seat (zc) (max (- (abs (- z zc)) 1) (- 8 (len2 x y))))
(study drive (shape shaft-field) (elastic (material steel))
  (fix brg-l (seat 10) (pivot (at 0 0 10)))
  (fix brg-r (seat 190) (pivot (at 0 0 190)) (slide 0 0 1))
  (fix coupling (- z -19) (dofs y) (support))
  (load belt (force (* 1N (neg P)) 0 0) (region (- (abs (- z 100)) 5)))
  (load torque-a (force 0 (* 1N Ft) 0) (region (max (- (abs (- z 100)) 5) (- seg-c x))))
  (load torque-b (force 0 (* 1N (neg Ft)) 0) (region (max (- (abs (- z 100)) 5) (+ x seg-c))))
  (res 1))


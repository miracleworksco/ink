; Benchmark 07: bolted flange joint, ISO threads, placed by frames.
;
; Task: two plain pipe flanges (OD 100, bore 40, 12.2 thick) clamped by four
; M10x1.5 hex bolts with hex nuts on a 75 mm pitch circle. Clearance holes
; are ISO 273 medium (11 mm). A bolt and its nut form a set. The bolt's
; under-head frame is the set's `head` frame. The nut sits on a frame one
; grip above it, turned so its thread engages: a nut moved dz along the bolt
; turns 360 dz / p degrees. Each set is seated by `head` on the lower
; flange's hole frame, turned about the flange axis per hole.
;
; Reference values:
;   - ISO 273 medium clearance: hole 11.0 vs M10 major 10.0 -> 0.5 mm radial.
;   - M10x1.5 minor dia (ISO 724, external) d3 = d - 1.2269 p = 8.160;
;     thread-iso uses d - 1.0825 p = 8.376 (the basic profile), so the
;     core radius should be 4.188.
;   - Grip = 2 x 12.2 = 24.4 mm (16.27 pitches: the nut phase is non-trivial);
;     bolt 40 long, nut 8 high -> 7.6 mm protrusion.
;   - Nut tapped with 0.1 mm clearance on the thread's distance field.
;
; Pass criteria (ink assert --res 0.025, ink interfere):
;   - bolt-in-hole clearance >= 0.45 mm and bolt thread vs nut thread
;     clearance >= 0.05 mm proven;
;   - interfere at --res 0.2: no overlapping pairs, none left `unknown`.
;     The flanges' faces touch over 6600 mm2; the search's budget covers
;     that at cells of 0.2 to 0.4 mm, not at the model's finest;
;   - bolt protrudes past the nut; hole radius 5.5 +- 0.01 (probes: open
;     to r 5.49, metal from r 5.51, proven).

(import "std")

(param d 10)          ; bolt major diameter
(param p 1.5)         ; pitch
(param t-fl 12.2)      ; flange thickness
(param pcd 75)        ; bolt circle
(param hole-d 11)     ; clearance hole
(param bolt-l 40)
(param thread-l 26)   ; threaded length from the tip
(param nut-h 8)
(param af 16)         ; hex across flats (head and nut)
(param tap-clear 0.1) ; nut thread clearance

; ---------------------------------------------------------------- flanges
(def flange-disc (max (cylinder-z 50 t-fl) (- 20 (len2 x y))))
(def holes (pattern-polar (move (cylinder-z (* 0.5 hole-d) 100) (* 0.5 pcd) 0 0) 4))
(def flange-field (move (subtract flange-disc holes) 0 0 (* 0.5 t-fl)))   ; z in [0, t-fl]
(body flange (shape flange-field) (attr material "steel") (attr color "#8FA3B8")
  (frame hole-under (at (rot-x 180) (move (* 0.5 pcd) 0 0)))   ; underside at the +x hole, z down
  (frame hole-top (at (move (* 0.5 pcd) 0 t-fl))))             ; top face at the +x hole, z up

; ---------------------------------------------------------------- bolt and nut
(def head (extrude-z (hexagon-2d af) -6.4 0))
(def shank (move (cylinder-z (* 0.5 d) (- bolt-l thread-l)) 0 0 (* 0.5 (- bolt-l thread-l))))
(def threads (max (thread-iso d p) (slab-z (- bolt-l thread-l) bolt-l)))
(body bolt (shape (min head shank threads)) (attr material "steel") (attr color "#606060")
  (frame seat (at)))                                  ; under the head, z along the shank

(def nut-block (extrude-z (hexagon-2d af) 0 nut-h))
(body nut (shape (subtract nut-block (- (thread-iso d p) tap-clear))) (attr material "steel") (attr color "#909090")
  (frame seat (at (rot-x 180))))                      ; bottom face, z down

; ---------------------------------------------------------------- assembly
(def grip (* 2 t-fl))
(assembly bolt-set
  (instance bolt bolt)
  (instance nut nut (at (move 0 0 grip)) seat (turn z (* 360 (/ grip p))))
  (frame head (frame-of bolt seat)))

(assembly root
  (instance lower flange)
  (instance upper flange (at (move 0 0 t-fl)))
  (instance set0 bolt-set (frame-of lower hole-under) head)
  (instance set1 bolt-set (at (frame-of lower hole-under) (rot-z 90)) head)
  (instance set2 bolt-set (at (frame-of lower hole-under) (rot-z 180)) head)
  (instance set3 bolt-set (at (frame-of lower hole-under) (rot-z 270)) head))

; ---------------------------------------------------------------- claims
(def bolt-below-head (max (shape-of root/set0/bolt) (- 0.1 z)))   ; the head seats: leave it out
(assert bolt-clears-hole (clearance (shape-of root/lower) bolt-below-head 0.45))
(assert bolt-clears-hole-upper (clearance (shape-of root/upper) bolt-below-head 0.45))
(assert thread-engages-with-clearance (clearance (shape-of root/set0/nut) (shape-of root/set0/bolt) 0.05))
; The nut really engages: its material reaches inside the bolt's major radius.
; There is no must-overlap claim (no negation), so this is a disjoint claim
; that must come back DISPROVEN:
; expect-status: control-nut-reaches-thread-disproven fails
(assert control-nut-reaches-thread-disproven
  (disjoint (shape-of root/set0/nut) (move (cylinder-z (- (* 0.5 d) 0.3) 100) (* 0.5 pcd) 0 0)))
; bolt tip protrudes past the nut: no nut material above the plane of the tip minus 1 pitch
(assert tip-protrudes (disjoint (shape-of root/set0/nut) (half-space 0 0 -1 (- 0 (+ bolt-l (neg p))))))
; the +x hole of the lower flange: open to r 5.49, metal from r 5.51
(assert hole-open (disjoint (shape-of root/lower) (move (cylinder-z 5.49 10) 37.5 0 6.1)))
(assert hole-metal (subset (move (subtract (cylinder-z 6 10) (cylinder-z 5.51 11)) 37.5 0 6.1) (shape-of root/lower)))

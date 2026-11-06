; Benchmark 02: spur gear pair with involute teeth.
;
; Task: module 2, 20 and 30 teeth, 20 deg pressure angle, 10 mm thick, 0.1 mm
; backlash at the pitch circle, centre distance m (z1 + z2) / 2 = 50 mm. Show
; that the pair meshes without interference over one tooth pitch. The flanks
; stay about 0.025 mm apart, so the proof needs a fine res.
;
; Pass criteria (ink assert --res 0.02):
;   - mesh-sweep: the pair is disjoint at every 0.5 deg over a full pinion
;     pitch (18 deg) (proven);
;   - the wrong-phase control (the gear not turned by half a pitch) disproven.
; expect-warning: the file reads t but states no t-span

(import "std")

(param m 2)                 ; module
(param z1 20)               ; teeth, pinion
(param z2 30)               ; teeth, gear
(param pa 20)               ; pressure angle (deg)
(param width 10)            ; face width
(param backlash 0.1)        ; total circumferential backlash at the pitch circle (mm)
(param bore 4)              ; shaft bore radius
(param phi 0)               ; pinion rotation (deg); the gear follows at -phi z1/z2

(def centre (* 0.5 m (+ z1 z2)))

; std's involute-gear-2d: tight under interval arithmetic, so the mesh
; clearance proves.
(fn gear-2d (z) (max (involute-gear-2d z m pa (* 0.5 backlash)) (- bore (len2 x y))))

(def pinion-field (extrude-z (gear-2d z1) 0 width))
(def gear-field (extrude-z (gear-2d z2) 0 width))

(body pinion (shape pinion-field) (attr color "#D08A2E") (attr density 7.85g/cm3))
(body gear (shape gear-field) (attr color "#3A7BC8") (attr density 7.85g/cm3))

; Pinion at the origin, gear at (centre, 0, 0). With z2 = 30, a gear tooth sits
; at 180 deg; turning it by half a pitch (180/z2) puts a gap facing the pinion's
; tooth at 0 deg. Rotation is driven by the param phi plus time t (1 rev/s), so
; `--param phi=...` or `--t ...` both sweep the mesh.
(assembly root
  (instance pinion pinion (at (rot-z (+ phi (* 360 t)))))
  (instance gear gear (at (rot-z (+ (/ 180 z2) (* (/ z1 z2) (neg (+ phi (* 360 t)))))) (move centre 0 0))))

; Meshing at several angles inside one tooth pitch (18 deg for the pinion),
; written as claims so a single `ink assert` checks them all.
(fn pinion-at (a) (rotate-z pinion-field a))
(fn gear-at (a) (move (rotate-z gear-field (- (/ 180 z2) (* a (/ z1 z2)))) centre 0 0))
(assert mesh-0 (disjoint (pinion-at 0) (gear-at 0)))
(assert mesh-3 (disjoint (pinion-at 3) (gear-at 3)))
(assert mesh-6 (disjoint (pinion-at 6) (gear-at 6)))
(assert mesh-9 (disjoint (pinion-at 9) (gear-at 9)))
(assert mesh-12 (disjoint (pinion-at 12) (gear-at 12)))
(assert mesh-15 (disjoint (pinion-at 15) (gear-at 15)))
; The same check over a whole tooth pitch in one claim: the union over 36
; angles (0, 0.5, ... 17.5 deg) of the pair's overlap is empty, one search.
; (for a (range 0 18 0.5) (assert mesh-$a ...)) would write it as 36 claims,
; each answered on its own.
(assert mesh-sweep
  (empty (reduce min (map (lambda (a) (max (pinion-at a) (gear-at a))) (range 0 18 0.5)))))
; Negative control, EXPECTED TO BE DISPROVEN: with the wrong phase (tooth facing
; tooth) the gears collide, so their disjoint claim fails with a witness.
; (overlaps a b) would claim the collision instead.
; expect-status: control-wrong-phase-disproven fails
(assert control-wrong-phase-disproven
  (disjoint (pinion-at 0) (move (rotate-z gear-field 0) centre 0 0)))
; Tooth thickness at the pitch circle: the arc from the tooth centre line to the
; flank is (pi m / 4) - backlash/4 = 1.5458 mm, i.e. the flank crosses r = 20 at
; 4.4285 deg (checked with ink eval in the notes). Tip and root radii:
; The pinion's faces lie on its envelope's: no search separates the faces at
; the model's floor within its budget, so it is decided at cells of 0.1 (to
; 0.2 mm).
(assert pinion-envelope (within pinion-field (-22 -22 0) (22 22 10)) (res 0.1))

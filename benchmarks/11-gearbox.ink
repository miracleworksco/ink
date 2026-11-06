; Benchmark 11: two-stage spur gearbox with proven meshing.
;
; Task: input pinion (18 T, m 1.5) drives a cluster gear (36 T) whose second
; pinion (18 T, m 2) drives the output gear (54 T): ratio 2 x 3 = 6. The
; output shaft sits at 90 deg from the first stage's line of centres, so the
; box is L-shaped. Stage 1 at z 0..10, stage 2 at z 12..24. Every gear is the
; analytic std involute-gear-2d, 20 deg pressure angle, 0.1 mm backlash per
; pair. Time t turns the input by one tooth pitch (20 deg) per unit time.
;
; Reference values:
;   - centre distances a1 = 1.5 (18 + 36) / 2 = 40.5, a2 = 2 (18 + 54) / 2 = 72;
;   - phases: a gear at bearing beta from its pinion turns as
;     phi_g = beta + 180 + 180/z_g - (z_p/z_g) (phi_p - beta);
;   - ratio: the output turns 1/6 of the input, in the same sense (two
;     external meshes): +3.333 deg per 20 deg at the input.
;
; Pass criteria (ink assert --res 0.02, ink interfere over the file's t-span):
;   - stage 1 and stage 2 mesh proven interference-free at 21 angles over a
;     full input tooth pitch (folded claims, like benchmark 02);
;   - wrong-phase controls disproven;
;   - interfere over t in [0, 1] (all 4 bodies moving): no overlaps and no
;     pair left unknown.

(import "std")
(t-span 0 1)

(param m1 1.5) (param m2 2)
(param zp1 18) (param zg1 36) (param zp2 18) (param zg2 54)
(param backlash 0.1)
(def a1 (* 0.5 m1 (+ zp1 zg1)))
(def a2 (* 0.5 m2 (+ zp2 zg2)))
(def pa 20)

(fn gear (z mo z0 z1 bore-r)
  (extrude-z (max (involute-gear-2d z mo pa (* 0.5 backlash)) (- bore-r (len2 x y))) z0 z1))
(fn shaft (r z0 z1) (extrude-z (circle-2d r) z0 z1))

(def pinion1-field (min (gear zp1 m1 0 10 0) (shaft 5 -12 0)))
(def cluster-field (min (gear zg1 m1 0 10 0) (gear zp2 m2 12 24 0) (shaft 6 -12 30)))
(def output-field (min (gear zg2 m2 12 24 0) (shaft 8 -12 30)))
(body pinion1 (shape pinion1-field) (attr material "steel") (attr color "#D08A2E"))
(body cluster (shape cluster-field) (attr material "steel") (attr color "#3A7BC8"))
(body output (shape output-field) (attr material "steel") (attr color "#6A9F58"))

; angles (deg) as functions of the input angle
(def phi-in (* (/ 360 zp1) t))
(fn follow (beta zp zg phi-p) (+ beta 180 (/ 180 zg) (* (neg (/ zp zg)) (- phi-p beta))))
(def phi-c (follow 0 zp1 zg1 phi-in))
(def phi-out (follow 90 zp2 zg2 phi-c))           ; pinion 2 shares the cluster's angle

(assembly root
  (instance input pinion1 (at (rot-z phi-in)))
  (instance cluster cluster (at (rot-z phi-c) (move a1 0 0)))
  (instance output output (at (rot-z phi-out) (move a1 a2 0))))

; ---------------------------------------------------------------- meshing at sampled angles
(def stage1-p (extrude-z (involute-gear-2d zp1 m1 pa (* 0.5 backlash)) 0 10))
(def stage1-g (extrude-z (involute-gear-2d zg1 m1 pa (* 0.5 backlash)) 0 10))
(def stage2-p (extrude-z (involute-gear-2d zp2 m2 pa (* 0.5 backlash)) 12 24))
(def stage2-g (extrude-z (involute-gear-2d zg2 m2 pa (* 0.5 backlash)) 12 24))
(fn s1-at (p) (max (rotate-z stage1-p p) (move (rotate-z stage1-g (follow 0 zp1 zg1 p)) a1 0 0)))
(fn s2-at (p) (let ((c (follow 0 zp1 zg1 p)))
  (max (move (rotate-z stage2-p c) a1 0 0) (move (rotate-z stage2-g (follow 90 zp2 zg2 c)) a1 a2 0))))
(assert stage1-mesh (empty (reduce min (map s1-at (linspace 0 20 21)))))
(assert stage2-mesh (empty (reduce min (map s2-at (linspace 0 20 21)))))
; negative controls (must be DISPROVEN): gear turned half a pitch
; expect-status: control-stage1-disproven fails
(assert control-stage1-disproven (empty (max (rotate-z stage1-p 0) (move (rotate-z stage1-g 0) a1 0 0))))
; expect-status: control-stage2-disproven fails
(assert control-stage2-disproven
  (empty (max (move stage2-p a1 0 0) (move (rotate-z stage2-g (+ 180 (follow 90 zp2 zg2 5))) a1 a2 0))))

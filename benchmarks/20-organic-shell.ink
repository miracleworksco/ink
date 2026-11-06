; Benchmark 20: a sculpted, hollow part with smooth blends (a mouse shell).
;
; Task: an organic palm body blended from ellipsoids (smooth-union), a
; thumb scoop cut with smooth-subtract, a flat base with a rounded edge,
; then hollowed to a 2 mm wall. Ellipsoids are non-uniformly scaled
; spheres: not distance fields, so a naive (shell f 2) wall varies in
; thickness; (redistance f) first should make it exact.
;
; Reference: the wall is 2.0 mm everywhere (by construction intent).
;
; Pass criteria (ink mesh, ink components, ink ray):
;   - solid body: watertight, 1 solid;
;   - hollow (redistanced) shell: 1 solid + 1 cavity, watertight;
;   - wall thickness along the probe rays within 0.1 mm of 2.0 (redistanced);
;   - the naive shell's wall is off by > 0.25 somewhere (the control).
; (the benchmark's test probes 2 ray directions: b20 in
; crates/ink-cli/tests/cli/benchmarks.rs.)

; The hump and the thumb scoop lie wholly inside their smooth blends, so the
; surface there is the blend's, not theirs. `check` proves it for the hump;
; the thumb scoop is unknown, not reported: interval bounds cannot rule out
; its plain-subtract side of the blend over any cell of the scoop's surface.
; expect-warning: body 'mouse': 'hump' is in its field but forms none of its surface
; expect-warning: body 'naive-shell': 'hump' is in its field but forms none of its surface

; The shells stand beside the solid body they are cut from.

(import "std")

(fn ellipsoid (a b c) (* (min a b c) (warp (sphere 1) (/ x a) (/ y b) (/ z c))))
(def palm (move (ellipsoid 55 32 24) 0 0 0))
(def hump (move (ellipsoid 28 24 20) -12 0 10))
(def thumb-scoop (move (ellipsoid 30 10 14) 5 34 8))
(def top (smooth-subtract (smooth-union palm hump 12) thumb-scoop 6))
(def base-plane (- 2 z))                                    ; keep z >= 2
(def sculpt (smooth-intersect top base-plane 3))
(body mouse (shape sculpt) (attr color "#5A6B7C"))

(def exact (redistance sculpt (res 0.4)))
(body mouse-shell (shape (move (shell exact 2) 0 120 0)) (attr color "#7C8B9A"))
(body naive-shell (shape (move (shell sculpt 2) 0 240 0)) (attr color "#9A7C7C"))

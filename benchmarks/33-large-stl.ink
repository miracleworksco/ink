; Benchmark 33: a large imported mesh (58 080 triangles) checked against its
; exact geometry and two mating parts.
;
; Task: a supplier sends a binary STL of a knurled, finned knob. Import it
; with import-stl and check its volume. Check clearances to a shaft and a
; shroud and a press fit with a ring. Probe how much of the knurl
; survives sampling, and mesh it back out.
; inputs/33-knob.stl (2.9 MB) is written by
; python3 tools/33-gen-knob.py inputs/33-knob.stl 240. The benchmark's test
; (b33 in crates/ink-cli/tests/cli/benchmarks.rs) checks it writes the same
; file. Its default n 2160 (523k triangles, 26 MB) is the same part on a
; finer grid.
;
; The part (mm, axis z, z 0 .. 60), every section a polygon on one angular
; grid of 240 samples:
;   - bore: r 4 through;
;   - grip z 0 .. 20: straight knurl, 120 teeth, tips r 15.0, roots r 14.4
;     (a 240-gon; pitch 0.785 mm at the tips, 0.6 deep);
;   - hub z 20 .. 60: r 8;
;   - 19 fins r 22, 1 thick at a 2 pitch (z 22..23, ..., 58..59): 1 mm gaps.
; Hand calculation (33-gen-knob.py prints it):
;   knurl section = 240 triangles 0.5 15 14.4 sin(1.5 deg) = 678.5065;
;   V = 20 (678.5065) + pi 8^2 (2 + 19) + 19 pi 22^2 - 60 pi 4^2
;     = 13 570.13 + 4 222.30 + 28 890.09 - 3 015.93 = 43 666.59 for true
;   circles; the mesh's 240-gons: 43 663.15 (-7.9e-5).
;   Mating parts: shaft r 3.9 in the r 4 bore (0.1 radial); shroud bore r 22.4
;   over the fins (0.4 radial); ring bore r 14.8 on the knurl (0.2 bite).
;   Grid sampling error: the mesh's signed distance d is 1-Lipschitz and the
;   trilinear sample I matches it at the grid points, so
;   |I - d| <= sqrt(sum w_c |p - c|^2) <= (sqrt 3 / 2) h = 0.217 at h 0.25.
;
;   Claims on the import are measured on its samples (docs/IR.md, "Sampled
;   fields"): a margin under that bound is within the sampling error. The
;   knurl's 0.05 margins are checked again on a region sampled at 0.025
;   (bound 0.022), and every claim on the exact part, which is guaranteed.
;
; Pass criteria (import res 0.25):
;   - ink check: no warnings but the exact part's, the shaft's and the
;     shroud's, declared below. The import reports a measured deviation
;     below the sqrt 3 / 2 h above (0.2165); its samples, in floats, are
;     exact to far less;
;   - ink mass --res 0.3: estimate within 0.2% of the mesh's 43 663.15:
;     the samples' own volume, 43 619.5 at res 0.15 and 0.1 (proven within
;     43 538 .. 43 706), is 0.1% under the mesh's, the trilinear samples
;     rounding the knurl's tips and the fins' edges;
;   - estimates, holding (holds, accuracy null): the clearances, the press
;     fit, the knurl probes on the import and fine-core-full-14.35;
;   - estimates, failing (fails, accuracy null): the three controls on the
;     import, and fine-control-core-full-14.6 (0.2 short at the roots);
;   - the exact-* claims hold to their accuracy, exact-control-core-full-14.5
;     fails to it;
;   - the knurl probes that test a tooth (ring-grips-knurl, tips-reach-14.95,
;     control-core-full-14.5 and their exact-* twins) say (accuracy 0.01):
;     a sharp tooth overlaps the ring or r 14.95, and a root's notch cuts
;     into r 14.5, less deep than the default 0.05 (on the exact part 0.018
;     and 0.035 mm), so at it they read as unknown or holding;
;   - ink mesh: watertight, 1 shell, volume within 0.1%;
;   - rays (ink ray --rays, from r 20 at z 10 toward the axis): tips
;     >= 14.75, roots <= 14.65 (within the sampling bound), fins, gaps and
;     bore within 0.005.

; The exact part stands beside the import, for comparison.

(import "std")

(param import-res 0.25)
(def knob-field (import-stl "inputs/33-knob.stl" (res import-res)))
(body knob (shape knob-field) (attr color "#B0B8C0"))

; ---------------------------------------------------------------- the exact part, for comparison (at x = 100)
(def knurl-pts (map (lambda (k) (let ((r (if (== (mod k 2) 0) 15 14.4)) (an (/ (* 2 pi k) 240)))
                                  (vec (* r (cos an)) (* r (sin an)))))
                    (range 240)))
(def grip (extrude-z (polygon knurl-pts) 0 20))
(def hub (extrude-z (circle-2d 8) 20 60))
(def fins (reduce min (map (lambda (k) (extrude-z (circle-2d 22) (+ 22 (* 2 k)) (+ 23 (* 2 k)))) (range 19))))
(def bore (cylinder-z 4 200))
(body exact (shape (move (subtract (min grip hub fins) bore) 100 0 0)) (attr color "#8FA3B8"))

; ---------------------------------------------------------------- mating parts
; The shaft and shroud stand in their clearances, touching nothing.
(body shaft (shape (move (cylinder-z 3.9 70) 0 0 30)) (attr color "#D79B00"))      ; 0.1 radial clearance in the bore
(def shroud-field (extrude-z (ring-2d 22.4 25) 21 61))     ; 0.4 radial clearance to the fins
(body shroud (shape shroud-field) (attr color "#6C8EBF"))
(def ring-field (extrude-z (ring-2d 14.8 18) 2 18))       ; inner r 14.8: tips 15.0 bite 0.2
(body ring (shape ring-field) (attr color "#82B366"))

; probes of the knurl: annuli r 6 .. ro and a band r ri .. 16, z 2 .. 18
(fn core (ro) (extrude-z (ring-2d 6 ro) 2 18))
(fn band (ri) (extrude-z (ring-2d ri 16) 2 18))
; the 1 mm gap between fins 1 and 2 (z 23 .. 24), r 10 .. 20
(def gap-probe (extrude-z (ring-2d 10 20) 23.2 23.8))

; ---------------------------------------------------------------- claims
; mating parts
(assert shaft-disjoint (disjoint (shape-of knob) (shape-of shaft)))
; expect-status: control-shaft-0.15 fails estimate
(assert control-shaft-0.15 (clearance (shape-of knob) (shape-of shaft) 0.15))
(assert shaft-clear-0.05 (clearance (shape-of knob) (shape-of shaft) 0.05))
(assert fins-clear-shroud (clearance (shape-of knob) (shape-of shroud) 0.1))
(assert shroud-clear-0.2 (clearance (shape-of knob) (shape-of shroud) 0.2))
; expect-status: control-fins-0.45 fails estimate
(assert control-fins-0.45 (clearance (shape-of knob) (shape-of shroud) 0.45))
(assert ring-grips-knurl (overlaps (shape-of knob) (shape-of ring)) (accuracy 0.01))
; knurl detail (tips r 15.0, roots r 14.4, 0.785 mm pitch): true of the mesh
(assert core-full-14.35 (subset (core 14.35) (shape-of knob)))
(assert tips-reach-14.95 (overlaps (shape-of knob) (band 14.95)) (accuracy 0.01))
(assert fin-gap-open (disjoint (shape-of knob) gap-probe))
; false of the mesh: the roots cut in to r 14.4
; expect-status: control-core-full-14.5 fails estimate
(assert control-core-full-14.5 (subset (core 14.5) (shape-of knob)) (accuracy 0.01))
; A region sampled finer: the knurl's teeth near +x, at res
; 0.025 (error bound 0.022), where its 0.05 margins clear the error.
(def knurl-fine (import-stl "inputs/33-knob.stl" (res 0.025) (region (13 -3 2) (16 3 18))))
(def teeth (move (box 2 4 14) 14.5 0 10))     ; x 13.5 .. 15.5, y -2 .. 2, z 3 .. 17: inside the region
(assert fine-core-full-14.35 (subset (max (core 14.35) teeth) knurl-fine))
; expect-status: fine-control-core-full-14.6 fails estimate
(assert fine-control-core-full-14.6 (subset (max (core 14.6) teeth) knurl-fine))   ; 0.2 short at the roots
; the same claims on the exact part
(assert exact-shaft-clear (clearance (shape-of exact) (move (shape-of shaft) 100 0 0) 0.05))
(assert exact-core-full-14.35 (subset (move (core 14.35) 100 0 0) (shape-of exact)))
(assert exact-tips-reach-14.95 (overlaps (shape-of exact) (move (band 14.95) 100 0 0)) (accuracy 0.01))
(assert exact-fin-gap-open (disjoint (shape-of exact) (move gap-probe 100 0 0)))
; expect-status: exact-control-core-full-14.5 fails
(assert exact-control-core-full-14.5 (subset (move (core 14.5) 100 0 0) (shape-of exact)) (accuracy 0.01))

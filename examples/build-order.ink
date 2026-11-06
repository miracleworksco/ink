; Proves a build order step by step (docs/IR.md, Assembly order): each part
; comes in from outside the parts placed before it, clear of them. A block
; is lowered onto a plate. Two M3 screws go up through the plate into the
; block's 5 mm blind taps, touching the block only in their thread's band,
; where thread meets tap. The 8 mm screw ends 4 mm in. The 12 mm one would
; end 8 mm in, past its tap into the block:
; expect-status: long-screw-into-block fails
(import "std")
; An ISO 4762 M3 socket head cap screw l long: its head (5.5 across, 3 tall,
; a 2.5 hex socket 1.3 deep) above its bearing face at z = 0 (seat, z out of
; the head, down the screw), its thread (M3 x 0.5, std's thread-iso, in
; phase with a tap's at along-thread) to its tip at z = -l; drive at the
; socket's bottom. Its thread is an attr, read where it turns into a part.
(fn m3-screw (l)
  (min (subtract (extrude-z (circle-2d 2.75) 0 3) (extrude-z (hexagon-2d 2.5) 1.7 4))
       (max (thread-iso 3 0.5) (slab-z (neg l) 0.1))))
(body m3x8 (shape (m3-screw 8)) (attr thread (vec 3 0.5)) (attr color "#2A2A2E")
  (frame seat (at (rot-x 180))) (frame drive (at (move 0 0 1.7))))
(body m3x12 (shape (m3-screw 12)) (attr thread (vec 3 0.5)) (attr color "#2A2A2E")
  (frame seat (at (rot-x 180))) (frame drive (at (move 0 0 1.7))))
(def tap (max (- (thread-iso 3 0.5) 0.05) (slab-z -1 5)))
(body plate (shape (subtract (move (box 30 30 4) 0 0 -2) (min (move (cylinder-z 1.7 20) -3 0 0) (move (cylinder-z 1.7 20) 3 0 0))))
  (frame top (at)))
(body block (shape (subtract (move (box 12 12 8) 0 0 4) (min (move tap -3 0 0) (move tap 3 0 0))))
  (frame bottom (at (rot-x 180))) (frame tap-0 (at (rot-x 180) (move -3 0 0))) (frame tap-1 (at (rot-x 180) (move 3 0 0))))
(def out (ranged 0 1))                       ; along each path: 0 at rest, 1 at its start
(assembly built (instance plate plate) (instance block block (frame-of plate top) bottom)
  (instance short-screw m3x8 (at (along-thread 4 .5) (frame-of block tap-0)) seat)    ; up through the plate
  (instance long-screw m3x12 (at (along-thread 4 .5) (frame-of block tap-1)) seat))
(def order (vec (build-laid "plate")
                (build-step "block" (vec 0 0 1))
                (build-step "short-screw" (vec 0 0 -1) "block")      ; from below, into the block's tap
                (build-step "long-screw" (vec 0 0 -1) "block")))
(for k (range 1 (len order))
  (for p (vec (nth (nth order k) 0))
    (def path-$k (let ((part (shape-of built/$p)) (s (nth order k)))
      (build-along part s (build-length s part (bounds (shape-of built (except (for q (build-later order k) built/$q))))) out)))
    (assert $p (disjoint path-$k (shape-of built (except (for q (build-left-out order k) built/$q)))))
    (for m (build-into order k)
      (assert $p-into-$m (disjoint path-$k (max (shape-of built/$m)
        (neg (place (apply thread-band (attr-of built/$p thread)) (frame-of built/$p drive)))))))))

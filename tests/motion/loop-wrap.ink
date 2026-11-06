; Links carried round a closed path over t: each link's station wraps with
; `mod` at the loop's start, where the frame closes continuously. Bounded
; as a choice of the wrap's two sides (eval::forked), a step crossing the
; start is one link either side of it, not the link anywhere on the loop.
(import "std")
(def P (path-loop (vec (path-line (vec 0 -10 0) (vec 50 -10 0))
                       (path-arc (vec 50 0 0) (vec 0 0 1) 10 (vec 0 -1 0) pi)
                       (path-line (vec 50 10 0) (vec 0 10 0))
                       (path-arc (vec 0 0 0) (vec 0 0 1) 10 (vec 0 1 0) pi)) (vec 0 0 1)))
(body link (shape (box 4 2 2)))
(t-span 0 2)
(assembly root
  ; The last two of twenty, which cross the loop's start.
  (for i (range 18 20) (instance link-$i link (path-at P (+ (* i (/ (path-length P) 20)) (* 10 t))))))
(assert links-apart (disjoint (shape-of root/link-18) (shape-of root/link-19)))

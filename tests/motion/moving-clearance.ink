; A clearance read from a solve that moves: q = 100t outgrows the gap
; between the balls (8 - t), so the clearance fails once the gap is short.
(import "std")
(t-span 0 1)
(solve distance (unknowns (q 0)) (== q (* 100 t)))
(body ball (shape (sphere 1)))
(assembly root (instance a ball (at (move t 0 0))) (instance b ball (at (move 10 0 0))))
(assert clear (clearance (shape-of root/a) (shape-of root/b) q))
; expect-status: clear fails estimate

; A promise is made for its constants as computed (docs/IR.md, The one
; promise). This cone's length, 1e308 - -1e308, overflows: norm2-up of it
; was the smallest float (inf / inf is 0), the side's term had slope 1e308
; under a promise of 1, and gap held, guaranteed, between a's and b's
; points 0.042 apart. The promise is dropped now
; (not guaranteed: unknown), and cone-z's check says why.
; expect-warning: (value (- z1 z0) < inf)
; expect-status: gap unknown
(import "std")
(def a (max (- (cone-z 1 1 -1e308 1e308) 0.1) (box 4 4 4)))
(def b (move (sphere 0.01) 1.05 0 0))
(assert gap (clearance a b 0.1))

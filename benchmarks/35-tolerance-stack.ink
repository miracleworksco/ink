; Benchmark 35: a tolerance stack, proven at its worst case.
;
; Task: three blocks, each 10 +-0.1 long, stacked along x from a datum
; wall at x = 0 (each placed after the one before: explicit placements
; that read the lengths), and an end wall whose face is at G = 30.5 +-0.1.
; A tighter design moves its wall to G2 = 30.2 +-0.1. Every claim covers
; every length in every band at once (the declared tolerances), so the
; stack is checked at its worst corner, wherever that is.
;
; Hand calculation (worst case):
;   - stack length L = a + b + c in [29.7, 30.3];
;   - end wall: gap G - L in [30.4 - 30.3, 30.6 - 29.7] = [0.1, 0.9], so
;     a 0.09 clearance holds, and G - L > 0.05;
;   - tight wall: gap G2 - L in [30.1 - 30.3, 30.3 - 29.7] = [-0.2, 0.6]:
;     at a = b = c = 10.1, G2 = 30.1 the last block runs 0.2 into the
;     wall, though at the nominal sizes the gap is 0.2.
;   For comparison only (not a proof): a root-sum-square
;   stack of the four +-0.1 tolerances is +-0.2, so the tight gap of 0.2 is
;   "statistically fine" where the worst case fails.
;
; Pass criteria (ink assert, bands applied by default):
;   - wall-clear, wall-gap: proven;
;   - tight-apart and tight-gap disproven (deliberate), each witness (what
;     shows it fails) with every length in its band and a + b + c > G2
;     (>= for tight-gap);
;     tight-apart's witness point lies between G2 and a + b + c along x;
;   - with --nominal (every length at its value) tight-apart and tight-gap
;     are proven: the tolerances are what fails them.

(import "std")

(param a 10 (tol -0.1 0.1))
(param b 10 (tol -0.1 0.1))
(param c 10 (tol -0.1 0.1))
(param G 30.5 (tol -0.1 0.1))      ; the end wall's face
(param G2 30.2 (tol -0.1 0.1))     ; a tighter design's end wall

(body block-a (shape (box a 10 10)))
(body block-b (shape (box b 10 10)))
(body block-c (shape (box c 10 10)))
; The walls are one housing's: the end wall's face is G from the datum's.
(body datum (shape (move (box 2 20 20) -1 0 0)) (frame end-face (at (move G 0 0))))
(body end-wall (shape (move (box 2 20 20) 1 0 0)))
(body tight-wall (shape (move (box 2 20 20) (+ G2 1) 0 0)))

(assembly root
  (instance datum datum)
  (instance a block-a (at (move (/ a 2) 0 0)))
  (instance b block-b (at (move (+ a (/ b 2)) 0 0)))
  (instance c block-c (at (move (+ a b (/ c 2)) 0 0)))
  (instance wall end-wall (frame-of datum end-face)))

(assert wall-clear (clearance (shape-of root/c) (shape-of root/wall) 0.09))
(assert wall-gap (value (- G (+ a b c)) > 0.05))

; The tight design fails at its worst corner (and passes at the values).
; expect-status: tight-apart fails
(assert tight-apart (disjoint (shape-of root/c) (shape-of tight-wall)))
; expect-status: tight-gap fails
(assert tight-gap (value (- G2 (+ a b c)) > 0))

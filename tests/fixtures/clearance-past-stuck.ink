; A clearance no cell size decides near one block (its field, sin^2 +
; cos^2 - 1, straddles 0 throughout it), 0.1 mm off a wall, and a second
; block 0.4 mm off it, closer than the 0.5 claimed. Best first, the proof
; refines the undecided block (its bounds the lowest) and ends at its
; first stuck cell; the hunt, largest cells first, reaches the second
; block and finds the pair (it ended unknown at that cell; hunting lowest
; bounds first stayed in the undecided block).
; expect-status: c fails
(import "std")
(def zero (- (+ (square (sin x)) (square (cos x))) 1))
(body wall (shape (move (box 12 6 1) 0 0 -0.5)))
(body near (shape (min (move (max (box 4 4 1) zero) -3 0 0.6) (move (box 1 1 1) 3 0 0.9))))
(assert c (clearance (shape-of wall) (shape-of near) 0.5))

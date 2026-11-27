; A soft band (E 2000 in 200000) from x 8.1 to 16.1 in a bar pulled 0.01 mm
; (nu 0) on 0.3 mm cells: uniform 2.451 MPa (its length 8, wherever its
; faces), a safety factor of 102, exactly. Cells of one material beside the
; face at 16.1 read 5.02 MPa; the grids read the factor 49.8, 61.7 and 86.4
; on 0.3, 0.6 and 1.2 mm cells, converging away from the truth (grid error
; 35.7). The cells on its faces each its two materials in series (a
; laminate), it reads 102.0. True: band-sf-over-100. False: band-sf-under-95-9.
(import "std")
(def part (max (neg x) (- x 24) (- (abs (- y 5)) 5) (- (abs (- z 5)) 5)))
(study s (shape part) (elastic (E (* 1MPa (select (max (- 8.1 x) (- x 16.1)) 2000 200000))) (nu 0) (yield 250MPa))
 (fix a (- x 0.001) (move 0 0 0))
 (fix b (- 23.999 x) (move 0.01 0 0)) (res 0.3))
(require band-sf-under-95-9 (min-safety-factor s < 95.88 (excluding support-margin)))
(require band-sf-over-100 (min-safety-factor s > 99.96 (excluding support-margin)))

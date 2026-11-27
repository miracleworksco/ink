; A finned heat sink: a 10 W chip under an aluminium base, cooled by a small
; fan (film coefficient 25 W/(m^2 K)), the chip claimed below 60 C.
(import "std")

(param fins 7)            ; fin count (load-time: it shapes the pattern)
(param fin-t 1.6)         ; fin thickness
(param fin-h 24)          ; fin height above the base
(param base-t 4)          ; base thickness
(def size 40)

(def base (move (box size size base-t) 0 0 (* 0.5 base-t)))
; Fins stop 1 mm short of the base's edges: faces flush with the base's
; would make the joint fillet grow a bead along them.
(def fin (move (box fin-t (- size 2) fin-h) 0 0 (+ base-t (* 0.5 fin-h))))
(def span (- size (+ 2 fin-t)))       ; fins 1 mm in from the sides too
(def fin-array (pattern-linear (move fin (* -0.5 span) 0 0) fins (/ span (- fins 1)) 0 0))
(body heatsink (shape (fillet-union base fin-array 1 90)) (attr color "#B8C2CC") (attr density 2.7g/cm3))

(study cooling (shape (shape-of heatsink)) (heat (k 167W/m/K))
  (source chip (power 10W) (region (max (- (len2 x y) 8) (- z 1))))
  (convection air (face heatsink) (coefficient 25W/m2/K) (ambient 25degC))
  (res 0.8))

(require chip-cool (max-temperature cooling < 60degC))

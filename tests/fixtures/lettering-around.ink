; Text round a circle (std lettering-around-2d): read from outside at the
; top, its tops outward; at the foot, read from inside, its tops inward.
(import "std")
(def F (font "sans" (wght 700)))
(def top (lettering-around-2d (lettering "LT" F 2) 10 90 1))
(def foot (lettering-around-2d (lettering "LT" F 2) 10 270 0))
(assert top-l-stem (subset (probe -1.2 10.6 0) top))
(assert top-t-bar (subset (probe 1.3 10.9 0) top))
(assert top-t-no-foot (disjoint top (probe 1.3 9.1 0)))
(assert foot-l-stem (subset (probe -1.2 -10.6 0) foot))
(assert foot-t-bar (subset (probe 1.3 -9.1 0) foot))
(assert foot-t-no-foot (disjoint foot (probe 1.3 -10.9 0)))

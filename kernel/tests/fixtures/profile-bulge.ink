; Both profile paths keep a repeated corner's outgoing semicircle.
; Reversing a translated 1 mm profile keeps the same solid.
(import "std")
(def pts (vec (vec 0 0) (vec 0.000000001 0 1) (vec 10 0) (vec 10 10) (vec 0 10)))
(def rounded-pts (vec (vec 0 0) (vec 0.000000001 0 1) (vec 10 0) (vec 10 10) (vec 0 10 0 1)))
(require plain-has-arc (subset (probe 5 -2 0) (polygon pts)))
(require rounded-has-arc (subset (probe 5 -2 0) (polygon rounded-pts)))

(def cw (polygon (vec (vec 100000000 100000000 -1) (vec 100000000 99999999) (vec 99999999 99999999) (vec 99999999 100000000))))
(def ccw (polygon (vec (vec 99999999 100000000) (vec 99999999 99999999) (vec 100000000 99999999 1) (vec 100000000 100000000))))
(def centre (probe 99999999.5 99999999.5 0))
(require clockwise-solid (subset centre cw))
(require counterclockwise-solid (subset centre ccw))

; The existing round belongs to the corner; only its outgoing edge moves
; to the last copy. On a 2 mm round this point is outside the quarter arc.
(def rounded-first (polygon (vec (vec 0 0 0 2) (vec 0 0) (vec 10 0) (vec 10 10) (vec 0 10))))
(require corner-still-round (disjoint (probe 0.25 0.25 0) rounded-first))

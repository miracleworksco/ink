; Slopes across switches. A bar leaning at 60 degrees, its underside and
; end no more than 30 degrees past vertical (the slope of the field along z
; at least -0.496), cut square by the extrusion's ends.
; - Along the line where its underside meets its end cap, the segment's
;   clamp switches; the slope reads the clamp's rate in several terms,
;   each case of the switch bounded apart (its rate one number in each).
; - Where it meets the extrusion's end, a rounded edge (emax: outside both
;   faces, the hypot of their distances): there the slope is the faces'
;   weighted by each distance over the hypot, each weight at most 1.
; Both inside the claims' -0.707 (45 degrees).
(import "std")
(def bar2 (- (segment-distance-at here2 (vec 0 -5) (vec -20 30)) 3))
(def part (max (extrude-x bar2 -20 20) (neg z)))
(def shell (max (neg part) (- part 0.1) (- 0.1 z)))
(assert clamp-line (range (diff part z) (- (sin (deg->rad 45))) 1.01 (max shell (- (abs x) 15))))
(assert rounded-edge (range (diff part z) (- (sin (deg->rad 45))) 1.01 shell))

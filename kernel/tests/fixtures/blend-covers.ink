; What a smooth blend covers forms none of the surface: the blend lies
; below what it rounds wherever its operands are within its width, and
; check proves it there, though the blend is narrow beside the cells.
(import "std")

(def body-ball (sphere 10))
(def bump (move (sphere 0.6) 10 0 0))
; expect-warning: body 'blended': 'bump' is in its field but forms none of its surface
(body blended (shape (smooth-union body-ball bump 1)))

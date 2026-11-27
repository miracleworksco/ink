; The speaker grille's halftone as one field in time: the hole sizes blend
; from the logo into two sources' ripples and back over 4 s. Every one of
; the 2,713 holes changes size at once; nothing is rebuilt but the field.
(import "std")
(import "speaker-grille.ink" as g)

(t-span 0 4s)

(def phase (clamp (* 1.25 (- (* 0.5 (- 1 (cos (* 90deg/s t)))) 0.1)) 0 1))
(body grille-morph
  (shape (g.grille-with (mix g.logo-k g.waves phase) (mix 0.5 0.3 phase)))
  (attr color (vec 0.09 0.09 0.10))
  (attr metallic 0.4)
  (attr roughness 0.5))

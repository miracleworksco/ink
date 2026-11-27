; The lattice bracket as printed: one material, bead-blasted titanium grey,
; cut open along y = 0 so the graded gyroid core shows by its geometry
; alone: its walls thicken and its pores close toward the pin and the
; bolts, and open out far from them. The cut face is the same metal a
; shade lighter, so the section reads.
(import "std")
(import "lattice-bracket.ink" as b)

(body bracket-plain
  (shape (max b.bracket-field (neg y)))
  (attr material-note "Ti-6Al-4V, laser powder bed fusion, bead-blasted")
  (attr color (mix (vec 0.50 0.50 0.51) (vec 0.70 0.70 0.70) (clamp (- 1 (/ (abs y) 0.25)) 0 1)))
  (attr metallic 0.8)
  (attr roughness 0.5))

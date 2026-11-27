; nafems-le5
; Source for all dimensions, loads, material and printed target: https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-le5.htm
; Axial stress sigma_xx at lower free flange edge A, x=2.5 m; -108 MPa.
; Z section length 10 m, web span 2 m along y, flanges 1 m along z,
; uniform thickness.1 m. Upper flange at y 1 m, z 0..1 m; lower at
; y-1 m, z-1 m..0; web z 0. E 210 GPa, nu.3, root fully clamped.
; End flanges carry opposite z shears of .6MN, generating 1.2MNm
; aboutx as in public drawing. PointA=(2.5, -1, -1)m on lower free edge.
; True finite plates are united through their overlapping corners.
; Reference is shell warping stress; finite thickness changes junction
; and warping stiffness, of scale t/flange width=.1. No independently
; accepted sign or bound on local 3D correction; target remains ungraded.
; scoreboard: {"id":"nafems-le5","physics":"static","target":-108,"target_kind":"idealisation","units":"MPa","study":"s","qualified":false,"sources":["https://docs.software.vt.edu/abaqusv2025/English/SIMACAEBMKRefMap/simabmk-c-le5.htm"],"quantity":"Axial stress sigma_xx at lower free flange edge A, x=2.5m; -108MPa.","units_override":"MPa","qualification":"Shell warping stress, true junction and free-edge correction unqualified","probe":"stress_mpa.xx","points":[[2500,-1000,-1000]]}
(import "std")
(def part (min (move (box 10000 2000 100) 5000 0 0) (move (box 10000 100 1000) 5000 1000 500) (move (box 10000 100 1000) 5000 -1000 -500)))
(body part (shape part))
(study s (shape part)
  (elastic (E 210000MPa) (nu 0.3))
  (fix root (- x 0.000001))
  (load upper (force 0 0 600000N) (region (max (- 9999.999999 x) (- 949.999999 y))))
  (load lower (force 0 0 -600000N) (region (max (- 9999.999999 x) (+ y 949.999999)))))

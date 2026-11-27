; fusion-thermal02
; Inputs and target: https://help.autodesk.com/cloudhelp/ENU/Fusion-Simulate/files/GUID-FD915CA4-97F6-4E21-9422-E274F637CB40.htm
; 80 x 40 x 40 inches; front/back 70 F; other faces insulated;
; k=0.5 BTU/(h in F), generation 1 BTU/(h in^3). BTU is IT (1055.05585262 J).
; Asked: maximum interior temperature, 470 F = 243.3333333333333 C.
; Exact 3D solution: T(z)=70 F + q ((20 in)^2-z^2)/(2k).
; Insulated sides make the 1D solution exact for this finite 3D block.
; Two very thin temperature regions express boundary faces, not held layers.
; scoreboard: {"id":"fusion-thermal02","physics":"thermal steady","target":243.33333333333334,"target_kind":"exact","units":"C","metric":"max-temperature","study":"s","qualified":true,"vendor_error":{"Fusion":"0.00% (hand-set 7 in quadratic tet mesh)"}}
(import "std")
(def wall (box 80in 40in 40in))
(body wall (shape wall))
(study s (shape wall)
  (heat (k (* 1W/m/K (/ (* 0.5 1055.05585262 1.8) (* 3600 0.0254)))))
  (temperature front (- z -507.999999) 70degF)
  (temperature back (- 507.999999 z) 70degF)
  (source generation (power (* 1W (/ (* 1055.05585262 80 40 40) 3600))) (region -1)))
; A broad envelope requests ink's three-grid estimate; it is not a pass criterion.
(require measurement (max-temperature s < 1000000degC))

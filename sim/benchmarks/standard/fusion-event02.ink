; fusion-event02
; All loading, initial velocity, duration and published results:
; https://help.autodesk.com/cloudhelp/ENU/Fusion-Simulate/files/GUID-3134F9C7-CBC9-4D39-BB88-E265880FBF66.htm
; Initial vx=vy=3.2174 ft/s, gravity y=-32.174 ft/s^2, duration .2 s.
; Exact x travel vx*.2=7.72176 in; source prints 7.72166 in despite those
; inputs. This arithmetic conflict quarantines the published x target.
; A 1 kg cube replaces the arbitrary no-contact tetrahedron: mass and shape
; cancel exactly in translation under uniform gravity, with no rotation.
; The source explicitly states these body parameters do not affect the test.
; Horizontal travel increases over the span, so its maximum is its final value.
; A broad field claim requests maximum horizontal travel over the span.
; scoreboard: {"id":"fusion-event02","physics":"dynamics","target":196.130164,"target_kind":"exact","units":"mm","study":"trajectory","qualified":false,"qualification":"Source horizontal target7.72166in conflicts with3.2174ft/s times.2s=7.72176in; quarantined","metric":"maximum horizontal displacement","vendor_error":{"Fusion":"0.0017% (hand-set one quadratic tetrahedron; printed target arithmetic conflict)"}}
; expect-warning: benchmarks/standard/fusion-event02.ink:18: projectile.mass-density: 1000 g published mass implies density
; expect-warning: projectile.mass-density: 1000 g published mass implies density
(import "std")
(gravity 0 -9.8066352m/s2 0)
(t-span 0 0.2s)
(body projectile (shape (box 1 1 1)) (attr mass 1kg))
(def dx (simulate 0 980.66352))
(def dy (simulate 0 980.66352))
(assembly (vertical (ypos 0))
  (instance p projectile (at) (slide y ypos)))
(assembly (trajectory (xpos 0) (ypos 0))
  (instance v (vertical ypos) (at) (slide x xpos)))
(assembly root (instance shot (trajectory dx dy)))
(require measurement (value dx < 1000000))

; A compression spring under a press, whose shape compresses as the joint
; moves. The carriage slides down its guide by `travel`. The spring between
; the base and the carriage is the body (spring LENGTH), given the length
; the carriage leaves it. A body takes arguments as an assembly does, so
; the spring's shape reads the press's travel, and through it t.
; It is a closed-end compression spring: 7.62 mm OD, 0.813 mm wire, 30 mm
; free, 8 coils, 6 active. Its coils close up at 8 x 0.813 = 6.5 mm, which
; limits the carriage's travel to 23.5 mm.
; The claims, over the stroke (preloaded 1 mm, then 7 more down and back):
; the spring never enters the carriage, and its top coil stays seated
; under it (a point in the coil's end, placed from the carriage, is in the
; wire). A spring that kept its free length would fail both.
(import "std")

(t-span 0 2)
(def free 30)
(def d 0.8128)                                   ; wire
(def coils 8) (def active 6)
(def solid (* d coils))
(def r (* 0.5 (- 7.62 d)))                       ; the wire's centre radius
; The wire's centre line, length long from z = 0: a dead coil at each end,
; its turns touching (rising a wire's width a turn), the active coils
; between at the pitch the rest leaves; 32 points a turn. The wire is
; std's tube round it.
(fn compression (length)
  (let ((dead (* 0.5 (- coils active))) (pitch (/ (- length d (* 2 dead d)) active)))
    (tube (map (lambda (u) (vec (* r (cos (* 2 pi u))) (* r (sin (* 2 pi u)))
                                (+ (* 0.5 d) (* d (min u dead)) (* pitch (clamp (- u dead) 0 active)) (* d (max 0 (- u dead active))))))
               (linspace 0 coils (+ 1 (floor (+ 1 (* 32 coils))))))
          (* 0.5 d))))

; A plate with a guide rod up the spring's bore (6 mm; the rod 5 mm).
(body base
  (shape (min (move (box 40 40 6) 0 0 -3) (move (cylinder-z 2.5 40) 0 0 20)))
  (attr color (vec 0.3 0.32 0.35)))
; Its underside at z = 0, a clearance hole for the rod.
(body carriage
  (shape (subtract (move (box 30 30 6) 0 0 3) (cylinder-z 3 20)))
  (attr color (vec 0.85 0.45 0.1)))
; Seated on z = 0, its top coil's end at z = length.
(body (spring (length free))
  (shape (compression length))
  (attr color (vec 0.67 0.7 0.74)) (attr metallic 0.8) (attr roughness 0.3))

(assembly (press (travel 0))
  (instance base base)
  (instance spring (spring (- free travel)))
  (instance carriage carriage (at (move 0 0 free)) (slide z (neg travel)) (limits (- solid free) 0)))

(assembly root
  (instance p (press (+ 1 (* 3.5 (- 1 (cos (* pi t))))))))

; The seat may touch, within the claim's 0.05 mm accuracy.
(assert spring-clear-of-carriage (disjoint (shape-of root/p/spring) (shape-of root/p/carriage)) (res 0.2))
(assert spring-seated
  (subset (place (box .2 .2 .2) (at (move r 0 (* -0.5 d)) (frame-of root/p/carriage)))
          (shape-of root/p/spring)))

; A #25 roller chain on two 16-tooth sprockets, every link a part. Link i
; rides the chain's pitch path at i pitches plus how far the chain has moved
; (std's roller-chain-path and path-link-at). So the links circulate, the
; wheels turn with them, interference sees each link, and the parts list
; counts 16 inner and 16 pin links. t is seconds; the chain runs at 20 mm/s.
(import "std")

; ANSI B29.1 #25, rollerless: 6.35 pitch, 3.3 bushes 3.1 apart inside,
; plates 5.9 tall and 0.76 thick, 2.3 pins 7.9 long. A link from its first
; pin at the origin to its second at (p, 0), pins along z, the chain's
; plane z = 0.
(def p 6.35)
(fn plate-2d () (min (circle-2d 2.95) (move (circle-2d 2.95) p 0 0) (move (rect-2d p 3.245) (* 0.5 p) 0 0)))
(fn plates (z0) (warp (extrude-z (plate-2d) z0 (+ z0 0.76)) x y (abs z)))   ; a pair, z0 out from the middle
(fn at-pins (f) (min f (move f p 0 0)))
(def inner-link (subtract (min (plates 1.55) (at-pins (cylinder-z 1.65 4.64))) (at-pins (cylinder-z 1.14 8))))
(def pin-link (min (plates 2.39) (at-pins (cylinder-z 1.15 7.9))))
; A 16T sprocket (ANSI B29.1): its pins' pitch circle p / (2 sin (180 / N))
; in radius; each gap a seat 1.005 d + 0.076 across (d the bush's 3.3) on
; the pitch circle, its flanks arcs of radius p less the seat's about the
; seats either side (tangent to the seat), its first gap on +y; tips p (0.6
; + cot (180 / N)) across; 2.794 thick; a hub 18.7 across 7.94 long on +z;
; a 6 bore.
(def teeth 16)
(def step (/ (* 2 pi) teeth))
(def rp (/ p (* 2 (sin (* 0.5 step)))))
(def rs (* 0.5 (+ (* 1.005 3.3) .076)))
(def teeth-2d
  (let ((ax (- (* (cos (* 0.5 step)) x) (* (sin (* 0.5 step)) y)))
        (ay (+ (* (sin (* 0.5 step)) x) (* (cos (* 0.5 step)) y)))
        (a (- (mod (+ (atan2 ax ay) step) step) (* 0.5 step)))   ; from the nearest gap; atan2's cut at a tooth's middle
        (u (abs (* (len2 x y) (sin a)))) (v (* (len2 x y) (cos a)))
        (gap (min (- (len2 u (- v rp)) rs)
                  (max (- (- p rs) (len2 (- u (* rp (sin step))) (- v (* rp (cos step))))) (- rp v)))))
    (max (circle-2d (* 0.5 p (+ 0.6 (/ 1 (tan (* 0.5 step)))))) (neg gap))))
(def sprocket-shape
  (subtract (min (extrude-z teeth-2d -1.397 1.397) (extrude-z (circle-2d 9.35) 1.35 9.34))
            (cylinder-z 3 30)))

(def links 32)
(def center (* 8 p))                            ; equal wheels: 2 * 8 + 16 = 32 links
(def P (roller-chain-path p teeth teeth center))
(def moved (* 20 t))                            ; mm the chain has run
(t-span 0 2)

(def steel (vec .40 .43 .47))
(body inner-link (shape inner-link) (attr color steel) (attr metallic .85) (attr roughness .3))
(body pin-link (shape pin-link) (attr color (vec .65 .68 .72)) (attr metallic .85) (attr roughness .3))
(body sprocket (shape sprocket-shape)
  (attr color (vec .64 .65 .66)) (attr metallic .9) (attr roughness .24))

; Neighbouring links share a pin: each one's far pin, through the next
; one's bush, is intended overlap.
(def pin-zone (move (cylinder-z 2.15 12) p 0 0))   ; round the bush, 0.5 out
(assembly root
  (instance a sprocket (at (rot-z (roller-chain-phase-a p teeth teeth center moved))))
  (instance b sprocket (at (rot-z (roller-chain-phase-b p teeth teeth center moved)) (move center 0 0)))
  (for i (range 0 links 2) (instance link-$i inner-link (path-link-at P (+ (* i p) moved) p)))
  (for i (range 1 links 2) (instance link-$i pin-link (path-link-at P (+ (* i p) moved) p)))
  (for i (range links)
    (for j (vec (mod (+ i 1) links))
      (contact link-$i link-$j pin-zone))))

; The chain closes: its pitch path is exactly 32 pitches round.
(assert chain-closes (value (abs (- (path-length P) (* links p))) < 1e-9))

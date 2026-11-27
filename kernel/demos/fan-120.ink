; A 120 mm axial fan, to the standard frame (120 x 120 x 25 mm, mounting
; holes on 105 mm centres), turning in it. The seven blades are one blade,
; written as a warp of space: a constant-pitch helicoid, so each blade
; twists from 53 degrees at the hub to 25 at the tip, swept forward and
; cambered, thinning outward, repeated about the axis by folding the
; angle, and blended into the hub with a smooth fillet. No loft, no
; guide curves, nothing that can fail to close.
(import "std")

(param blades 7)
(param hub-r 20)       ; hub radius, mm
(param tip-r 56.5)     ; blade tip radius, mm (2 mm from the bore)
(param lead 167)       ; the helicoid's lead: axial advance per turn, mm
(param sweep 0.38)     ; forward sweep of the tip, rad
(param camber 0.006)   ; camber, rad per mm^2 of chord
(param rate 90deg/s)   ; display speed (the real one is about 1,500 rpm)

(t-span 0 4s)

(def r (len2 x y))
(def span (clamp (/ (- r hub-r) (- tip-r hub-r)) 0 1))

; The blade's mean surface: angle = z 2 pi / lead + sweep span^2 + camber z^2.
; Unrolled at radius r, the surface is a line at the pitch angle b
; (tan b = lead / 2 pi r); the point's arc a from it, folded into one
; blade's sector, gives its distance n = a sin b from the line and its
; station u = a cos b + (z - zc) / sin b along it. The blade's section is
; a capsule in (u, n), so its leading and trailing edges are round, square
; to the blade, and its axial depth is the same at every radius.
(def mean (+ (* z (/ (* 2 pi) lead)) (* sweep (square span)) (* camber (square z))))
(def sector (/ (* 2 pi) blades))
(def off (- (mod (+ (- (atan2 y x) mean) (* 0.5 sector)) sector) (* 0.5 sector)))
(def hyp (len2 lead (* 2 pi r)))
(def sin-b (/ lead hyp))
(def cos-b (/ (* 2 pi r) hyp))
(def a (* r off))
(def n (* a sin-b))
(def u (+ (* a cos-b) (/ (- z 0.75) sin-b)))
(def thick (mix 2.0 1.3 span))
(def chord (/ 15 sin-b))
(def blade (max (- (len2 (max 0 (- (abs u) (* 0.5 (- chord thick)))) n) (* 0.5 thick))
                (- r tip-r)))

(def hub (fillet-intersect (- r hub-r) (- (abs (- z 1.5)) 9.5) 3 90))
(def impeller-field (max (smooth-union hub (max blade (- (- hub-r 3) r)) 3)
                         (cylinder-z 60 30)))

; The frame: two square flanges joined by the venturi ring and the corner
; bosses, with four struts and the motor mount on the outlet side (-z).
(def bore (+ 58.5 (* 0.35 (square (max 0 (- (abs z) 9))))))
(def ring (max (- r (+ bore 1.8)) (- bore r) (- (abs z) 12.5)))
(def flange (max (extrude-z (rounded-rect-2d 120 120 5) 9 12.5) (- bore r)))
(def flanges (min flange (mirror-z flange)))
(def corners (pattern-polar (move (cylinder-z 5 25) 52.5 52.5 0) 4))
(def mount (move (cylinder-z 22 3) 0 0 -11))
(def strut (max (move (rotate-z (box 40 3.2 3) -12) 40 0 -11) (- r (+ bore 0.5))))
(def struts (pattern-polar (rotate-z strut 45) 4))
(def holes (pattern-polar (move (cylinder-z 2.15 40) 52.5 52.5 0) 4))
(def frame-field
  (subtract (smooth-union (smooth-union (min ring flanges) corners 2) (min mount struts) 1.5)
            holes))

(body frame
  (shape (max frame-field (box 121 121 26)))
  (attr material "abs")
  (attr color (vec 0.13 0.13 0.14))
  (attr roughness 0.6))

(body impeller
  (shape impeller-field)
  (attr material "abs")
  (attr color (vec 0.55 0.62 0.72))
  (attr roughness 0.5))

(assembly root
  (instance frame frame)
  (instance impeller impeller (at) (turn z (* rate t))))

; The impeller clears the frame by 0.9 mm or more as it turns: the tips
; the bore by 2 mm, the hub's foot the motor mount by about 1.
; It holds at each of the 17 moments read. Over the motion the frame is
; asked member by member (ring, flanges, corners, mount, struts: each clear
; of the turning impeller), all but the last smooth union's blend where it
; adds material, whose search reads the whole frame and ends at its limit
; of boxes; so the answer is not settled, not a yes from the moments.
; expect-status: runs-clear unknown
(require runs-clear (clearance (shape-of root/impeller) (shape-of root/frame) 0.9))

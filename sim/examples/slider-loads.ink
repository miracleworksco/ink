; A camera slider's loads against its parts' ratings: an MGN12H carriage
; pulled along its rail by a belt carries a pan axis, a shaft on two 6000
; bearings and a 51100 thrust bearing, turning a 1.5 kg camera held 30 mm
; off its axis. Each claim reads a share of a reaction (docs/DYNAMICS.md,
; Supports, Ratings).
(import "std")
(gravity 0 0 -9.81m/s2)
(t-span 0 2s)

; The parts, from their makers' drawings, their ratings from their
; catalogues (slider-loads.json): an MGNR12 rail 400 long (12 wide, 8 tall,
; its bottom at z = 0, along +x from x = 0); its MGN12H carriage (27 wide,
; 45.4 long, from 3 to 13 above the rail's bottom, about its rail at x = 0);
; 6000 bearings (10 x 26 x 8) and a 51100 thrust bearing (10 x 24 x 9),
; each about its centre. A support's record is its (attr spec ...); its
; load frame is where its rating acts (a bearing's centre, z along its
; bore; the carriage's x along its rail).
(def records (read-json "slider-loads.json"))
; ISO 76 static load and ISO 281 basic rating life as derived facts on a
; bearing body, read by its ratings (docs/DYNAMICS.md, Ratings).
(fn bearing-derives (spec axial)
  (let ((r (get spec "load-ratings"))
        (c (* 1N (get r "C-N"))) (c0 (* 1N (get r "C0-N"))) (p (get r "life-exponent")))
    (vec
      (vec "static-load-rating" (lambda () c0))
      (vec "static-equivalent-load"
        (lambda (radial-load axial-load)
          (let ((fr (abs radial-load)) (fa (abs axial-load)))
            (if axial
              (vec fa (<= fr (* 1e-12 fa)) "radial force on a thrust bearing")
              (max fr (+ (* .6 fr) (* .5 fa))))))
        1N 1N)
      (vec "equivalent-load"
        (lambda (radial-load axial-load)
          (let ((fr (abs radial-load)) (fa (abs axial-load)))
            (if axial
              (vec fa (<= fr (* 1e-12 fa)) "radial force on a thrust bearing")
              (vec fr (<= fa (* .19 fr)) "no combined-load factors for thrust over 0.19 of radial load"))))
        1N 1N)
      (vec "life-exponent" (lambda () p))
      (vec "rated-life"
        (lambda (duty-load mean-speed)
          (vec (/ (* 1e6 360deg (pow (/ c (max duty-load (* 1e-6 c))) p)) (abs mean-speed))
               (> (abs mean-speed) 0rpm) "at rest"))
        1N 1rpm))))
(body rail (shape (move (box 400 12 8) 200 0 4)) (attr mass 260g) (attr color (vec .7 .72 .75)))
(body carriage (shape (subtract (move (box 45.4 27 10) 0 0 8) (box 50 12.4 16.4))) (attr mass 54g) (attr color (vec .75 .77 .8))
  (attr spec (get records "carriage")) (frame load (at)))
(body radial (shape (subtract (cylinder-z 13 8) (cylinder-z 5 9))) (attr mass 18g)
  (attr spec (get records "radial")) (attr ratings (get (get records "radial") "ratings"))
  (attr derives (bearing-derives (get records "radial") 0)) (frame load (at)))
(body thrust (shape (subtract (cylinder-z 12 9) (cylinder-z 5 10))) (attr mass 19g)
  (attr spec (get records "thrust")) (attr ratings (get (get records "thrust") "ratings"))
  (attr derives (bearing-derives (get records "thrust") 1)) (frame load (at)))

(body base (shape (subtract (move (box 60 40 40) 0 0 20) (cylinder-z 13 100))) (attr material "aluminum-6061")
  (frame low (at (move 0 0 8)))                  ; the radial bearings' centres, on the pan axis
  (frame high (at (move 0 0 28)))
  (frame top (at (move 0 0 44.5))))              ; the thrust bearing's centre
(body shaft (shape (move (cylinder-z 5 50) 0 0 25)) (attr material "steel")
  (frame mount (at (move 30 0 75))))             ; the camera's centre: its underside on the shaft's end
(body camera (shape (box 60 40 50)) (attr mass 1500g))

(assembly (slider (s 200) (pan 0))
  (instance rail rail)
  (instance c carriage (frame-of rail) (slide x s) (effort-limits -60N 60N))   ; what the belt may pull, N
  (instance base base (at (move 0 0 13) (frame-of c)))
  (instance ba radial (frame-of base low) (fact required-life 10000h))
  (instance bb radial (frame-of base high) (fact required-life 10000h))
  (instance th thrust (frame-of base top) (fact required-life 10000h))
  (instance pan shaft (frame-of base) (turn z pan) (effort-limits -2000N*mm 2000N*mm))
  (instance camera camera (frame-of pan mount))
  (supports pan ba bb th))
(assembly root (instance s (slider (+ 200 (* 150 (sin (* 180deg/s t)))) (* 90deg (sin (* 180deg/s t))))))

; The camera turns 1 mm over the thrust bearing's shaft washer.
(require camera-over-thrust (clearance (shape-of root/s/th) (shape-of root/s/camera) 0.5))

; The carriage: the whole reaction but the belt's pull along the rail.
(def carriage (nth (reaction root/s/c (by root/s/c)) 0))
(require carriage-static (value (carriage-static-ratio carriage (get records "carriage")) <= 0.5))   ; fSL, fSM >= 2

; The belt: its pull on the carriage, on 30 N of pretension.
(require belt-tension (value (belt-tension-ratio (belt-tension (nth (reaction root/s/c) 0) 30N) (get records "belt")) <= 1))

; The pan axis's bearings: the radial bearings share its tilting moment
; over their 20 mm, the thrust bearing its weight.
(require pan-bearings (within-ratings root/s (excluding root/s/c root/s/pan)))

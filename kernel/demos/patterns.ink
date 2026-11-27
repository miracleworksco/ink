; One part, twelve fills: a 40 mm lattice coupon (the cube a lattice is
; tested in) filled twelve ways. Every fill is one field over the same
; cube, so any two can be crossfaded: (param a) and (param b) pick two of
; them and (param w) blends a into b (0: all a, 1: all b). A blend of two
; fields is a field, so the morph between, say, a gyroid and a honeycomb
; is a real solid at every w, with its own walls and tunnels.
;   ink mesh patterns.ink                             gyroid half into Schwarz P
;   ink mesh patterns.ink --param a=3 --param w=0     the graded gyroid alone
(import "std")

(param a 0)            ; the fill shown (0 to 11, as in `fills` below)
(param b 1)            ; the fill it blends into
(param w 0.5)          ; how far: 0 all a, 1 all b; by default the morph
                       ; halfway, where both fills form the walls

(def side 40)          ; the coupon, mm
(def h (* 0.5 side))
(def cube (move (box side side side) 0 0 h))
(def zc (- z h))       ; height from the coupon's centre

; ---------------------------------------------------------------- TPMS
; Triply periodic minimal surfaces, each divided by its greatest gradient
; (sqrt(3) k for all three, k = 2 pi / period), so a shell of one is at
; least as thick as asked.
(def period 10)
(def k (/ (* 2 pi) period))
(def norm (* k (sqrt 3)))
(def gyr (gyroid period))
(def schwarz-p (/ (+ (cos (* k x)) (cos (* k y)) (cos (* k z))) norm))
(def diamond
  (let ((sx (sin (* k x))) (sy (sin (* k y))) (sz (sin (* k z)))
        (cx (cos (* k x))) (cy (cos (* k y))) (cz (cos (* k z))))
    (/ (+ (* sx sy sz) (* sx cy cz) (* cx sy cz) (* cx cy sz)) norm)))

; A sheet th thick about a surface's zero: |F| under th/2. Written here,
; not as std's shell (an offset, which takes a distance and a fixed
; thickness): F is not a distance, and a graded sheet's th is a field.
(fn sheet-of (f th) (- (abs f) (* 0.5 th)))
(def fill-gyroid (sheet-of gyr 1.0))
(def fill-schwarz (sheet-of schwarz-p 1.0))
(def fill-diamond (sheet-of diamond 0.9))
; Graded: the same gyroid, its wall thickness a field: 75 % solid at the
; base falling to 12 % at the top (a heavier lattice where the load is).
(def rho (mix 0.75 0.12 (clamp (/ z side) 0 1)))
(def fill-graded (sheet-of gyr (gyroid-thickness rho period)))

; ---------------------------------------------------------------- struts
; Octet truss (struts joining each node of a face-centred cubic lattice to
; its 12 neighbours), 10 mm cells, struts 0.8 mm in radius. Folded into
; one octant of a cell, the cell's struts are six segments.
(fn seg (px py pz ax ay az bx by bz)
  (let ((ux (- bx ax)) (uy (- by ay)) (uz (- bz az))
        (s (clamp (/ (+ (* (- px ax) ux) (* (- py ay) uy) (* (- pz az) uz))
                     (+ (square ux) (square uy) (square uz))) 0 1)))
    (len3 (- px (+ ax (* s ux))) (- py (+ ay (* s uy))) (- pz (+ az (* s uz))))))
(def fill-octet
  (let ((c period) (hc (* 0.5 period))
        (qx (abs (- (mod x c) hc))) (qy (abs (- (mod y c) hc))) (qz (abs (- (mod z c) hc))))
    (- (min (seg qx qy qz 0 0 0 hc hc 0) (seg qx qy qz 0 0 0 hc 0 hc) (seg qx qy qz 0 0 0 0 hc hc)
            (seg qx qy qz hc hc 0 hc 0 hc) (seg qx qy qz hc hc 0 0 hc hc) (seg qx qy qz hc 0 hc 0 hc hc))
       0.8)))

; ---------------------------------------------------------------- 2D cells
; Honeycomb: hexagonal cells 6 mm across flats, 0.7 mm walls, on a 1 mm
; skin, the core machined to a dome 38 mm high at the centre and 22 at the
; corners (a contoured sandwich core, as under a radome or a wing skin).
; The nearest cell centre is the nearer of two rectangular lattices, the
; second offset by half a cell.
(def hex-pitch 6)
(def hex-py (* hex-pitch (sqrt 3)))
(fn hex-local (ox oy)
  (let ((cx (- (* hex-pitch (floor (+ 0.5 (/ (+ x ox) hex-pitch)))) ox))
        (cy (- (* hex-py (floor (+ 0.5 (/ (+ y oy) hex-py)))) oy)))
    (vec (- x cx) (- y cy))))
(fn hex-cell-distance (p) (warp (hexagon-2d hex-pitch) (nth p 0) (nth p 1) 0))
(def hex-a (hex-local 0 0))
(def hex-b (hex-local (* 0.5 hex-pitch) (* 0.5 hex-py)))
(def hex-near (select (- (len2 (nth hex-a 0) (nth hex-a 1)) (len2 (nth hex-b 0) (nth hex-b 1))) hex-a hex-b))
(def hex-dome (- z (- 38 (* 0.02 (+ (square x) (square y))))))
(def fill-honeycomb (max (min (- (abs (hex-cell-distance hex-near)) 0.35) (- z 1)) hex-dome))

; Halftone perforation: a hollow cube, 1.6 mm walls, its faces
; perforated on a 2.6 mm hex grid, each hole's diameter a field of where it
; is on its face (0.6 to 2.1 mm, in rings 12 mm apart): a vent or speaker
; housing whose open area is designed. One hole field per axis, each kept
; to the faces square to it.
(def ht-pitch 2.6)
(def ht-py (* ht-pitch (sqrt 3)))
(fn ht-holes (u v)
  (let ((r (+ 0.3 (* 0.75 (+ 0.5 (* 0.5 (cos (/ (* 2 pi (len2 u v)) 12)))))))
        (hole (lambda (ou ov)
          (let ((cu (- (* ht-pitch (floor (+ 0.5 (/ (+ u ou) ht-pitch)))) ou))
                (cv (- (* ht-py (floor (+ 0.5 (/ (+ v ov) ht-py)))) ov)))
            (- (len2 (- u cu) (- v cv)) r)))))
    (max (min (hole 0 0) (hole (* 0.5 ht-pitch) (* 0.5 ht-py))) (- (max (abs u) (abs v)) 17.5))))
(def fill-halftone
  (max (- (abs (warp (box (- side 1.6) (- side 1.6) (- side 1.6)) x y zc)) 0.8)
       (neg (min (ht-holes x y) (ht-holes y zc) (ht-holes x zc)))))

; Phyllotaxis: a sunflower pin-fin heatsink, 170 pins 1.7 mm across on a
; 3 mm disc, 32 mm tall at the centre falling to 19 at the rim, at the
; seeds of a sunflower (Vogel's model: seed n at radius c sqrt(n), turned
; n golden angles). The pins are one union of 170 distances, the seeds
; computed as the file loads.
(def phy-c 1.45)
(def golden (* pi (- 3 (sqrt 5))))
(def phy-d (reduce min (map (lambda (n)
  (len2 (- x (* phy-c (sqrt n) (cos (* golden n)))) (- y (* phy-c (sqrt n) (sin (* golden n))))))
  (range 1 171))))
(def phy-pins (max (- phy-d 0.85) (- z (- 33 (* 0.036 (square (len2 x y)))))))
(def phy-disc (max (- (len2 x y) 19.8) (- z 3)))
(def fill-phyllotaxis (smooth-union phy-pins phy-disc 1.2))

; ---------------------------------------------------------------- 3D cells
; Voronoi foam: the walls between the cells of 125 seeds jittered off an
; 8 mm grid (a hash of each cell's index places its seed), 0.8 mm thick.
; F1 and F2 are the distances to the nearest and second-nearest seed over
; the 27 neighbouring cells; the wall is where they are nearly equal.
(def vc 8)
(fn fract (v) (- v (floor v)))
(fn hash (i j kk s) (fract (* 43758.5453 (sin (+ (* 127.1 i) (* 311.7 j) (* 74.7 kk) s)))))
(def vi (floor (/ x vc)))
(def vj (floor (/ y vc)))
(def vk (floor (/ z vc)))
(def vds
  (flatten (map (lambda (di) (flatten (map (lambda (dj) (map (lambda (dk)
    (let ((i (+ vi di)) (j (+ vj dj)) (kk (+ vk dk)))
      (len3 (- x (* vc (+ i 0.5 (* 0.8 (- (hash i j kk 0) 0.5)))))
            (- y (* vc (+ j 0.5 (* 0.8 (- (hash i j kk 17.3) 0.5)))))
            (- z (* vc (+ kk 0.5 (* 0.8 (- (hash i j kk 41.9) 0.5))))))))
    (range -1 2))) (range -1 2)))) (range -1 2))))
(def vf1 (reduce min vds))
(def vf2 (reduce min (map (lambda (d) (select (- vf1 d) d 1e3)) vds)))
(def fill-voronoi (- (* 0.5 (- vf2 vf1)) 0.4))

; ---------------------------------------------------------------- the rest
; Menger sponge, three levels: the cube, each third's centre cross bored
; out, again in each of the 20 cubes left, and again (holes down to
; 1.5 mm). Self-similar, and three folds of space.
(def fill-menger
  (let ((px (/ x h)) (py (/ y h)) (pz (/ zc h)))
    (* h (apply max (map (lambda (m)
      (let ((s (pow 3 m))
            (rx (abs (- 1 (* 3 (abs (- (mod (* px s) 2) 1))))))
            (ry (abs (- 1 (* 3 (abs (- (mod (* py s) 2) 1))))))
            (rz (abs (- 1 (* 3 (abs (- (mod (* pz s) 2) 1)))))))
        (/ (- (min (max rx ry) (max ry rz) (max rz rx)) 1) (* 3 s))))
      (range 3))))))

; Diamond knurl on a 36 mm knob, the coupon's height: two families of helical grooves at +-30
; degrees, 72 around, 0.5 mm deep, the ends chamfered.
(def kn-r 18)
(def kn-n 72)
(def kn-turn (/ (* kn-n (atan2 y x)) (* 2 pi)))
(def kn-lead (/ z (* (/ (* 2 pi kn-r) kn-n) (sqrt 3))))
(fn tri (s) (abs (- (fract s) 0.5)))
(def kn-depth (* 0.5 2 (min (tri (+ kn-turn kn-lead)) (tri (- kn-turn kn-lead)))))
(def fill-knurl
  (max (- (len2 x y) (- kn-r 0.5 (neg kn-depth)))
       (- (+ (len2 x y) (abs zc)) (+ kn-r 18.5))))

; A woven block: eight plain-weave fabrics stacked 5 mm apart, warp along
; x and weft along y, each strand 1.15 mm in radius passing over one and
; under the next, 5 mm apart (a 3D-woven composite preform).
(def wv-p 5)
(def wv-a 0.9)
(def wv-r 1.15)
(def wv-layer (* wv-p (floor (+ 0.5 (/ zc wv-p)))))
(def wv-zl (- zc wv-layer))
(def wv-j (floor (+ 0.5 (/ y wv-p))))
(def wv-i (floor (/ x wv-p)))
(def wv-warp (len2 (- y (* wv-p wv-j)) (- wv-zl (* wv-a (sin (/ (* pi x) wv-p)) (cos (* pi wv-j))))))
(def wv-weft (len2 (- x (* wv-p (+ wv-i 0.5))) (+ wv-zl (* wv-a (cos (/ (* pi y) wv-p)) (cos (* pi wv-i))))))
(def fill-weave (* 0.8 (- (min wv-warp wv-weft) wv-r)))

; ---------------------------------------------------------------- all twelve
(def fills (vec fill-gyroid fill-schwarz fill-diamond fill-graded
                fill-octet fill-honeycomb fill-voronoi fill-halftone
                fill-phyllotaxis fill-menger fill-knurl fill-weave))
(def names (vec "gyroid" "schwarz-p" "diamond" "graded" "octet" "honeycomb"
                "voronoi" "halftone" "phyllotaxis" "menger" "knurl" "weave"))
(def colors (vec (vec 0.95 0.42 0.12) (vec 0.13 0.55 0.92) (vec 0.55 0.35 0.85) (vec 0.9 0.2 0.3)
                 (vec 0.2 0.7 0.55) (vec 0.95 0.72 0.15) (vec 0.85 0.4 0.6) (vec 0.55 0.6 0.68)
                 (vec 0.8 0.5 0.2) (vec 0.6 0.45 0.3) (vec 0.72 0.74 0.77) (vec 0.15 0.45 0.65)))

(def fill (blend (nth fills a) (nth fills b) w))
(def paint (+ (* (- 1 w) (nth colors a)) (* w (nth colors b))))

(body coupon
  (shape (max cube fill))
  (attr color paint)
  (attr metallic 0.25)
  (attr roughness 0.45))

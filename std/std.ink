; Implicit Kernel standard library.
;
; Plain IR, math over x y z t: the kernel never sees these names, only
; arithmetic, min and max. The language reads three std layouts (a path's
; and `carried`'s, a loop's closing gap from path-closes, and materials.ink's
; (data ...) records), and its hints name std's functions. Read it, copy it,
; replace it.
;
; Conventions:
;   - A shape is a field < 0 inside, > 0 outside, 0 on the surface.
;   - Lengths are mm. Angles are degrees, except path-arc's sweep and the
;     path tools' angles (radians).
;   - "Exact": the field is the true Euclidean distance. Others bound it or
;     are right only on the surface. Clearance and meshing work on them;
;     offsets need a distance field. Thickness is sampled (measured).
;   - Exact primitives promise a slope of at most 1 (1-Lipschitz) with
;     (lipschitz 1 f), which clearance proofs trust in place of interval
;     arithmetic's looser slope bound. Warps keep the promise
;     (scales rescale it), as do min, max and adding a constant.
;   - The promise holds only for fixed sizes: a size that reads x, y or z
;     (a radius varying along the part, a sweep's section at s) drops it
;     (`ink spec --section "The one promise"`), and so does a size that
;     overflows or is undefined (inf / inf), where it does.

; ---------------------------------------------------------------- scalar helpers

(fn clamp (v lo hi) (min (max v lo) hi))
; n numbers evenly from a to b, both ends included (n a constant >= 2).
(fn linspace (a b n)
  (if (< n 2) (error "linspace takes n >= 2, got " n)
    (map (lambda (i) (+ a (* (- b a) (/ i (- n 1))))) (range n))))
; The gradient of a field, (vec df/dx df/dy df/dz).
(fn grad (f) (vec (diff f x) (diff f y) (diff f z)))
(fn mix (a b k) (+ a (* (- b a) k)))
; mix for k in [0, 1], clamped to the range of a and b. The value is the
; same, but interval bounds see the clamp, so a blend of bounded shapes
; stays bounded far out.
(fn blend (a b k) (min (max (mix a b k) (min a b)) (max a b)))
; A DC motor's (effort-limits LO HI) at full voltage and speed `rate`: the
; line from `stall` torque at rest to none at `free` speed, either way,
; capped at current `amps` (stall draws `stall-amps`). rate and free share
; a unit: deg/s, or mm/s for a linear drive (torques N mm, forces N).
(fn dc-torque-max (stall free stall-amps amps rate) (min (* stall (- 1 (/ rate free))) (* stall (/ amps stall-amps))))
(fn dc-torque-min (stall free stall-amps amps rate) (max (neg (* stall (+ 1 (/ rate free)))) (neg (* stall (/ amps stall-amps)))))
; Where a support's force acts in its plane: the zero-moment point, (vec x
; y) in mm, about which the moment is only a twist about z. r is the
; support's reaction, (vec fx fy fz mx my mz), in a frame on its plane with
; z into what stands on it. Where fz <= 0 what stands there lifts off and
; the point is at infinity, so a claim that it lies within a footprint fails.
(fn zmp (r)
  (let ((fz (nth r 2)))
    (vec (select (neg fz) (/ (neg (nth r 4)) fz) inf) (select (neg fz) (/ (nth r 3) fz) inf))))
; The piecewise-linear function through `pts` ((x y) ..., x ascending) at
; s, flat past the ends.
(fn curve-at (pts s)
  (fold + (nth (nth pts 0) 1)
    (map (lambda (a b) (* (/ (- (nth b 1) (nth a 1)) (- (nth b 0) (nth a 0)))
                          (- (clamp s (nth a 0) (nth b 0)) (nth a 0))))
         (slice pts 0 -1) (slice pts 1))))
(fn len2 (a b) (sqrt (+ (square a) (square b))))
(fn len3 (a b c) (sqrt (+ (square a) (square b) (square c))))
; How near a support is to sliding, with r as for zmp: |(fx, fy)| / fz,
; the force across its plane over the force into it. What stands there
; slides once this passes the support's friction coefficient, which the
; model states: (value (slip r) < 0.8). Infinite where fz <= 0 (lift-off).
(fn slip (r)
  (let ((fz (nth r 2)))
    (select (neg fz) (/ (len2 (nth r 0) (nth r 1)) fz) inf)))
(fn deg->rad (d) (* d (/ pi 180)))
(fn rad->deg (a) (* a (/ 180 pi)))
; Inverse tangent, cosine and sine, in radians; acos and asin take |v| <= 1.
(fn atan (v) (atan2 v 1))
(fn acos (v) (atan2 (sqrt (- 1 (square v))) v))
(fn asin (v) (atan2 v (sqrt (- 1 (square v)))))
(fn mod (a m) (- a (* m (floor (/ a m)))))
(fn sign (a) (select a -1 1))

; Synchronous belt pulleys, with grooves derived from the belt tooth: how
; far they are from a maker's groove is unknown.
; GT: (pitch PLD height back root-radius flank-radius crown-radius offset),
; pitch 2 or 3, where PLD (pitch line differential) is pitch radius less tip
; radius. From Unitta/MISUMI belt drawings, NOT a released pulley
; generating profile.
(fn gt-profile (pitch)
  (nth (vec (vec 2 .254 .75 .63 .15 1 .555 .4)
            (vec 3 .381 1.14 1.26 .25 1.52 .85 .61))
    (if (== pitch 2) 0 (if (== pitch 3) 1 2))))
; A GT tooth's right half as three arcs, (end-x end-y signed-radius
; centre-x centre-y) each, keeping the drawn tooth height.
(fn gt-tooth-half (pitch)
  (let ((p (gt-profile pitch)) (h (nth p 2)) (rr (nth p 4))
        (rf (nth p 5)) (rt (nth p 6)) (b (nth p 7))
        (cv (- h rt)) (d (sqrt (- (square (- rf rt)) (square b))))
        (fv (- cv d)) (ru (+ (neg b) (sqrt (- (square (+ rf rr)) (square (- rr fv)))))))
    (vec (vec (+ (neg b) (/ (* rf b) (- rf rt))) (+ fv (/ (* rf d) (- rf rt))) rt 0 cv)
         (vec (+ (neg b) (/ (* rf (+ ru b)) (+ rf rr))) (+ fv (/ (* rf (- rr fv)) (+ rf rr))) rf (neg b) fv)
         (vec ru 0 (neg rr) ru rr))))
; A tooth's 2D field from exact line and circle distances, not sampled.
; Its overlap into the belt's backing (y >= -back/2) is hidden there.
(fn timing-tooth-2d (segments height back)
  (let ((u (abs x))
        (ends (concat (vec (vec 0 height)) (map (lambda (g) (slice g 0 2)) segments)))
        (side (lambda (k)
          (let ((g (nth segments k)) (r (nth g 2)) (cu (nth g 3)) (d (len2 (- u cu) (- y (nth g 4)))))
            (if (> r 0) (if (> cu 0) (min (- d r) (- u cu)) (- d r))
              (if (< r 0) (max (- (neg r) d) (- u cu))
                (let ((a (nth ends k)) (e (- (nth ends (+ k 1)) a)))
                  (/ (- (* (- y (nth a 1)) (nth e 0)) (* (- u (nth a 0)) (nth e 1))) (norm2-up (nth e 0) (nth e 1)))))))))
        (above (fold (lambda (f k) (let ((h (nth (nth ends k) 1))) (min (max f (- h y)) (max (side k) (- y h)))))
                     (side 0) (range 1 (len segments))))
        (u0 (nth (nth ends -1) 0)))
    (min (max above (neg y)) (max y (- u u0) (- (* -.5 back) y)))))
(fn gt-tooth-2d (pitch)
  (let ((p (gt-profile pitch))) (timing-tooth-2d (gt-tooth-half pitch) (nth p 2) (nth p 3))))
; n `groove`s (x across, y inward from the tip) cut by rigid sector folds
; into radius `tip`, not past `root`, turned by `turn`, with lands (the rim
; between grooves) centred on +x. No slope is promised, and a fractional n
; is not rounded.
(fn timing-teeth-2d (n tip root groove turn)
  (let ((r (len2 x y)) (sector (/ (* 2 pi) n))
        (a (- (atan2 y x) (deg->rad turn)))
        (phase (* n a))
        (fractional (- (mod (+ a pi) sector) (* .5 sector)))
        (g (if (== n (floor n))
             (/ (+ (* .5 pi) (atan2 (cos phase) (abs (sin phase)))) n)
             (abs fractional)))
        (side (if (== n (floor n)) (sin phase) fractional))
        (u (* r (sin g))) (v (- tip (* r (cos g))))
        (cut (select side (warp groove (neg u) v z) (warp groove u v z))))
    (min (max (tag tooth-tips (- r tip)) (neg cut)) (- r root))))
; A groove drawn in the pulley's x, y at `angle`, moved rigidly into ISO/GT
; groove coordinates: x clockwise across, y inward from `tip`.
(fn timing-radial-groove-2d (tip angle groove)
  (let ((a (deg->rad angle)) (c (cos a)) (s (sin a)) (r (- tip y)))
    (warp groove (+ (* r c) (* x s)) (- (* r s) (* x c)) z)))
; GT pulley teeth, pitch 2 or 3. No clearance added: the separately drawn
; belt tooth's outline, placed by its PLD. Deviation from a maker's groove
; is UNKNOWN.
(fn gt-pulley-teeth-2d (n pitch turn)
  (let ((p (gt-profile pitch)) (tip (- (/ (* n pitch) (* 2 pi)) (nth p 1))))
    ; 1-Lipschitz for whole n only: other n do not close the ring.
    (lipschitz (if (> n 0) (if (== n (floor n)) 1 0) 0)
      (timing-teeth-2d n tip (- tip (nth p 2)) (gt-tooth-2d pitch) turn))))

; ISO 13050:2022 H Table 14: (n-min n-max depth X R1 phi rt), pitch 3, 5
; or 8, one row per tooth-count range; (vec) where no row covers n.
(fn htd-profile (pitch n)
  (let ((rows (nth (vec
    (vec (vec 10 13 1.190 .029 .991 15 .181) (vec 14 25 1.179 .112 .889 9 .229)
         (vec 26 80 1.219 .028 .927 8 .191) (vec 81 200 1.234 .074 .925 4 .301))
    (vec (vec 12 16 1.989 .307 1.265 10 .432) (vec 17 25 2.009 .320 1.270 6 .508)
         (vec 26 80 2.052 .081 1.438 2 .488) (vec 81 200 2.056 .028 1.552 5 .569))
    (vec (vec 22 27 3.295 0 2.675 11.3 .874) (vec 28 89 3.604 0 2.629 7 1.024)
         (vec 90 200 3.630 0 2.639 6.6 1.008)))
     (if (== pitch 3) 0 (if (== pitch 5) 1 (if (== pitch 8) 2 3))))))
    (fold (lambda (found row) (if (>= n (nth row 0)) (if (<= n (nth row 1)) row found) found)) (vec) rows)))
; ISO H Table 9 belt pitch-line differentials (not an added running gap).
(fn htd-pld (pitch) (nth (vec .381 .572 .686) (if (== pitch 3) 0 (if (== pitch 5) 1 (if (== pitch 8) 2 3)))))
; polygon's (x y bulge) at p: the minor arc about c to q.
(fn timing-arc (p q c)
  (let ((a (- p c)) (b (- q c)))
    (vec (nth p 0) (nth p 1)
      (tan (/ (atan2 (- (* (nth a 0) (nth b 1)) (* (nth a 1) (nth b 0))) (dot a b)) 4)))))
; Points for polygon: a stand-in, from Table 14, for missing Figure 11.
; X is half the bottom flat, R1's centre is (X, depth - R1), and phi is the
; straight flank's angle to the radial axis. The flank touches R1, and rt
; touches the flank and the actual tip circle. All joins are tangent, the
; bottom flat's too. Not ISO cutter certification.
(fn htd-groove-polygon (pitch n tip)
  (let ((p (htd-profile pitch n)) (h (nth p 2)) (xx (nth p 3))
        (r (nth p 4)) (a (deg->rad (nth p 5))) (rt (nth p 6))
        (ca (cos a)) (sa (sin a)) (cv (- h r))
        (d (- (+ r rt (* ca xx) (* sa cv)) (* sa tip)))
        (q (sqrt (- (square (- tip rt)) (square d))))
        (c (vec (+ (* d ca) (* q sa)) (+ tip (* d sa) (neg (* q ca)))))
        (lip (- c (vec (* rt ca) (* rt sa))))
        (top (+ (vec 0 tip) (* (- c (vec 0 tip)) (/ tip (- tip rt)))))
        (flank (vec (+ xx (* r ca)) (+ cv (* r sa))))
        (bottom (vec xx h)) (mirror (lambda (v) (vec (neg (nth v 0)) (nth v 1)))))
    (concat
      (vec (timing-arc top lip c) lip (timing-arc flank bottom (vec xx cv)))
      (if (> xx 0) (vec bottom) (vec))
      (vec (timing-arc (mirror bottom) (mirror flank) (vec (neg xx) cv))
           (mirror flank) (timing-arc (mirror lip) (mirror top) (mirror c))
           (mirror top) (vec (neg (nth top 0)) (neg pitch)) (vec (nth top 0) (neg pitch))))))
; HTD teeth at an explicit tip radius, for ODs (outside diameters) a
; source gives. The groove is its polygon read at |x|, so mirrored halves
; meet exactly. The slope promise is 1 + 4e-15, for whole n.
(fn htd-pulley-teeth-at-tip-2d (n pitch tip turn)
  (lipschitz (if (> n 0) (if (== n (floor n)) (+ 1 4e-15) 0) 0)
    (timing-teeth-2d n tip (- tip (nth (htd-profile pitch n) 2))
      (warp (polygon (htd-groove-polygon pitch n tip)) (abs x) y z) turn)))
; The tip radius: half the OD in Pfeifer's published 8 mm HTD table, n
; 22..200. It keeps the table's rounding and small-count corrections, not
; OD = PD (pitch diameter) - 2 PLD.
(fn htd8-tip-radius (n)
  (* .5 (nth (vec
    54.65 57.2 59.74 62.29 64.84 67.38 70.08 72.62 75.13 77.65 80.16 82.68 85.21 87.76
    90.3 92.85 95.39 97.94 100.49 103.03 105.58 108.13 110.67 113.22 115.77 118.31 120.86 123.41
    125.95 128.5 131.05 133.59 136.14 138.68 141.23 143.78 146.32 148.87 151.42 153.96 156.51 159.06
    161.6 164.15 166.7 169.24 171.79 174.34 176.88 179.43 181.97 184.52 187.07 189.61 192.16 194.71
    197.25 199.8 202.35 204.89 207.44 209.99 212.53 215.08 217.63 220.17 222.72 225.27 227.81 230.36
    232.9 235.45 238.0 240.54 243.09 245.64 248.18 250.73 253.28 255.82 258.37 260.92 263.46 266.01
    268.56 271.1 273.65 276.19 278.74 281.29 283.83 286.38 288.93 291.47 294.02 296.57 299.11 301.66
    304.21 306.75 309.3 311.85 314.39 316.94 319.48 322.03 324.58 327.12 329.67 332.22 334.76 337.31
    339.86 342.4 344.95 347.5 350.04 352.59 355.14 357.68 360.23 362.77 365.32 367.87 370.41 372.96
    375.51 378.05 380.6 383.15 385.69 388.24 390.79 393.33 395.88 398.43 400.97 403.52 406.07 408.61
    411.16 413.7 416.25 418.8 421.34 423.89 426.44 428.98 431.53 434.08 436.62 439.17 441.72 444.26
    446.81 449.36 451.9 454.45 456.99 459.54 462.09 464.63 467.18 469.73 472.27 474.82 477.37 479.91
    482.46 485.01 487.55 490.1 492.65 495.19 497.74 500.28 502.83 505.38 507.92
  ) (- n 22))))
(fn htd-pulley-teeth-2d (n pitch turn)
  (let ((tip (if (== pitch 8) (htd8-tip-radius n)
                (- (/ (* n pitch) (* 2 pi)) (htd-pld pitch)))))
    (htd-pulley-teeth-at-tip-2d n pitch tip turn)))

; ---------------------------------------------------------------- combination

; Union of any number of shapes: (min a b c ...). Intersection:
; (max a b c ...). f grown by d (shrunk for d < 0): (- f d), an exact
; distance only if f is one. (subtract a b c ...) is a with b, c, ... cut
; away.
(fn subtract (a b more ...) (apply max (concat (vec a (neg b)) (map neg more))))

; Smooth union of a and b, blended over k mm (Quilez's polynomial smooth
; minimum). It lies between min(a, b) - k/4 and min(a, b). Where a and b
; are within k of each other the blend is its own construction, named
; "blend". Elsewhere it is min(a, b), whose faces are a's and b's.
(fn smooth-union (a b k)
  (let ((h (clamp (+ 0.5 (/ (* 0.5 (- b a)) k)) 0 1))
        (m (min a b)))
    (select (- (abs (- a b)) k)
      (tag "blend" (min m (max (- m (* 0.25 k)) (- (mix b a h) (* k h (- 1 h))))) opaque)
      m)))
(fn smooth-intersect (a b k) (neg (smooth-union (neg a) (neg b) k)))
(fn smooth-subtract (a b k) (neg (smooth-union (neg a) b k)))

; Concave fillet of radius r filling an inside corner; union it with the
; part. u and v are distances from the two walls into the open corner
; (each > 0 away from its wall), e.g. u = (- y 5) for a wall face at
; y = 5. Unlike smooth-union it adds material only in the corner, never on
; coplanar faces.
(fn corner-fillet (u v r)
  (max (- u r) (- v r) (neg u) (neg v) (- r (len2 (- u r) (- v r)))))

; Fillets: a ball of radius r rolled along the edge where two shapes meet.
; - (fillet-union a b r angle) joins a and b and fills their concave edge
;   wherever the ball, rolled outside both, cannot reach (a closing).
; - (fillet-intersect a b r angle) intersects them and rounds their convex
;   edge with the ball rolled inside (an opening).
; - chamfer-union and chamfer-intersect (a b r) cut a flat bevel instead,
;   meeting each face where the other is r away.
; angle is the corner's opening in degrees: the air between the faces for
; a union, the material for an intersection. It is 180 less the angle
; between the faces' outward normals, 90 for a square corner. It may be a
; field, where the angle varies along the edge.
; Exact, inside and out, where a and b are exact distances to faces whose
; sections square to the edge are straight lines meeting at angle: planes,
; a boss on a floor, any extrusion or revolve of a straight-edged corner.
; A curved section (a rib into a round boss) is off by up to about
; r^2 / (6 R) at a square corner, R the surface's radius across the edge.
; A point is filled only where corner-reach finds no ball of radius r,
; touching one face straight behind the point, that is clear of the other
; shape. So two faces that face the same way (a part buried in another
; near that one's far face, faces coplanar along a seam) are not read as a
; corner. Facing faces of a and b closer than about 0.6 r (at angle 90)
; are bridged, as a closing does.
; The edge is its own construction, named "fillet" or "chamfer", where the
; ball is the nearest thing. Elsewhere the result is min(a, b) (max(a, b)
; for an intersection), whose faces are a's and b's. r = 0 is no edge.
(fn fillet-union (a b r angle)
  (let ((th (deg->rad (- 90 angle))) (c (neg (sin th))) (s (cos th))
        (u (- r a)) (v (- r b)) (w (/ (- v (* c u)) s)))
    (select (neg r)
      (select (neg (min (- u (* c v)) w))
        (select (neg (corner-reach a b r))
          (tag "fillet" (- (max r (min a b)) (len2 (max u (* c v)) (max w 0))) opaque)
          (min a b))
        (min a b))
      (min a b))))
(fn fillet-intersect (a b r angle) (neg (fillet-union (neg a) (neg b) r angle)))
(fn chamfer-union (a b r)
  (let ((sharp (min a b)) (cut (* (sqrt 0.5) (- (+ a b) r))))
    (select (neg r)
      (select (neg (min (- r a) (- r b)))
        (select (neg (corner-reach a b r)) (min (tag "chamfer" cut opaque) sharp) sharp)
        sharp)
      sharp)))
(fn chamfer-intersect (a b r) (neg (chamfer-union (neg a) (neg b) r)))
; > 0 where a corner of a and b with radius r may be filled at a point:
; r less the larger of a at the centre of the ball touching b's face
; straight behind the point, and b at the centre of the one touching a's.
; Where it is <= 0, one such ball is clear of both (exactly, for distance
; fields) and holds the point.
(fn corner-reach (a b r)
  (let ((ga (grad a)) (gb (grad b))
        (ka (/ (- r a) (sqrt (dot ga ga)))) (kb (/ (- r b) (sqrt (dot gb gb)))))
    (- r (max (warp a (+ x (* kb (nth gb 0))) (+ y (* kb (nth gb 1))) (+ z (* kb (nth gb 2))))
              (warp b (+ x (* ka (nth ga 0))) (+ y (* ka (nth ga 1))) (+ z (* ka (nth ga 2))))))))
;
; Morphological opening by radius r: removes features too thin for a ball
; of radius r to fit inside (thin fins, specks) and rounds convex edges to
; radius r. Approximate, as redistance is. To check wall thickness use
; (min-thickness f > t): it measures by rays and leaves sharp edges alone.
(fn opening (f r) (- (redistance (+ (redistance f) r)) r))

; Hollow a shape: a wall of the given thickness centred on its surface.
(fn shell (f thickness) (- (abs f) (* 0.5 thickness)))

; ---------------------------------------------------------------- appearance

; A body of pieces, each with its own finish: (component f color metallic
; roughness). For cs, a vec of components, paint-shape is their union and
; (paint cs i) is attribute i (1 color, 2 metallic, 3 roughness) of the
; piece the point is in, or nearest. E.g. (attr color (paint cs 1)).
(fn component (f color metallic roughness) (vec f color metallic roughness))
(fn paint-shape (cs) (reduce min (map (lambda (c) (nth c 0)) cs)))
(fn paint (cs i)
  (nth (fold (lambda (a c) (select (- (nth c 0) (nth a 0)) c a)) (nth cs 0) (slice cs 1)) i))
; A scalar field f (a height, a package's result) as an sRGB colour, from
; lo (dark blue) through cyan, green and yellow to hi (dark red), held
; beyond both. It is Turbo (Google, 2019), the usual simulation contour
; map, by its published degree-5 fit.
; (attr color (color-map z 0 50)) colours a part by height. The colour
; depends on f alone, so a viewer can draw a legend.
(fn color-map (f lo hi)
  (let ((s (clamp (/ (- f lo) (- hi lo)) 0 1)))
    (let ((poly (lambda (cs) (fold (lambda (a c) (+ (* a s) c)) 0 cs))))
      (vec (poly (vec 59.28637943 -152.94239396 132.13108234 -42.66032258 4.61539260 0.13572138))
           (poly (vec 2.82956604 4.27729857 -14.18503333 4.84296658 2.19418839 0.09140261))
           (poly (vec 27.34824973 -89.90310912 110.36276771 -60.58204836 12.64194608 0.10667330))))))

; ---------------------------------------------------------------- transforms of space (warps)

(fn move (f dx dy dz) (warp f (- x dx) (- y dy) (- z dz)))

; Turn a shape by deg about an axis through the origin, right-handed
; (counter-clockwise seen from the axis's + end), as at's rot-x/y/z do.
; 90 about x takes +y to +z, about y +z to +x, about z +x to +y.
(fn rotate-x (f deg)
  (let ((c (cos (deg->rad deg))) (s (sin (deg->rad deg))))
    (warp f x (+ (* c y) (* s z)) (- (* c z) (* s y)))))
(fn rotate-y (f deg)
  (let ((c (cos (deg->rad deg))) (s (sin (deg->rad deg))))
    (warp f (- (* c x) (* s z)) y (+ (* s x) (* c z)))))
(fn rotate-z (f deg)
  (let ((c (cos (deg->rad deg))) (s (sin (deg->rad deg))))
    (warp f (+ (* c x) (* s y)) (- (* c y) (* s x)) z)))

; Uniform scale by s about the origin. A negative s also reflects through
; the origin. s = 0 is degenerate: the field is 0 everywhere.
(fn scale (f s) (* (abs s) (warp f (/ x s) (/ y s) (/ z s))))

; Reflection across a plane through the origin (the shape flipped).
(fn mirror-x (f) (warp f (neg x) y z))
(fn mirror-y (f) (warp f x (neg y) z))
(fn mirror-z (f) (warp f x y (neg z)))
; A shape together with its reflection, the usual "mirror feature".
(fn mirror-union-x (f) (min f (mirror-x f)))
(fn mirror-union-y (f) (min f (mirror-y f)))
(fn mirror-union-z (f) (min f (mirror-z f)))
; Symmetric from the + side: the part at x >= 0 is kept and reflected onto
; x < 0, replacing what the shape had there. Exact and cheap.
; To keep both sides of an asymmetric shape use mirror-union.
(fn symmetric-x (f) (warp f (abs x) y z))
(fn symmetric-y (f) (warp f x (abs y) z))
(fn symmetric-z (f) (warp f x y (abs z)))

; Infinite repetition with period p: copies centred at 0, +-p, +-2p, ...
; The copy at the origin is the shape itself. Clip with a region to bound
; it. For n copies starting at the shape use pattern-linear.
(fn repeat-x (f p) (warp f (- (mod (+ x (* 0.5 p)) p) (* 0.5 p)) y z))
(fn repeat-y (f p) (warp f x (- (mod (+ y (* 0.5 p)) p) (* 0.5 p)) z))
(fn repeat-z (f p) (warp f x y (- (mod (+ z (* 0.5 p)) p) (* 0.5 p))))
; n copies, each (dx dy dz) on from the last, as pattern-linear makes but
; for the cost of one (a row of pins, holes or teeth). Each point is taken
; to its nearest copy along the step, copy 0 centred on the origin. So build
; f about the origin and move the row. An f off the origin along the step
; shows once, or not at all. If each copy is symmetric about its own plane
; square to the step, the field is pattern-linear's. Otherwise a copy
; beyond a point's cell is not seen.
(fn repeat-linear (f n dx dy dz)
  (let ((i (clamp (floor (+ 0.5 (/ (+ (* dx x) (* dy y) (* dz z)) (+ (square dx) (square dy) (square dz)))))
                  0 (- n 1))))
    (warp f (- x (* i dx)) (- y (* i dy)) (- z (* i dz)))))
; n copies about the z axis for the cost of one (a toothed rim from one
; tooth). Each copy is the shape in its wedge, the wedge within 180/n
; degrees of +x, turned by multiples of 360/n. The shape beyond its wedge
; is never drawn but decides distances near the wedge's edges, so a
; feature that reaches into its neighbours' wedges should include them.
; For n copies of any shape at their own places use pattern-polar.
(fn repeat-polar (f n)
  (let ((p (/ (* 2 pi) n))
        (a (- (mod (+ (atan2 y x) (* 0.5 p)) p) (* 0.5 p)))
        (r (len2 x y)))
    (warp f (* r (cos a)) (* r (sin a)) z)))

; Twist about z by deg-per-mm degrees per mm of height. Not
; distance-preserving. Intersect it with a containing shape such as a
; cylinder-z: interval arithmetic cannot bound a twist far out (`ink spec
; --section "Claims/Where analysis looks"`).
(fn twist-z (f deg-per-mm)
  (let ((a (deg->rad (* deg-per-mm z))))
    (warp f (+ (* (cos a) x) (* (sin a) y)) (- (* (cos a) y) (* (sin a) x)) z)))

; Extrude a 2D field g(x, y) along z between z0 and z1 (exact if g is).
; A cylinder between two heights is (extrude-z (circle-2d r) z0 z1), a ring
; (extrude-z (ring-2d ri ro) z0 z1), a box (extrude-z (rect-2d sx sy) z0 z1).
(fn extrude-z (g z0 z1)
  (emax g (- (abs (- z (* 0.5 (+ z0 z1)))) (* 0.5 (- z1 z0)))))
; The same along x (g's x along world y, its y along z) and along y (g's x
; along z, its y along x).
(fn extrude-x (g x0 x1) (warp (extrude-z g x0 x1) y z x))
(fn extrude-y (g y0 y1) (warp (extrude-z g y0 y1) z x y))

; (revolve-z g), built in, turns a 2D profile g(x, y) about the z axis:
; profile x is the radius, y the height. A profile may start at x = 0: its
; edge on the axis is no surface, and the axis is inside where it touches.

; ---------------------------------------------------------------- basic fields (primitives)

(fn sphere (r) (lipschitz 1 (- (len3 x y z) r)))                ; exact

; Box of full sizes sx sy sz, centred on the origin. Exact.
(fn box (sx sy sz)
  (let ((qx (- (abs x) (* 0.5 sx)))
        (qy (- (abs y) (* 0.5 sy)))
        (qz (- (abs z) (* 0.5 sz))))
    (lipschitz 1 (emax qx qy qz))))

(fn rounded-box (sx sy sz r)
  (- (box (- sx (* 2 r)) (- sy (* 2 r)) (- sz (* 2 r))) r))

; Cylinder along z, radius r, full height h, centred. Exact.
(fn cylinder-z (r h)
  (let ((dr (- (len2 x y) r)) (dz (- (abs z) (* 0.5 h))))
    (lipschitz 1 (emax dr dz))))
(fn cylinder-x (r h) (warp (cylinder-z r h) y z x))
(fn cylinder-y (r h) (warp (cylinder-z r h) z x y))

; Torus about z. Exact for major >= minor (ring and horn tori). For a
; spindle torus (minor > major), inside near its axis the value only
; bounds the distance (1-Lipschitz: slope at most 1).
(fn torus-z (major minor) (lipschitz 1 (- (len2 (- (len2 x y) major) z) minor)))

; Cone or frustum along z: radius r0 at z0, r1 at z1. z0 < z1 and z1 - z0
; finite, or ink check warns; a radius may be 0 (a point). Exact, with
; slope at most 1 even with its constants rounded to floats, however long
; and thin the cone.
(fn cone-z (r0 r1 z0 z1)
  (let ((dr (- r1 r0)) (dz (- z1 z0)) (l (norm2-up dr dz))
        (side (vec (/ dr l) (/ dz l)))
        (r (len2 x y)) (p (vec r z)))
    (check (value (- z1 z0) > 0) (check (value (- z1 z0) < inf)
      (lipschitz 1 (max (slab-z z0 z1)
                        (/ (* dz (- r (mix r0 r1 (clamp (/ (- z z0) dz) 0 1)))) l)
                        (corner-distance-at p (vec r0 z0) (vec -1 0) side)
                        (corner-distance-at p (vec r1 z1) (vec -1 0) (neg side))))))))

; Capsule from (0 0 0) to (0 0 h) with radius r (h may be negative). Exact.
(fn capsule-z (h r) (lipschitz 1 (- (len3 x y (- z (clamp z (min 0 h) (max 0 h)))) r)))

; Half-space n . p < d. n need not be unit; the field is exact when it is.
; No slope promise: interval arithmetic already gives a linear field's
; slope exactly.
(fn half-space (nx ny nz d) (- (+ (* nx x) (* ny y) (* nz z)) d))
(fn slab-x (x0 x1) (max (- x0 x) (- x x1)))
(fn slab-y (y0 y1) (max (- y0 y) (- y y1)))
(fn slab-z (z0 z1) (max (- z0 z) (- z z1)))
; A box by its corners, exact: (emax (slab-x x0 x1) (slab-y y0 y1) (slab-z z0 z1)).


; 2D fields in (x, y), for extrude-z and revolve-z.
(fn circle-2d (r) (lipschitz 1 (- (len2 x y) r)))
(fn rect-2d (sx sy)
  (let ((qx (- (abs x) (* 0.5 sx))) (qy (- (abs y) (* 0.5 sy))))
    (lipschitz 1 (emax qx qy))))
(fn rounded-rect-2d (sx sy r) (- (rect-2d (- sx (* 2 r)) (- sy (* 2 r))) r))
; Ring (annulus) between radii ri and ro (0 < ri < ro), centred. Exact.
(fn ring-2d (ri ro) (lipschitz 1 (max (- (len2 x y) ro) (- ri (len2 x y)))))
; Regular hexagon s across flats (a nut, a hex shaft or key), centred, its
; flats facing +-x and its corners on the y axis. Exact.
(fn hexagon-2d (s)
  (let ((h (* 0.5 s)) (a (/ h (sqrt 3))) (k (* 0.5 (sqrt 3)))
        (p (vec (abs x) (abs y))))
    (lipschitz 1 (max (corner-distance-at p (vec h a) (vec 0 -1) (vec (neg k) 0.5))
                      (corner-distance-at p (vec 0 (* 2 a)) (vec k -0.5) (vec (neg k) -0.5))))))

; Signed distance in the plane from point p to the convex wedge at corner
; c between edge directions e1 and e2, at an angle between 0 and 180
; degrees. Exact, with slope at most 1 even with its constants rounded to
; floats. A convex profile's distance is the largest of its corners'
; (cone-z, hexagon-2d). An obtuse wedge whose cosine is within
; 2e-15 m^2 / (|e1| |e2|) of 0, m the longer edge's length, is taken as
; acute.
(fn corner-distance-at (p c e1 e2)
  (let ((v (- p c))
        (cr (cross2 e1 e2)) (cs (dot e1 e2))
        (s (/ cr (abs cr)))
        (n1 (* (neg s) (vec (neg (nth e1 1)) (nth e1 0))))
        (n2 (* s (vec (neg (nth e2 1)) (nth e2 0))))
        (g1 (neg (dot v e1))) (g2 (neg (dot v e2)))
        (h1 (dot v n1)) (h2 (dot v n2))
        (m (max (norm2-up (nth e1 0) (nth e1 1)) (norm2-up (nth e2 0) (nth e2 1))))
        (err (* 2e-15 m m))
        (k (max 1 (* (/ (+ (abs cr) err) (- (neg cs) err)) (+ 1 2e-15))))
        (acute (max (emax g1 h1) (emax g2 h2)))
        (obtuse (max h1 h2 (min (emax (min g1 (* k h1)) h1) (emax (min g2 (* k h2)) h2)))))
    (/ (select (+ cs err) obtuse acute) m)))

; |(a, b)| of two numbers, rounded up, at every scale: a field divided by
; it has slope at most 1 even with its constants rounded to floats. 0
; gives the smallest float. An infinite a or b, or a length past the
; largest float, gives inf.
(fn norm2-up (a b)
  (let ((m (max (abs a) (abs b))) (n (max m 5e-324))
        (big (select (- m 1e-290) 1.0715086071862673e301 1)))
    (select (- m 1.7976931348623157e308)
      (+ (/ (* (* m big) (len2 (/ a n) (/ b n)) (+ 1 2e-15)) big) 5e-324) inf)))

;; ---------------------------------------------------------------- vectors
;
; A vector is a (vec ...) of fields. Arithmetic is element-wise, so
; (- p a), (* v 2) and (+ a b) work on points.

(def here (vec x y z))          ; the evaluation point
(def here2 (vec x y))           ; the same in the xy plane, for 2D profiles
(fn dot (a b) (reduce + (* a b)))
(fn norm (a) (sqrt (dot a a)))
(fn cross (a b)
  (vec (- (* (nth a 1) (nth b 2)) (* (nth a 2) (nth b 1)))
       (- (* (nth a 2) (nth b 0)) (* (nth a 0) (nth b 2)))
       (- (* (nth a 0) (nth b 1)) (* (nth a 1) (nth b 0)))))
(fn cross2 (a b) (- (* (nth a 0) (nth b 1)) (* (nth a 1) (nth b 0))))
(fn distance (a b) (norm (- a b)))

;; ---------------------------------------------------------------- curves
;
; A curve is a list of points: (vec u v) for profiles, (vec x y z) for paths.

; Exact distance from point p to segment a-b.
(fn segment-distance-at (p a b)
  (let ((pa (- p a)) (ba (- b a)))
    (norm (- pa (* ba (clamp (/ (dot pa ba) (dot ba ba)) 0 1))))))

; Exact distance to an open polyline, in 3D or in the xy plane (-2d).
(fn polyline-distance-at (p pts) (reduce min (map (lambda (a b) (segment-distance-at p a b)) (slice pts 0 -1) (slice pts 1))))
(fn polyline-distance (pts) (polyline-distance-at here pts))
(fn polyline-distance-2d (pts) (polyline-distance-at here2 pts))

; Close a point list into a loop (repeat the first point at the end).
(fn closed (pts) (concat pts (slice pts 0 1)))

; (polygon pts), built in, is a closed xy profile of lines and arcs
; (ink spec --section curves): exact along its edges and a bound past its
; corners, so offsets mitre them.

; A polygon of corners pts (vec u v), exact distance outside, so an outward
; offset rounds convex corners: (- (polygon-distance pts) r). Its sign, and
; its value inside (a bound), are polygon's. Each edge costs a segment
; distance, not a half-plane. A corner with a bulge or rounding (a third or
; fourth number not 0) is refused.
; It makes no slope promise: where inside meets outside the field may jump
; by rounding error. Proofs bound its slope by interval arithmetic, and cells
; straddling a jump get no bound.
(fn polygon-distance (pts)
  (let ((p (polygon pts))
        (q (map (lambda (a) (vec (nth a 0) (nth a 1))) pts))
        (d (reduce min (map (lambda (a b) (segment-distance-at here2 a b)) q (concat (slice q 1) (slice q 0 1)))))
        (curved (reduce max (map (lambda (a) (max (if (> (len a) 2) (abs (nth a 2)) 0) (if (> (len a) 3) (abs (nth a 3)) 0))) pts))))
    (if (> curved 0)
      (error "polygon-distance takes corners (vec u v), its edges straight: a corner with a bulge or a rounding (a third or fourth number not 0) is not its profile. Round convex corners by an outward offset, (- (polygon-distance pts) r), or use (polygon pts)")
      (select p p d))))

(fn regular-polygon-points (n r) (map (lambda (i) (vec (* r (cos (/ (* 2 pi i) n))) (* r (sin (/ (* 2 pi i) n))))) (range n)))

; Cubic Bezier: the point at parameter s, and n segments' points along it.
(fn bezier3 (p0 p1 p2 p3 s)
  (let ((u (- 1 s)))
    (+ (* p0 (* u u u)) (* p1 (* 3 u u s)) (* p2 (* 3 u s s)) (* p3 (* s s s)))))
(fn bezier3-points (p0 p1 p2 p3 n) (map (lambda (s) (bezier3 p0 p1 p2 p3 s)) (linspace 0 1 (+ n 1))))

; A smooth curve through all the points, 2D or 3D (Catmull-Rom), with
; per-seg segments between consecutive points.
(fn catmull-rom (p0 p1 p2 p3 s)
  (* 0.5 (+ (* 2 p1) (* (- p2 p0) s)
            (* (+ (* 2 p0) (* -5 p1) (* 4 p2) (neg p3)) (* s s))
            (* (+ (neg p0) (* 3 p1) (* -3 p2) p3) (* s s s)))))
(fn spline-points (pts per-seg)
  (let ((last (- (len pts) 1)))
    (concat
      (flatten
        (map (lambda (i)
               (let ((p0 (nth pts (max 0 (- i 1)))) (p1 (nth pts i))
                     (p2 (nth pts (+ i 1))) (p3 (nth pts (min last (+ i 2)))))
                 (map (lambda (s) (catmull-rom p0 p1 p2 p3 s)) (slice (linspace 0 1 (+ per-seg 1)) 0 -1))))
             (range last)))
      (slice pts -1))))

; Helix points about z: radius r, pitch per turn, turns, n segments.
(fn helix-points (r pitch turns n)
  (map (lambda (s) (let ((a (* 2 pi turns s))) (vec (* r (cos a)) (* r (sin a)) (* pitch turns s))))
       (linspace 0 1 (+ n 1))))

;; ---------------------------------------------------------------- sweeps and lofts

;; Every sweep checks that its path is at least as thick as its section
;; reaches (check, path-thickness). A path turning tighter, or running
;; through itself, makes a self-overlapping body; ink check says where.
;; Build path-following bodies on these, or state the same check.

; A solid rod of radius r along a path, round ends. Exact outside; inside
; a bend, where segments' rods overlap, a 1-Lipschitz bound (slope at most
; 1, smaller than the distance in magnitude).
(fn tube (pts r) (check (path-thickness pts >= r) (- (polyline-distance pts) r)))

; A hollow pipe along a path: outer radius, wall thickness. Ends are cut
; flat, square to the path; a loop has none (by-nearest-segment). Exact on
; a straight path; at a bend's round joint, a 1-Lipschitz bound.
(fn pipe (pts outer-r wall)
  (let ((last (- (len pts) 2)) (mid (- outer-r (* 0.5 wall))) (closes (path-closes pts)))
    (check (path-thickness pts >= outer-r) (by-nearest-segment pts
      (lambda (i a b)
        (let ((pa (- here a)) (ba (- b a)) (tn (unit ba))
              (u (/ (dot pa ba) (dot ba ba)))
              (u (select (if (== i 0) closes -1) (max u 0) u))
              (u (select (if (== i last) closes -1) (min u 1) u))
              (wall-d (shell (- (norm (- pa (* ba u))) mid) wall))
              (cuts (concat (if (== i 0) (vec (neg (dot pa tn))) (vec)) (if (== i last) (vec (dot (- here b) tn)) (vec)))))
          (if (== (len cuts) 0) wall-d (select closes wall-d (emax wall-d (reduce max cuts))))))))))

; Frames along a 3D polyline that twist as little as possible
; (rotation-minimising), one per segment at its start: z along it, x the
; normal (`up` picks the first), y the binormal. Carried round a loop
; (path-closes) they come back turned about the tangent; that turn is
; spread along the loop by length, so they close. They are (at ...)
; frames: (place f F) puts f there, (get F "z") is the tangent.
(fn path-frames (pts up)
  (map (lambda (r) (axes (origin (nth r 0)) (x (nth r 2)) (z (nth r 1)))) (path-rmf pts up)))
; path-frames as vectors: (vec origin tangent normal binormal turn) per
; segment, turn the radians its frame turns on along it (a loop's share of
; its closing turn; 0 on an open path).
(fn path-rmf (pts up)
  (let ((tans (map (lambda (a b) (unit (- b a))) (slice pts 0 -1) (slice pts 1)))
        (n0 (unit (cross (cross (nth tans 0) up) (nth tans 0))))
        (rmf (fold (lambda (acc t1)
                     (let ((prev (nth acc -1)))
                       (concat acc (vec (vec t1 (rmf-step (nth prev 0) (nth prev 1) t1))))))
                   (vec (vec (nth tans 0) n0))
                   (slice tans 1)))
        (lens (map (lambda (a b) (norm (- b a))) (slice pts 0 -1) (slice pts 1)))
        (total (reduce + lens))
        (starts (fold (lambda (acc l) (concat acc (vec (+ (nth acc -1) l)))) (vec 0) (slice lens 0 -1)))
        ; The last frame carried on to the first tangent, and the turn
        ; about it to the first frame (below 1e-12, rounding: none).
        (back (rmf-step (nth (nth rmf -1) 0) (nth (nth rmf -1) 1) (nth tans 0)))
        (turn (select (path-closes pts) (atan2 (dot (cross back n0) (nth tans 0)) (dot back n0)) 0))
        (turn (select (- (abs turn) 1e-12) 0 turn)))
    (map (lambda (o f s l)
           (let ((a (* turn (/ s total))) (tn (nth f 0)) (n (nth f 1)) (b (cross (nth f 0) (nth f 1))))
             (vec o tn (+ (* n (cos a)) (* b (sin a))) (- (* b (cos a)) (* n (sin a))) (* turn (/ l total)))))
         (slice pts 0 -1) rmf starts lens)))

; Where `here` lies in the section of segment i (a to b) of a path with
; path-rmf's `frames`: its frame's x and y, turned on by the frame's turn as
; here goes from the mitre at a to the mitre at b (the planes bisecting the
; joints), so neighbouring sections meet on their mitre.
(fn path-section-xy (frames i a b)
  (let ((f (nth frames i)) (q (- here a))
        (ta (+ (nth (nth frames (- i 1)) 1) (nth f 1)))
        (tb (+ (nth f 1) (nth (nth frames (if (== i (- (len frames) 1)) 0 (+ i 1))) 1)))
        (da (/ (dot q ta) (max (norm ta) 1e-9)))
        (db (/ (dot (- b here) tb) (max (norm tb) 1e-9)))
        (th (* (nth f 4) (clamp (/ da (max (+ da db) (* 1e-9 (norm (- b a))))) 0 1)))
        (u (dot q (nth f 2))) (v (dot q (nth f 3))))
    (vec (+ (* u (cos th)) (* v (sin th))) (- (* v (cos th)) (* u (sin th))))))

; Parallel transport: normal r0 turned by the rotation taking tangent t0 to
; t1. Parallel tangents leave r0 as it is.
(fn rmf-step (t0 r0 t1)
  (let ((c (cross t0 t1)) (cs (dot t0 t1)))
    (+ (* r0 cs) (cross c r0) (* c (/ (dot c r0) (+ 1 cs))))))

; Negative on a loop: a path whose first and last points are within 1e-9
; mm. A loop has no ends; its last segment joins its first as any two do.
; A number, for select, so a params-shaped path is not fixed at load time.
(fn path-closes (pts) (- (norm (- (nth pts -1) (nth pts 0))) 1e-9))

; The value of the path's nearest segment: (per-seg i a b) is the field
; where segment i, a to b, is nearest. An open path's end segments are cut
; flat; on a loop the last segment precedes the first (path-closes).
; Neighbouring segments meet at the plane bisecting their joint (a mitre),
; outside the joint too. Interval bounds skip segments far from the point,
; so evaluation stays local.
(fn by-nearest-segment (pts per-seg)
  (let ((last (- (len pts) 2))
        (closes (path-closes pts))
        (tans (map (lambda (a b) (unit (- b a))) (slice pts 0 -1) (slice pts 1)))
        (pieces (map (lambda (i a b)
                       (let ((tn (nth tans i))
                             (v (per-seg i a b))
                             (before (max 0 (neg (dot (- here a) (+ (nth tans (- i 1)) tn)))))
                             (after (max 0 (dot (- here b) (+ tn (nth tans (if (== i last) 0 (+ i 1)))))))
                             ; An open path's ends are cut flat, not mitred.
                             (v (if (== i 0) (select closes v (max v (neg (dot (- here a) tn)))) v))
                             (v (if (== i last) (select closes v (max v (dot (- here b) tn))) v))
                             (before (if (== i 0) (select closes before 0) before))
                             (after (if (== i last) (select closes after 0) after)))
                         (vec (+ (segment-distance-at here a b) (+ before after)) v)))
                     (range (+ last 1)) (slice pts 0 -1) (slice pts 1))))
    (nth (reduce (lambda (a b) (select (- (nth a 0) (nth b 0)) a b)) pieces) 1)))

; Sweep a 2D profile (a field in x, y) along a 3D polyline: its x follows
; the rotation-minimising frame's normal, its y the binormal (path-frames:
; on a loop they turn along each segment by its share of the closing turn,
; from mitre to mitre). Joints are mitred, a sharp mitre's outer corner
; rounded at the profile's reach (as in tube and pipe). Ends are flat; a
; loop has none. Along a straight run the field is the profile's (exact if
; it is and its frame does not turn); past the flat ends' rims, a bound.
; Sample curves finely.
(fn sweep (profile pts up)
  (let ((frames (path-rmf pts up))
        (reach (profile-reach profile)))
    (check (path-thickness pts >= reach) (by-nearest-segment pts
      (lambda (i a b)
        (let ((xy (path-section-xy frames i a b)))
          (max (warp profile (nth xy 0) (nth xy 1) 0)
               (- (segment-distance-at here a b) reach))))))))

; Sweep a section that varies along the path. (section s) is the 2D profile
; at s, the fraction of the path's length (0 start, 1 end): it can scale,
; round, twist or morph. Frames, mitres and ends are as in sweep.
; Any field of s works, but s must not be used at load time (in bounds or
; if). A size that reads s drops primitives' slope promises (std's header):
; proofs bound the sweep's slope by interval arithmetic.
; All sections must be proven bounded at once, which interval bounds cannot
; do for a twist (twist-z): intersect a twisting section with a containing
; shape, (max (circle-2d R) ...).
(fn sweep-varying (section pts up)
  (let ((frames (path-rmf pts up))
        (lens (map (lambda (a b) (norm (- b a))) (slice pts 0 -1) (slice pts 1)))
        (starts (fold (lambda (acc l) (concat acc (vec (+ (nth acc -1) l)))) (vec 0) (slice lens 0 -1)))
        (total (reduce + lens))
        ; The reach of every section: s = (z + 1) / 2 over z in [-1, 1].
        (reach (profile-reach (warp (section z) x y (* 0.5 (+ z 1))))))
    (check (path-thickness pts >= reach) (by-nearest-segment pts
      (lambda (i a b)
        (let ((f (nth frames i)) (q (- here a)) (xy (path-section-xy frames i a b))
              (s (clamp (/ (+ (nth starts i) (dot q (nth f 1))) total) 0 1)))
          (max (warp (section z) (nth xy 0) (nth xy 1) s)
               (- (segment-distance-at here a b) reach))))))))


; Loft from 2D profile g0 at z0 to g1 at z1 (fields in x, y), blending
; linearly. Not an exact distance.
(fn loft-z (g0 g1 z0 z1)
  (max (blend g0 g1 (clamp (/ (- z z0) (- z1 z0)) 0 1)) (slab-z z0 z1)))

;; ---------------------------------------------------------------- paths of lines and arcs (belts, chains, cables)
;
; Lines and circular arcs end to end in 3D, each tangent to the next: a path
; a belt, chain or cable moves along. Its pieces may read t (a belt round a
; moving carriage). Station s is mm along it from its start. Belts and
; chains on paths: below (belts and chains).
;   (path-line a b)              from point a to point b
;   (path-arc c u r f sweep)     about centre c and unit axis u, radius r,
;                                from c + r f (f a unit vector square to u),
;                                sweep radians, right-handed about u
;   (path-arc c u r f sweep rho) ... with rho mm of station per radian, not
;                                r: a chain on an n-tooth sprocket advances
;                                n pitch / 2 pi per radian, though its pins
;                                sit on the larger pitch circle
;                                (roller-chain-path)
;   (path pieces up)             an open path; (path-loop pieces up) a closed
;                                one, its stations wrapping round its length
;   (path-length P)              mm
;   (path-at P s)                the frame at station s: origin on the path,
;                                x along it, y = z cross x, z up (a planar
;                                path's normal) made square to the path at
;                                its start and turned with it about each
;                                arc's axis
;   (path-point P s)             its origin, as a vec
;   (path-link-at P s len)       a rigid link's frame from station s to
;                                s + len (a chain link, a tread): origin at
;                                the point at s, x along the chord to the
;                                point at s + len, z as at s + len / 2
; Past an open path's ends its first and last pieces run on. up must be
; known when the model loads. ink check reports pieces that do not join end
; to end and tangent; round a corner with an arc. The frame moves without
; jumps as s changes; a loop's wrap jumps in s, not in the frame.
; A piece is (vec kind length ...): 0 line (a b), 1 arc (c u r f sweep rho).
(fn path-line (a b) (vec 0 (norm (- b a)) a b))
(fn path-arc (c u r f sweep rest ...)
  (let ((rho (if (== (len rest) 0) r (nth rest 0))))
    (vec 1 (* rho (abs sweep)) c u r f sweep rho)))
; The angle an arc has turned u mm along it, signed as its sweep.
(fn path-arc-angle (p u) (/ (* u (nth p 6)) (nth p 1)))
; A piece's displacement from its start u mm along it, and its tangent there.
(fn path-piece-move (p u)
  (if (== (nth p 0) 0)
    (* (/ (- (nth p 3) (nth p 2)) (nth p 1)) u)
    (let ((r (nth p 4)) (f (nth p 5)) (th (path-arc-angle p u)))
      (* (+ (* f (- (cos th) 1)) (* (cross (nth p 3) f) (sin th))) r))))
(fn path-piece-tangent (p u)
  (if (== (nth p 0) 0)
    (/ (- (nth p 3) (nth p 2)) (nth p 1))
    (let ((f (nth p 5)) (th (path-arc-angle p u)) (sense (/ (nth p 6) (abs (nth p 6)))))
      (* (- (* (cross (nth p 3) f) (cos th)) (* f (sin th))) sense))))
(fn path-piece-start (p) (if (== (nth p 0) 0) (nth p 2) (+ (nth p 2) (* (nth p 5) (nth p 4)))))
(fn path-of (pieces up closed)
  (let ((n (len pieces))
        (lens (map (lambda (p) (nth p 1)) pieces))
        (next (lambda (i) (nth pieces (if (== i (- n 1)) 0 (+ i 1)))))
        ; How far each piece's end is from the next's start, in place and
        ; direction. A loop's last piece joins its first.
        (gaps (map (lambda (i)
                     (let ((p (nth pieces i)) (q (next i)))
                       (max (norm (- (+ (path-piece-start p) (path-piece-move p (nth p 1))) (path-piece-start q)))
                            (norm (- (path-piece-tangent p (nth p 1)) (path-piece-tangent q 0))))))
                   (range (if (== closed 1) n (- n 1))))))
    (check (value (if (== (len gaps) 0) 0 (reduce max gaps)) < 1e-6)
      (vec pieces (unit up) closed
           (fold (lambda (acc l) (concat acc (vec (+ (nth acc -1) l)))) (vec 0) (slice lens 0 -1))
           (reduce + lens)))))
(fn path (pieces up) (path-of pieces up 0))
(fn path-loop (pieces up) (path-of pieces up 1))
(fn path-length (P) (nth P 4))
; How far station s is along piece i, clamped to the piece, except where
; an open path's first and last pieces run on.
(fn path-along (P i s)
  (let ((pieces (nth P 0)) (u (- s (nth (nth P 3) i))) (open (== (nth P 2) 0)))
    (if (* open (== i 0))
      (if (== (len pieces) 1) u (min u (nth (nth pieces i) 1)))
      (if (* open (== i (- (len pieces) 1))) (max u 0) (clamp u 0 (nth (nth pieces i) 1))))))
(fn path-station (P s) (if (== (nth P 2) 1) (mod s (nth P 4)) s))
(fn path-point (P s)
  (let ((pieces (nth P 0)) (s (path-station P s)))
    (fold (lambda (acc i) (+ acc (path-piece-move (nth pieces i) (path-along P i s))))
          (path-piece-start (nth pieces 0)) (range (len pieces)))))
; A frame with z along up (a unit vector, known when the model loads).
(fn path-up-frame (up)
  (let ((ax (cross (vec 0 0 1) up)) (sn (norm ax)))
    (if (< sn 1e-12)
      (if (> (nth up 2) 0) (at) (at (rot-x 180)))
      (at (rot (nth ax 0) (nth ax 1) (nth ax 2) (* (/ 180 pi) (atan2 sn (nth up 2))))))))
(fn path-at (P s)
  (let ((pieces (nth P 0)) (s (path-station P s)) (up (nth P 1))
        (t0 (path-piece-tangent (nth pieces 0) 0)) (base (path-up-frame up))
        (a (dot t0 (get base "x"))) (b (dot t0 (get base "y"))) (h (dot t0 up))
        (o (path-point P s))
        ; the frame at the start, turned about each arc's axis by the part
        ; of the arc before s
        (turns (fold (lambda (acc i)
                       (let ((p (nth pieces i)))
                         (if (== (nth p 0) 0) acc
                           (let ((u (nth p 3)))
                             (concat acc (vec (at (rot (nth u 0) (nth u 1) (nth u 2) (* (/ 180 pi) (path-arc-angle p (path-along P i s)))))))))))
                     (vec) (range (len pieces)))))
    (at (fold (lambda (f g) (at f g))
              (at (rot-y (neg (* (/ 180 pi) (atan2 h (len2 a b))))) (rot-z (* (/ 180 pi) (atan2 b a))) base)
              turns)
        (move (nth o 0) (nth o 1) (nth o 2)))))
(fn path-link-at (P s len)
  (let ((m (path-at P (+ s (* 0.5 len)))) (a (path-point P s)) (ch (- (path-point P (+ s len)) a))
        (d (- a (get m "origin"))))
    (at (rot-z (* (/ 180 pi) (atan2 (dot ch (get m "y")) (dot ch (get m "x"))))) m (move (nth d 0) (nth d 1) (nth d 2)))))

; Fields along a planar path, one in a plane square to its up (a belt, a
; track). Near the path a point has path coordinates: station s, offset n
; across (along the frame's y, up cross the tangent: left of travel) and
; height w along up. Each piece gives them in its region, between the lines
; through its ends square to it; its line or circle runs on past its ends.
; A field made of them is a union of one term per piece, read from its
; piece's coordinates, numbers and start station alone. Each term is cut at
; an open path's ends (its first piece at station 0, its last at its
; length; a line there is cut flat) and bounded by its piece's box (a
; line's ends, an arc's whole circle) grown by the term's reach.
;   (path-sweep SECTION P)     SECTION (a 2D profile: x across, y up) carried
;                              along P, bent with it: a band, a hose
;   (path-repeat F P PITCH SHIFT)  F (a 2D profile in the path's plane: x
;                              along, y across, symmetric about x = 0) placed
;                              whole, upright and unbent, at every station
;                              k PITCH + SHIFT: a belt's teeth, a track's
;                              cleats. F must reach under half a pitch along
;                              the path from its station (round an arc, at
;                              its radius).
;   (path-pieces P PER REACH)  the union of (PER c i), piece i's term at its
;                              coordinates c, cut and bounded so; REACH is
;                              how far a term reaches from its piece
;   (path-sweep-at SECTION c), (path-repeat-at R P i c) with R from
;   (path-features F PITCH SHIFT): the two's terms, to make both (a belt's
;   band and teeth) reading each piece once:
;     (path-pieces P (lambda (c i) (min (path-sweep-at S c) (path-repeat-at R P i c))) reach)
; Piece i's (vec s n w out) at here, out <= 0 in its region. An open path's
; first and last lines are not clipped at its ends.
(fn path-piece-coords (P i)
  (let ((p (nth (nth P 0) i)) (up (nth P 1)) (s0 (nth (nth P 3) i))
        (open (== (nth P 2) 0)) (n (len (nth P 0))))
    (if (== (nth p 0) 0)
      (let ((a (nth p 2)) (L (nth p 1)) (tn (/ (- (nth p 3) a) L)) (q (- here a)) (along (dot q tn))
            (cuts (concat (if (* open (== i 0)) (vec) (vec (neg along))) (if (* open (== i (- n 1))) (vec) (vec (- along L))))))
        (vec (+ s0 along) (dot q (cross up tn)) (dot q up) (if (== (len cuts) 0) -1 (reduce max cuts))))
      (let ((q (- here (nth p 2))) (u (nth p 3)) (sweep (nth p 6)) (sense (/ sweep (abs sweep)))
            (h (* 0.5 (abs sweep)))
            (mid (+ (* (nth p 5) (cos (* 0.5 sweep))) (* (cross u (nth p 5)) (sin (* 0.5 sweep)))))
            (side (cross u mid))
            ; the angle about u from the arc's middle; atan2's cut is opposite
            (a (atan2 (dot q side) (dot q mid)))
            ; past the radial lines through its ends: a wedge, reflex past pi
            (h1 (+ (* (neg (sin h)) (dot q mid)) (* (cos h) (dot q side))))
            (h2 (- (* (neg (sin h)) (dot q mid)) (* (cos h) (dot q side)))))
        (vec (+ s0 (* (nth p 7) (+ h (* sense a))))
             (* sense (dot u up) (- (nth p 4) (len2 (dot q mid) (dot q side))))
             (dot q up)
             (select (- h (* 0.5 pi)) (max h1 h2) (min h1 h2)))))))
(fn path-pieces (P per reach)
  (let ((pieces (nth P 0)) (n (len pieces)) (open (== (nth P 2) 0)))
    (reduce min (map (lambda (i)
      (let ((p (nth pieces i)) (c (path-piece-coords P i))
            (ends (if (== (nth p 0) 0) (vec (nth p 2) (nth p 3)) (vec (- (nth p 2) (nth p 4)) (+ (nth p 2) (nth p 4)))))
            (lo (- (map min (nth ends 0) (nth ends 1)) reach)) (hi (+ (map max (nth ends 0) (nth ends 1)) reach)))
        (max (per c i) (reduce max (concat (- lo here) (- here hi)))
             (if (* open (== i 0)) (neg (nth c 0)) -1e9)
             (if (* open (== i (- n 1))) (- (nth c 0) (nth (nth P 3) i) (nth p 1)) -1e9))))
      (range n)))))
; A section's term, clipped to its piece's region grown 0.2 mm so no crack
; shows where pieces join (a line parts from its tangent arc there by
; 0.2^2 / 2 r: microns).
(fn path-sweep-at (section c) (max (warp section (nth c 1) (nth c 2) z) (- (nth c 3) 0.2)))
(fn path-sweep (section P) (path-pieces P (lambda (c i) (path-sweep-at section c)) (profile-reach section)))
; Each piece places the features whose stations lie on it, turned with the
; path, unbent, and owns them out to the midpoints to their neighbours (an
; open path's end pieces to its ends, where features are cut). So a feature
; overhanging its piece's end (a belt tooth across a tangent point) stays
; whole.
(fn path-features (feature pitch shift)
  (let ((box (bounds (max feature (- (abs z) 1)))))
    (vec feature pitch shift (nth (nth box 0) 1) (nth (nth box 1) 1))))
(fn path-repeat-at (R P i c)
  (let ((feature (nth R 0)) (pitch (nth R 1)) (shift (nth R 2)) (lo (nth R 3)) (hi (nth R 4))
        (pieces (nth P 0)) (n (len pieces)) (open (== (nth P 2) 0)) (up (nth P 1))
        (p (nth pieces i)) (s (nth c 0)) (s0 (nth (nth P 3) i))
        ; the midpoint before the first feature at or past station s
        (gap (lambda (s) (- (+ shift (* pitch (neg (floor (/ (- shift s) pitch))))) (* 0.5 pitch))))
        (tri (+ 0.25 (/ (asin (neg (cos (* 2 pi (/ (- s shift) pitch))))) (* 2 pi))))   ; 0 at a feature, 1/2 between
        (local (if (== (nth p 0) 0)
                 (vec (* pitch tri) (nth c 1))
                 (let ((d (* (/ (nth p 6) (abs (nth p 6))) (dot (nth p 3) up))) (r (nth p 4))
                       (rho (- r (* d (nth c 1)))) (al (* (/ pitch (nth p 7)) tri)))
                   (vec (* rho (sin al)) (* d (- r (* rho (cos al)))))))))
    (max (warp feature (nth local 0) (nth local 1) z) (- lo (nth local 1)) (- (nth local 1) hi)
         ; owned features, cut only within a pitch of the piece's ends
         (min (max (if (* open (== i 0)) -1 (- (gap s0) s))
                   (if (* open (== i (- n 1))) -1 (- s (gap (+ s0 (nth p 1))))))
              (+ (nth c 3) pitch))
         (- (nth c 3) (* 0.5 pitch)))))
(fn path-repeat (feature P pitch shift)
  (let ((R (path-features feature pitch shift)))
    (path-pieces P (lambda (c i) (path-repeat-at R P i c)) (profile-reach feature))))

; A body that is a profile carried along a path (a belt, a cable, a hose):
; (attr carried (carried PROFILE P LENGTH SHIFT)) says its shape is
; PROFILE, a solid in material coordinates (x the station along the
; material, 0 to LENGTH; y and z across, along the frame's y and z), with
; material station m at path station m + SHIFT. SHIFT may read t (a loop
; circulating). The shape stays the body's field; this lets a viewer bend
; one mesh of the profile along the path at every moment (ink mesh --out
; FILE.glb reports how far the two disagree).
(fn carried (profile P length shift) (vec profile P length shift))

;; ---------------------------------------------------------------- belts and chains

; A belt profile b is a record, like the catalog's standards.json
; "belt-profiles": pitch, pitch-line-differential (land to pitch line),
; back (land to back), tooth-height, and its tooth, as "tooth" (its
; drawing's segments, for timing-tooth-2d) or, where those cannot draw it,
; "tooth-polygon". Paths built here lie in z = 0 with up +z; wheels'
; circles are in x and y.
;   (belt-along b P width side shift)      b along any planar pitch path P
;   (belt-along-pitched b p P width side shift)   ... its teeth p apart
;   (belt-carried b P width side shift)    the same as std's carried profile
; A closed belt on two pulleys, A (ta teeth) at the origin, B (tb) at (c, 0),
; p its pitch:
;   (sync-pitch-r p n) (sync-belt-length p ta tb c)
;   (sync-belt-center-distance p ta tb belt-teeth)  centres that close it
;   (sync-belt-phase-a p ta tb c) (sync-belt-phase-b p ta tb c)  wheel turns, deg
;   (sync-belt-path p ta tb c) (sync-belt b ta tb c width shift) (sync-belt-carried ...)
; A closed belt round any wheels, inside the loop or on its back:
;   (belt-loop-path circles) (sync-belt-around b circles teeth width shift)
; An open belt from start to end round circles (vec cx cy r dir face):
;   (belt-route b start circles end width shift)   its shape
;   (belt-route-path b start circles end)          its pitch path
;   (belt-route-length b start circles end)        mm
;   (belt-route-phase b start circles end i shift) circle i's turn, deg
;   (belt-route-span b start circles end k z)      its k-th run at height z
;   (belt-route-carried b start circles end width shift)
; A roller chain of pitch p on sprockets of ta and tb teeth, c apart:
;   (roller-chain-path p ta tb c)
;   (roller-chain-phase-a p ta tb c shift) (roller-chain-phase-b p ta tb c shift)

;; ---------------------------------------------------------------- a belt along a path

; One tooth of belt profile b in its plane: x along the belt from the
; tooth's middle, y from the land toward its tip. "tooth" is its right half
; from the tip (0, tooth-height) down to the land, as its drawing's segments
; (u v r cu cv). Each ends at (u, v): an arc of radius |r| about (cu, cv),
; the tooth inside (r > 0) or outside (r < 0, a root round), or a line
; (r = 0), each bounding the tooth over the heights it spans. The tooth is
; exactly as drawn: discs, their outsides and half-planes. Below the land
; it runs back / 2 into the band, as wide as where its last segment meets
; the land. A tooth the segments cannot draw (a concave tip: not one
; boundary from the axis down) is "tooth-polygon", a closed std polygon in
; the same coordinates, with an empty "tooth". It runs as far into the
; band, within its pitch.
(fn belt-tooth (b)
  (if (> (len (get b "tooth")) 0) (belt-tooth-segments b)
    (max (polygon (get b "tooth-polygon")) (- (* -0.5 (get b "back")) y)
         (- (abs x) (* 0.5 (get b "pitch"))))))
(fn belt-tooth-segments (b)
  (timing-tooth-2d (get b "tooth") (get b "tooth-height") (get b "back")))

; A synchronous belt of any tooth profile b along a planar pitch path P,
; width wide about the path's plane. Its teeth face the frame's y if side
; is 1 (the inside of a loop running counter-clockwise about up), -y if -1.
; The band is carried along the path, bent with it (path-sweep-at). Each
; tooth is placed whole and unbent, upright on the pitch line at station
; (k + 1/2) pitch + shift and turned with the path round a wrap
; (path-repeat-at). Lands are at k pitch + shift. On a pulley turned to put
; a groove under each tooth (sync-belt-phase-a and -b, belt-route-phase) it
; sits with the groove form's own backlash, its land on the tip circle.
(fn belt-along (b P width side shift) (belt-along-pitched b (get b "pitch") P width side shift))
; The same with teeth p apart along P, each as drawn (a belt stretched
; round its loop: p is its length over its tooth count).
(fn belt-along-pitched (b p P width side shift)
  (let ((pld (get b "pitch-line-differential"))
        (tooth (warp (belt-tooth b) x (- (* side y) pld) z))
        (band (max (- (* side x) pld) (- (- pld (get b "back")) (* side x)) (- (abs y) (* 0.5 width))))
        (teeth (path-features tooth p (+ shift (* 0.5 p)))))
    ; each piece's band and teeth from its coordinates c, read once
    (path-pieces P (lambda (c i) (min (tag belt-back (path-sweep-at band c))
                                      (tag belt-teeth (max (path-repeat-at teeth P i c) (- (abs (nth c 2)) (* 0.5 width))))))
                 (max (profile-reach band) (profile-reach tooth)))))
; The same belt as a std carried profile, for (attr carried ...) beside its
; shape: laid straight along x, teeth as on a run, unbounded (the carried
; length cuts it). shift circulates a loop's material and phases an open
; belt's teeth. A viewer bending it bends the teeth round a wrap, unlike
; the belt; the GLB export reports by how much.
(fn belt-carried (b P width side shift)
  (let ((p (get b "pitch")) (pld (get b "pitch-line-differential")) (n (* side y))
        (phase (if (== (nth P 2) 1) 0 shift))
        (straight (max (min (max (- n pld) (- (- pld (get b "back")) n))
                            (warp (belt-tooth b) (- (mod (- x phase) p) (* 0.5 p)) (- n pld) z))
                       (- (abs z) (* 0.5 width)))))
    (carried straight P (path-length P) (- shift phase))))

;; ---------------------------------------------------------------- a closed belt on two pulleys

; A closed belt on two pulleys, p its pitch: A (ta teeth) centred at the
; origin, B (tb) at (c, 0). Its pitch path runs counter-clockwise from A's
; upper tangent point, a land there at shift 0: round A, the lower run to
; B, round B and back. shift circulates the belt, mm along its pitch path
; (a pulley's turn times its pitch radius).
(fn sync-pitch-r (p n) (* 0.5 (/ (* p n) pi)))
(fn sync-belt-theta (p ta tb c)                  ; angle of A's upper tangent point, radians
  (acos (/ (- (sync-pitch-r p ta) (sync-pitch-r p tb)) c)))
(fn sync-belt-length (p ta tb c)
  (let ((ra (sync-pitch-r p ta)) (rb (sync-pitch-r p tb)) (th (sync-belt-theta p ta tb c)))
    (+ (* ra (- (* 2 pi) (* 2 th))) (* rb 2 th) (* 2 c (sin th)))))
; The c that closes a belt of belt-teeth teeth: Newton's method on
; L(c) = p belt-teeth from c = (L - pi (ra + rb)) / 2.
(fn sync-belt-center-distance (p ta tb belt-teeth)
  (let ((L (* p belt-teeth))
        (c0 (* 0.5 (- L (* pi (+ (sync-pitch-r p ta) (sync-pitch-r p tb)))))))
    (fold (lambda (c i) (- c (/ (- (sync-belt-length p ta tb c) L) (* 2 (sin (sync-belt-theta p ta tb c))))))
          c0 (range 8))))
; Pulley turns (degrees) that put a groove under each tooth at shift 0.
(fn sync-belt-phase-a (p ta tb c) (* (/ 180 pi) (sync-belt-theta p ta tb c)))
(fn sync-belt-phase-b (p ta tb c)
  (let ((ra (sync-pitch-r p ta)) (rb (sync-pitch-r p tb)) (th (sync-belt-theta p ta tb c))
        (s-b0 (+ (* ra (- (* 2 pi) (* 2 th))) (* c (sin th)))))   ; station at B's lower tangent point
    (* (/ 180 pi) (- (neg th) (/ s-b0 rb)))))
; Its pitch path: belt-loop-path round A and B.
(fn sync-belt-path (p ta tb c)
  (belt-loop-path (vec (vec 0 0 (sync-pitch-r p ta) 1) (vec c 0 (sync-pitch-r p tb) 1))))
; ink check reports a pulley too small to wrap, one whose pitch radius is
; under the belt's reach from its pitch line (pitch-line differential plus
; tooth height).
(fn sync-belt (b ta tb c width shift)
  (let ((p (get b "pitch")))
    (check (value (min (sync-pitch-r p ta) (sync-pitch-r p tb)) >= (+ (get b "pitch-line-differential") (get b "tooth-height")))
      (belt-along b (sync-belt-path p ta tb c) width 1 shift))))
(fn sync-belt-carried (b ta tb c width shift)
  (belt-carried b (sync-belt-path (get b "pitch") ta tb c) width 1 shift))

;; ---------------------------------------------------------------- a closed belt round any wheels

; A closed belt's pitch path round any number of wheels (z = 0, up +z),
; starting round the first circle. circles are (vec cx cy r dir), in the order a
; counter-clockwise belt meets them, r the pitch path's radius there.
; dir 1 is a pulley inside the loop, r its teeth times the pitch over 2 pi.
; dir -1 is an idler outside, on the belt's back, r its radius plus the
; back less the pitch-line differential.
(fn belt-loop-path (circles)
  (let ((n (len circles))
        (segs (map (lambda (i) (belt-route-tangent (nth circles i) (nth circles (mod (+ i 1) n)))) (range n)))
        (arc (lambda (i)
               (let ((node (nth circles i)) (prev (nth segs (mod (+ i (- n 1)) n))) (next (nth segs i)) (r (nth node 2)))
                 (path-arc (vec (nth node 0) (nth node 1) 0) (vec 0 0 (nth node 3)) r
                           (vec (/ (- (nth prev 2) (nth node 0)) r) (/ (- (nth prev 3) (nth node 1)) r) 0)
                           (belt-route-sweep node prev next)))))
        (line (lambda (sg) (path-line (vec (nth sg 0) (nth sg 1) 0) (vec (nth sg 2) (nth sg 3) 0)))))
    (path-loop (fold (lambda (acc i) (concat acc (vec (arc i) (line (nth segs i))))) (vec) (range n))
               (vec 0 0 1))))
; Belt b of teeth teeth round circles (belt-loop-path's), its teeth evenly
; spaced, the loop's length over teeth apart: stretched to its loop where
; that is not b's pitch, as an idler-tensioned belt is. ink check reports a
; wheel too small to wrap: a pulley under the teeth's reach, an idler under
; the back's.
(fn sync-belt-around (b circles teeth width shift)
  (let ((pld (get b "pitch-line-differential"))
        (room (fold (lambda (m c) (min m (- (nth c 2) (if (> (nth c 3) 0) (+ pld (get b "tooth-height")) (- (get b "back") pld)))))
                    1e9 circles)))
    (check (value room >= 0)
      (let ((P (belt-loop-path circles))) (belt-along-pitched b (/ (path-length P) teeth) P width 1 shift)))))

;; ---------------------------------------------------------------- an open belt

; An open belt from start to end, (vec x y) each, round circles (vec cx cy
; radius direction face) in path order. direction is 1 where the belt runs
; round it counter-clockwise, -1 clockwise. face 1: its teeth on the circle,
; radius the pitch radius. face 0: its back on a smooth drum, radius the
; drum's (the pitch path adds the back less the pitch-line differential).
; Adjacent circles need a nondegenerate common tangent, and direction *
; (2 face - 1) must be the same on every circle: a belt in one plane keeps
; its teeth on one edge. Ends and centres may read t, but direction and
; face must be fixed when the model loads.
(fn belt-route-node (b c)
  (vec (nth c 0) (nth c 1)
       (+ (nth c 2) (if (== (nth c 4) 0) (- (get b "back") (get b "pitch-line-differential")) 0))
       (nth c 3) (nth c 4) (nth c 2)))
(fn belt-route-end-node (p) (vec (nth p 0) (nth p 1) 0 0 0 0))
(fn belt-route-nodes (b start circles end)
  (concat (vec (belt-route-end-node start)) (map (lambda (c) (belt-route-node b c)) circles)
          (vec (belt-route-end-node end))))

; The tangent run from circle a to b, direct or crossed by the signs of
; their directions, (vec ax ay bx by tx ty L): where it leaves a and meets
; b, its unit direction and length.
(fn belt-route-tangent (a b)
  (let ((dx (- (nth b 0) (nth a 0))) (dy (- (nth b 1) (nth a 1)))
        (q (- (* (nth b 2) (nth b 3)) (* (nth a 2) (nth a 3))))
        (d2 (+ (square dx) (square dy))) (L (sqrt (- d2 (square q))))
        (tx (/ (+ (* L dx) (* q dy)) d2))
        (ty (/ (- (* L dy) (* q dx)) d2))
        (ax (+ (nth a 0) (* (nth a 2) (nth a 3) ty)))
        (ay (- (nth a 1) (* (nth a 2) (nth a 3) tx)))
        (bx (+ (nth b 0) (* (nth b 2) (nth b 3) ty)))
        (by (- (nth b 1) (* (nth b 2) (nth b 3) tx))))
    (vec ax ay bx by tx ty L)))
(fn belt-route-segments (nodes)
  (map (lambda (i) (belt-route-tangent (nth nodes i) (nth nodes (+ i 1))))
       (range (- (len nodes) 1))))
; Run k at height z (run 0 from start to the first circle), (vec p q) for
; std's belt-span: p where it leaves its earlier end, q where it meets its
; later.
(fn belt-route-span (b start circles end k z)
  (let ((sg (nth (belt-route-segments (belt-route-nodes b start circles end)) k)))
    (vec (vec (nth sg 0) (nth sg 1) z) (vec (nth sg 2) (nth sg 3) z))))
(fn belt-route-entry-angle (node prev)
  (atan2 (- (nth prev 3) (nth node 1)) (- (nth prev 2) (nth node 0))))
; The angle a node turns from entry to exit in its direction, in (0, 2 pi),
; with tight interval bounds while the circles move with t.
(fn belt-route-sweep (node prev next)
  (let ((ex (- (nth prev 2) (nth node 0))) (ey (- (nth prev 3) (nth node 1)))
        (xx (- (nth next 0) (nth node 0))) (xy (- (nth next 1) (nth node 1))))
    (+ pi (atan2 (neg (* (nth node 3) (- (* ex xy) (* ey xx)))) (neg (+ (* ex xx) (* ey xy)))))))
(fn belt-route-arc-length (node prev next)
  (* (nth node 2) (belt-route-sweep node prev next)))
(fn belt-route-station (nodes segs i)
  (fold (lambda (acc j)
          (+ acc (nth (nth segs j) 6)
             (belt-route-arc-length (nth nodes (+ j 1)) (nth segs j) (nth segs (+ j 1)))))
        0 (range i)))
(fn belt-route-length (b start circles end)
  (let ((nodes (belt-route-nodes b start circles end))
        (segs (belt-route-segments nodes))
        (last (- (len segs) 1)))
    (+ (belt-route-station nodes segs last) (nth (nth segs last) 6))))

; The turn (degrees) of the wheel at circle i that puts a groove under each
; tooth. A positive shift turns it in the circle's direction.
(fn belt-route-phase (b start circles end i shift)
  (let ((nodes (belt-route-nodes b start circles end)) (segs (belt-route-segments nodes))
        (node (nth nodes (+ i 1))) (prev (nth segs i))
        (station (+ (belt-route-station nodes segs i) (nth prev 6))))
    (* (/ 180 pi) (+ (belt-route-entry-angle node prev)
                     (/ (* (nth node 3) (- shift station)) (nth node 2))))))

; The pitch path as points: the ends, and each wrap in k segments from
; entry to exit tangent point, within rp (1 - cos (a / 2 k)) of an arc of
; radius rp turning a (0.03 mm for k = 16 on a half turn of a 20T GT2
; pulley).
(fn belt-route-points (b start circles end k)
  (let ((nodes (belt-route-nodes b start circles end)) (segs (belt-route-segments nodes))
        (arc (lambda (i)
               (let ((node (nth nodes (+ i 1))) (prev (nth segs i)) (next (nth segs (+ i 1)))
                     (a0 (belt-route-entry-angle node prev))
                     (da (* (nth node 3) (/ (belt-route-sweep node prev next) k))))
                 (map (lambda (j) (vec (+ (nth node 0) (* (nth node 2) (cos (+ a0 (* da j)))))
                                       (+ (nth node 1) (* (nth node 2) (sin (+ a0 (* da j)))))))
                      (range (+ k 1)))))))
    (concat (vec start) (flatten (map arc (range (len circles)))) (vec end))))

; The pitch path as a std path (z = 0, up +z): each run a line, each wrap
; an arc.
(fn belt-route-path (b start circles end)
  (let ((nodes (belt-route-nodes b start circles end)) (segs (belt-route-segments nodes))
        (line (lambda (sg) (path-line (vec (nth sg 0) (nth sg 1) 0) (vec (nth sg 2) (nth sg 3) 0))))
        (arc (lambda (i)
               (let ((node (nth nodes (+ i 1))) (prev (nth segs i)) (next (nth segs (+ i 1))) (r (nth node 2)))
                 (path-arc (vec (nth node 0) (nth node 1) 0) (vec 0 0 (nth node 3)) r
                           (vec (/ (- (nth prev 2) (nth node 0)) r) (/ (- (nth prev 3) (nth node 1)) r) 0)
                           (belt-route-sweep node prev next))))))
    (path (concat (fold (lambda (acc i) (concat acc (vec (line (nth segs i)) (arc i)))) (vec) (range (len circles)))
                  (vec (line (nth segs -1))))
          (vec 0 0 1))))

; The open belt, centred on z = 0: belt-along on its pitch path, teeth from
; start to end, free ends cut flat, square to the path. shift moves its
; teeth along it, not the path's length or envelope. ink check reports
; where the pitch path runs through itself or bends tighter than the belt
; reaches on its tooth side (pitch-line differential plus tooth height).
(fn belt-route-side (start circles end)
  (let ((c (nth circles 0))) (* (nth c 3) (if (> (nth c 4) 0) 1 -1))))
(fn belt-route (b start circles end width shift)
  (check (path-thickness (belt-route-points b start circles end 16) >= (+ (get b "pitch-line-differential") (get b "tooth-height")))
    (belt-along b (belt-route-path b start circles end) width (belt-route-side start circles end) shift)))
; The same as a std carried profile, for (attr carried ...) beside its shape.
(fn belt-route-carried (b start circles end width shift)
  (belt-carried b (belt-route-path b start circles end) width (belt-route-side start circles end) shift))

;; ---------------------------------------------------------------- a roller chain on two sprockets

; A roller chain of pitch p on two sprockets, A (ta teeth) at the origin
; and B (tb) at (c, 0). roller-chain-path is its pins' path, a std
; path-loop (z = 0, up +z) counter-clockwise from A's upper tangent point:
; round A, the lower run to B, round B and back. Round a sprocket of N
; teeth it follows the pitch circle, radius p / (2 sin (180 / N)),
; advancing N p / 2 pi mm along the path per radian, so pins p apart sit on
; the pitch polygon's corners. On the runs they stay on the tangent line,
; without chordal rise and fall. Link i of the chain, moved shift along
; it, is (path-link-at P (+ (* i p) shift) p). roller-chain-phase-a and
; -b turn the sprockets (degrees) to keep their gaps, on +y unturned, under
; the pins. Equal sprockets at c = m p close a chain of 2 m + N links
; exactly. Check others with path-length.
(fn roller-chain-r (p n) (/ (/ p (sin (/ pi n))) 2))
(fn roller-chain-theta (p ta tb c)                ; angle of A's upper tangent point, radians
  (acos (/ (- (roller-chain-r p ta) (roller-chain-r p tb)) c)))
(fn roller-chain-path (p ta tb c)
  (let ((ra (roller-chain-r p ta)) (rb (roller-chain-r p tb))
        (th (roller-chain-theta p ta tb c)) (ct (cos th)) (st (sin th)))
    (path-loop (vec (path-arc (vec 0 0 0) (vec 0 0 1) ra (vec ct st 0) (- (* 2 pi) (* 2 th)) (/ (* ta p) (* 2 pi)))
                    (path-line (vec (* ra ct) (* ra (neg st)) 0) (vec (+ c (* rb ct)) (* rb (neg st)) 0))
                    (path-arc (vec c 0 0) (vec 0 0 1) rb (vec ct (neg st) 0) (* 2 th) (/ (* tb p) (* 2 pi)))
                    (path-line (vec (+ c (* rb ct)) (* rb st) 0) (vec (* ra ct) (* ra st) 0)))
               (vec 0 0 1))))
(fn roller-chain-phase-a (p ta tb c shift)
  (* (/ 180 pi) (+ (- (roller-chain-theta p ta tb c) (* 0.5 pi)) (/ (* 2 pi shift) (* ta p)))))
(fn roller-chain-phase-b (p ta tb c shift)
  (let ((th (roller-chain-theta p ta tb c))
        (s-b0 (+ (* (/ (* ta p) (* 2 pi)) (- (* 2 pi) (* 2 th))) (* c (sin th)))))   ; station at B's lower tangent point
    (* (/ 180 pi) (- (/ (* 2 pi (- shift s-b0)) (* tb p)) th (* 0.5 pi)))))

;; ---------------------------------------------------------------- patterns

; n copies of f, turned evenly about z or stepped by (dx dy dz).
(fn pattern-polar (f n) (reduce min (map (lambda (i) (rotate-z f (/ (* 360 i) n))) (range n))))
(fn pattern-linear (f n dx dy dz) (reduce min (map (lambda (i) (move f (* i dx) (* i dy) (* i dz))) (range n))))
; n points evenly on a circle of radius r about z, the first deg from +x.
(fn polar-points (n r deg)
  (map (lambda (i) (let ((a (deg->rad (+ deg (/ (* 360 i) n))))) (vec (* r (cos a)) (* r (sin a)) 0))) (range n)))
; The four corners (+-hx, +-hy, 0) of a rectangle.
(fn rect-points (hx hy) (vec (vec hx hy 0) (vec (neg hx) hy 0) (vec (neg hx) (neg hy) 0) (vec hx (neg hy) 0)))

;; ---------------------------------------------------------------- thread frames
;
; A thread frame is one in which a thread is thread-iso's right-hand helix
; about z, in its phase. Turned half over about x it is still one. Tapped
; holes' frames (tap-0 ...) and nuts' and screws' seats are thread frames,
; so a screw seated face to face on a tap's frame, or a nut on a screw's,
; is in phase however the tapped part was turned. Seated h further out
; (past a washer or a plate), place it at (at (along-thread h p) F): moved
; h along z and turned with the helix. A part places its taps at its
; frames, (at-frames tap frames), and declares the same frames. Turned
; copies of a tap (a bolt circle) keep phase, which follows the frame. A
; mirror is not a frame: a mirrored tap is left-hand.
(fn along-thread (h p) (at (rot-z (* 360 (/ h p))) (move 0 0 h)))
; The ring about z where an ISO thread of size d and pitch p meets its
; mate: from just inside the external thread's root (d/2 - 0.613p) to past
; the internal one's (d/2). It is what turns into the mate, not what is
; pushed.
(fn thread-band (d p) (max (- (len2 x y) (+ (* 0.5 d) (* 0.5 p))) (- (- (* 0.5 d) (* 0.65 p)) (len2 x y))))
(fn at-frames (f frames) (reduce min (map (lambda (F) (place f F)) frames)))
; Unturned frames at each point (vec x y z): copies of f at points are
; (at-frames f (point-frames points)).
(fn point-frames (points) (map (lambda (q) (at (move (nth q 0) (nth q 1) (nth q 2)))) points))

;; ---------------------------------------------------------------- build order
;
; A build order is a vec of steps, the parts in the order they go in.
; (build-step part from into ...): the part (a path under the assembly, as
; text), the unit vec it comes in from, and the parts its thread turns into
; (a T-nut, a tapped part). (build-laid part): a part laid in place by hand
; (a belt), with no straight path. A step's path starts (build-length s
; part box) mm out along from, where the part's box is past that of the
; parts placed before it: as if from infinity. box is their (bounds ...),
; found once as the model loads: it is refused for parts that move with t
; and reads ranged params at their set value, so a build order is for an
; assembly at rest. Step k is buildable if its part, moved by every u in
; 0..1 of that path ((build-along f s length u), u a ranged param), is
; disjoint from the parts placed before it but those it turns into (the
; assembly except (build-left-out order k)), and from each of those but
; the thread's band about its axis (thread-band on its drive frame).
; Outside the band its core and head meet that part as any other part
; would, so a screw longer than its blind tap, or coming in off its
; thread's axis, fails.
(fn build-step (part from into ...) (vec part from into))
(fn build-laid (part) (vec part (vec 0 0 0) (vec)))
(fn build-later (order k) (map (lambda (s) (nth s 0)) (slice order k)))
(fn build-into (order k) (nth (nth order k) 2))
(fn build-left-out (order k) (concat (build-later order k) (build-into order k)))
(fn build-length (s part box)
  (let ((v (nth s 1)) (p (bounds part))
        (along (lambda (b f) (reduce + (map (lambda (i) (f (* (nth v i) (nth (nth b 0) i)) (* (nth v i) (nth (nth b 1) i)))) (range 3))))))
    (max 0 (- (along box max) (along p min)))))
(fn build-along (f s length u)
  (let ((d (* u length)) (v (nth s 1))) (move f (* d (nth v 0)) (* d (nth v 1)) (* d (nth v 2)))))

;; ---------------------------------------------------------------- text
;
; Single-stroke (engraving) text: a 2D field in (x, y), the exact distance
; to the round-capped strokes, so it extrudes, engraves, embosses and
; offsets like any exact 2D profile. Inside, where strokes cross, it is a
; lower bound on the distance to the edge. It promises slope 1
; (`ink spec --section "The one promise"`).
;
;   (text-2d "ABC 123" h)    cap height h, stroke h/10, centred on the origin
;   (text-2d-styled text h stroke spacing ax ay)
;       stroke: stroke width. spacing: extra gap between letters (0 keeps
;       the font's). ax, ay: the anchor, as fractions of the line's advance
;       width and cap height. (0 0) puts the baseline's left end at the
;       origin, (0.5 0.5) the centre.
;   (text-width text h stroke spacing)   the line's advance width
; Sizes and anchors may read x, y or z, without the promise.
;
; Each advance grows by the stroke, keeping the font's gaps. Ink stays
; within the advance width, and at most stroke/2 below the baseline or
; above the caps, except descenders (g p y, 0.33 h), brackets and some
; punctuation; `bounds` gives a proven box. One line, of printable ASCII
; (space to ~) and the degree sign. Any other character is an error when
; the model loads (index out of range in hershey-glyph).
;
; Font: Hershey simplex roman (futural.jhf, from
; https://github.com/kamalmostafa/hershey-fonts), glyphs verbatim in the
; Hurt format: letter pairs, R the origin, x = code - 82, y = code - 82
; downwards. The first pair is the left and right extent; " R" lifts the
; pen. Its terms ask that these acknowledgements go with the data:
;   - The Hershey Fonts were originally created by Dr. A. V. Hershey while
;     working at the U. S. National Bureau of Standards.
;   - The format of the Font data in this distribution was originally
;     created by James Hurt, Cognition, Inc., 900 Technology Park Drive,
;     Billerica, MA 01821 (mit-eddie!ci-dandelion!hurt).
; and that it is not converted into the U.S. NTIS's format. The degree
; sign is ours, not the font's: an octagon at the top of the caps.

(def hershey-simplex
  (vec
    "JZ"                                                            ; space
    "MWRFRT RRYQZR[SZRY"                                            ; !
    "JZNFNM RVFVM"                                                  ; "
    "H]SBLb RYBRb RLOZO RKUYU"                                      ; #
    "H\\PBP_ RTBT_ RYIWGTFPFMGKIKKLMMNOOUQWRXSYUYXWZT[P[MZKX"       ; $
    "F^[FI[ RNFPHPJOLMMKMIKIIJGLFNFPGSHVHYG[F RWTUUTWTYV[X[ZZ[X[VYTWT" ; %
    "E_\\O\\N[MZMYNXPVUTXRZP[L[JZIYHWHUISJRQNRMSKSIRGPFNGMIMKNNPQUXWZY[[[\\Z\\Y" ; &
    "MWRHQGRFSGSIRKQL"                                              ; '
    "KYVBTDRGPKOPOTPYR]T`Vb"                                        ; (
    "KYNBPDRGTKUPUTTYR]P`Nb"                                        ; )
    "JZRLRX RMOWU RWOMU"                                            ; *
    "E_RIR[ RIR[R"                                                  ; +
    "NVSWRXQWRVSWSYQ["                                              ; ,
    "E_IR[R"                                                        ; -
    "NVRVQWRXSWRV"                                                  ; .
    "G][BIb"                                                        ; /
    "H\\QFNGLJKOKRLWNZQ[S[VZXWYRYOXJVGSFQF"                         ; 0
    "H\\NJPISFS["                                                   ; 1
    "H\\LKLJMHNGPFTFVGWHXJXLWNUQK[Y["                               ; 2
    "H\\MFXFRNUNWOXPYSYUXXVZS[P[MZLYKW"                             ; 3
    "H\\UFKTZT RUFU["                                               ; 4
    "H\\WFMFLOMNPMSMVNXPYSYUXXVZS[P[MZLYKW"                         ; 5
    "H\\XIWGTFRFOGMJLOLTMXOZR[S[VZXXYUYTXQVOSNRNOOMQLT"             ; 6
    "H\\YFO[ RKFYF"                                                 ; 7
    "H\\PFMGLILKMMONSOVPXRYTYWXYWZT[P[MZLYKWKTLRNPQOUNWMXKXIWGTFPF" ; 8
    "H\\XMWPURRSQSNRLPKMKLLINGQFRFUGWIXMXRWWUZR[P[MZLX"             ; 9
    "NVROQPRQSPRO RRVQWRXSWRV"                                      ; :
    "NVROQPRQSPRO RSWRXQWRVSWSYQ["                                  ; ;
    "F^ZIJRZ["                                                      ; <
    "E_IO[O RIU[U"                                                  ; =
    "F^JIZRJ["                                                      ; >
    "I[LKLJMHNGPFTFVGWHXJXLWNVORQRT RRYQZR[SZRY"                    ; ?
    "E`WNVLTKQKOLNMMPMSNUPVSVUUVS RQKOMNPNSOUPV RWKVSVUXVZV\\T]Q]O\\L[JYHWGTFQFNGLHJJILHOHRIUJWLYNZQ[T[WZYYZX RXKWSWUXV" ; @
    "I[RFJ[ RRFZ[ RMTWT"                                            ; A
    "G\\KFK[ RKFTFWGXHYJYLXNWOTP RKPTPWQXRYTYWXYWZT[K["             ; B
    "H]ZKYIWGUFQFOGMILKKNKSLVMXOZQ[U[WZYXZV"                        ; C
    "G\\KFK[ RKFRFUGWIXKYNYSXVWXUZR[K["                             ; D
    "H[LFL[ RLFYF RLPTP RL[Y["                                      ; E
    "HZLFL[ RLFYF RLPTP"                                            ; F
    "H]ZKYIWGUFQFOGMILKKNKSLVMXOZQ[U[WZYXZVZS RUSZS"                ; G
    "G]KFK[ RYFY[ RKPYP"                                            ; H
    "NVRFR["                                                        ; I
    "JZVFVVUYTZR[P[NZMYLVLT"                                        ; J
    "G\\KFK[ RYFKT RPOY["                                           ; K
    "HYLFL[ RL[X["                                                  ; L
    "F^JFJ[ RJFR[ RZFR[ RZFZ["                                      ; M
    "G]KFK[ RKFY[ RYFY["                                            ; N
    "G]PFNGLIKKJNJSKVLXNZP[T[VZXXYVZSZNYKXIVGTFPF"                  ; O
    "G\\KFK[ RKFTFWGXHYJYMXOWPTQKQ"                                 ; P
    "G]PFNGLIKKJNJSKVLXNZP[T[VZXXYVZSZNYKXIVGTFPF RSWY]"            ; Q
    "G\\KFK[ RKFTFWGXHYJYLXNWOTPKP RRPY["                           ; R
    "H\\YIWGTFPFMGKIKKLMMNOOUQWRXSYUYXWZT[P[MZKX"                   ; S
    "JZRFR[ RKFYF"                                                  ; T
    "G]KFKULXNZQ[S[VZXXYUYF"                                        ; U
    "I[JFR[ RZFR["                                                  ; V
    "F^HFM[ RRFM[ RRFW[ R\\FW["                                     ; W
    "H\\KFY[ RYFK["                                                 ; X
    "I[JFRPR[ RZFRP"                                                ; Y
    "H\\YFK[ RKFYF RK[Y["                                           ; Z
    "KYOBOb RPBPb ROBVB RObVb"                                      ; [
    "KYKFY^"                                                        ; \
    "KYTBTb RUBUb RNBUB RNbUb"                                      ; ]
    "JZRDJR RRDZR"                                                  ; ^
    "I[Ib[b"                                                        ; _
    "NVSKQMQORPSORNQO"                                              ; `
    "I\\XMX[ RXPVNTMQMONMPLSLUMXOZQ[T[VZXX"                         ; a
    "H[LFL[ RLPNNPMSMUNWPXSXUWXUZS[P[NZLX"                          ; b
    "I[XPVNTMQMONMPLSLUMXOZQ[T[VZXX"                                ; c
    "I\\XFX[ RXPVNTMQMONMPLSLUMXOZQ[T[VZXX"                         ; d
    "I[LSXSXQWOVNTMQMONMPLSLUMXOZQ[T[VZXX"                          ; e
    "MYWFUFSGRJR[ ROMVM"                                            ; f
    "I\\XMX]W`VaTbQbOa RXPVNTMQMONMPLSLUMXOZQ[T[VZXX"               ; g
    "I\\MFM[ RMQPNRMUMWNXQX["                                       ; h
    "NVQFRGSFREQF RRMR["                                            ; i
    "MWRFSGTFSERF RSMS^RaPbNb"                                      ; j
    "IZMFM[ RWMMW RQSX["                                            ; k
    "NVRFR["                                                        ; l
    "CaGMG[ RGQJNLMOMQNRQR[ RRQUNWMZM\\N]Q]["                       ; m
    "I\\MMM[ RMQPNRMUMWNXQX["                                       ; n
    "I\\QMONMPLSLUMXOZQ[T[VZXXYUYSXPVNTMQM"                         ; o
    "H[LMLb RLPNNPMSMUNWPXSXUWXUZS[P[NZLX"                          ; p
    "I\\XMXb RXPVNTMQMONMPLSLUMXOZQ[T[VZXX"                         ; q
    "KXOMO[ ROSPPRNTMWM"                                            ; r
    "J[XPWNTMQMNNMPNRPSUTWUXWXXWZT[Q[NZMX"                          ; s
    "MYRFRWSZU[W[ ROMVM"                                            ; t
    "I\\MMMWNZP[S[UZXW RXMX["                                       ; u
    "JZLMR[ RXMR["                                                  ; v
    "G]JMN[ RRMN[ RRMV[ RZMV["                                      ; w
    "J[MMX[ RXMM["                                                  ; x
    "JZLMR[ RXMR[P_NaLbKb"                                          ; y
    "J[XMM[ RMMXM RM[X["                                            ; z
    "KYTBRCQDPFPHQJRKSMSOQQ RRCQEQGRISJTLTNSPORSTTVTXSZR[Q]Q_Ra RQSSUSWRYQZP\\P^Q`RaTb" ; {
    "NVRBRb"                                                        ; |
    "KYPBRCSDTFTHSJRKQMQOSQ RRCSESGRIQJPLPNQPURQTPVPXQZR[S]S_Ra RSSQUQWRYSZT\\T^S`RaPb" ; }
    "F^IUISJPLONOPPTSVTXTZS[Q RISJQLPNPPQTTVUXUZT[Q[O"))            ; ~
(def hershey-degree "MWUITGRFPGOIPKRLTKUI")       ; the degree sign (U+00B0)

; A glyph by code point: (vec left right segments), each segment (vec a b),
; in font units (cap height 21, baseline y = 0, y up). Control characters
; (below 32) fail like any missing character, not wrap to the table's end.
(fn hershey-glyph (c)
  (let ((g (char-codes (if (== c 176) hershey-degree (nth hershey-simplex (if (< c 32) 1000 (- c 32))))))
        (pairs (map (lambda (i) (vec (nth g (* 2 i)) (nth g (+ 1 (* 2 i))))) (range 1 (/ (len g) 2))))
        (strokes (fold (lambda (acc p)
                         (if (== (nth p 0) 32)
                           (vec (concat (nth acc 0) (vec (nth acc 1))) (vec))
                           (vec (nth acc 0) (concat (nth acc 1) (vec (vec (- (nth p 0) 82) (- 91 (nth p 1))))))))
                       (vec (vec) (vec)) pairs))
        (strokes (concat (nth strokes 0) (vec (nth strokes 1)))))
    (vec (- (nth g 0) 82) (- (nth g 1) 82)
         (flatten (map (lambda (st) (map (lambda (a b) (vec a b)) (slice st 0 -1) (slice st 1))) strokes)))))

(fn text-advances (gs h stroke spacing)
  (map (lambda (g) (+ (* (/ h 21) (- (nth g 1) (nth g 0))) stroke spacing)) gs))
(fn text-width (s h stroke spacing)
  (- (reduce + (text-advances (map hershey-glyph (char-codes s)) h stroke spacing)) spacing))

(fn text-2d-styled (s h stroke spacing ax ay)
  (let ((k (/ h 21))
        (gs (map hershey-glyph (char-codes s)))
        (advs (text-advances gs h stroke spacing))
        (w (- (reduce + advs) spacing))
        (x0s (fold (lambda (acc a) (concat acc (vec (+ (nth acc -1) a)))) (vec (* (neg ax) w)) (slice advs 0 -1)))
        (y0 (* (neg ay) h))
        (ds (flatten (map (lambda (g x0)
                            (let ((o (vec (+ x0 (* 0.5 stroke) (* (neg k) (nth g 0))) y0))
                                  (segs (nth g 2)))
                              (if (== (len segs) 0) (vec)
                                (vec (reduce min (map (lambda (sg) (segment-distance-at here2 (+ o (* k (nth sg 0))) (+ o (* k (nth sg 1))))) segs))))))
                          gs x0s))))
    (lipschitz 1 (- (reduce min ds) (* 0.5 stroke)))))

(fn text-2d (s h) (text-2d-styled s h (* 0.1 h) 0 0.5 0.5))

; ---------------------------------------------------------------- lattices (TPMS)
;
; Fields zero on a triply periodic surface of the given cell period,
; divided by their greatest gradient, sqrt(3) k (k = 2 pi / period), so
; their slope is at most 1. A `shell` of one is at least as thick as asked
; (thicker where the gradient is less).

; It promises slope 1 + 1e-15, allowing for float rounding.
(fn gyroid (period)
  (let ((k (/ (* 2 pi) period)))
    (lipschitz (+ 1 1e-15)
      (/ (+ (* (sin (* k x)) (cos (* k y)))
            (* (sin (* k y)) (cos (* k z)))
            (* (sin (* k z)) (cos (* k x))))
         (* k (sqrt 3))))))

; A sheet lattice with walls th thick, filling region.
(fn lattice-fill (region lattice th) (max region (shell lattice th)))

; ---------------------------------------------------------------- manufacturing
;
; Claims that a part can be made. Each process's published minimums, the
; numbers to claim, are in `ink spec --section "Manufacturing checks"`.
;
; The air around f, out to margin past its bounding box.
; (min-thickness (air-around part 2) > 0.5) holds where no gap, slot or
; hole is narrower than 0.5 mm. Only gaps read thinner than margin, so
; choose a margin above the gap claimed. An estimate, by rays like
; min-thickness.
(fn air-around (f margin)
  (let ((b (bounds f)) (lo (nth b 0)) (hi (nth b 1)))
    (subtract (emax (slab-x (- (nth lo 0) margin) (+ (nth hi 0) margin))
                    (slab-y (- (nth lo 1) margin) (+ (nth hi 1) margin))
                    (slab-z (- (nth lo 2) margin) (+ (nth hi 2) margin)))
              f)))
; The air in region at least lo clear of f along each axis: where f and
; its six copies moved lo along +-x, +-y and +-z are all outside. A hidden
; pocket that holds a ball of radius lo has its centre in it. Only signs
; are read, so any field of the shape gives the same set.
(fn clear-air (f region lo)
  (max region (neg (min f (move f lo 0 0) (move f (neg lo) 0 0) (move f 0 lo 0) (move f 0 (neg lo) 0)
                            (move f 0 0 lo) (move f 0 0 (neg lo))))))
; f read s ahead along the unit direction (dx dy dz): inside where the
; point s along it from here is material.
(fn ahead (f dx dy dz s) (move f (* -1 s dx) (* -1 s dy) (* -1 s dz)))
; Line of sight: the points of band with material of f ahead along every
; view, a direction (vec dx dy dz). stock holds the part (a mill's block, a
; mould's box). A view looks s times the stock across (its size dotted
; with |dx dy dz|), so at s = 1 it passes the part's far side from
; anywhere in band. Material of f outside the stock blocks too, within
; reach. Each view's s is a param the claim searches over 0 to 1, named
; after it: t-slot-reaches 0, t-slot-reaches 1, ...
; (empty (hidden part stock (clear-air part stock 0.5) views)) holds where
; every point of the stock's air sees out along some view, for every s.
; Then nothing is undercut from those directions (a mill's setups, a
; mould's pull, a spray, a probe). A failure's witness is a hidden point
; and each s where its view is blocked. Only signs are read, so it holds
; for any field of the shape, toleranced too. Undercuts shallower than lo
; plus twice the claim's accuracy (how far off it may be) pass. Sight
; lines have no width: whether a cutter fits is not checked. Each view's
; param multiplies the search: a claim per setup, its band cut to what
; that setup reaches, may decide where one over all views does not.
(fn hidden (f stock band views)
  (let ((b (bounds stock))
        (size (map (lambda (k) (- (nth (nth b 1) k) (nth (nth b 0) k))) (range 3)))
        (reaches (map (lambda (v) (if (== (len v) 3) (ranged 0 1)
                                      (error "a view of hidden is a direction (vec dx dy dz): hidden ranges how far it looks itself; drop its 4th element")))
                      views)))
    (apply max (concat (vec band)
      (map (lambda (v s) (let ((d (unit v)))
                           (ahead f (nth d 0) (nth d 1) (nth d 2) (* (dot (map abs d) size) s))))
           views reaches)))))

; Graded lattices: the wall thickness th giving a sheet lattice relative
; density rho (material over cell volume, 0 to 1; may be a field), for
; (lattice-fill region (gyroid period) (gyroid-thickness rho period)).
; The shell is |F| <= c, F the surface's function (sin a cos b + ... for
; the gyroid), c = pi sqrt(3) th / period. c is a cubic in rho through 0
; (tpms-thickness's a1 a2 a3), fit to the share of a cell where |F| <= c on
; a 240^3 grid, up to rho = most (0.95); rho is clamped to 0..most. The fit
; is within 0.0011 of rho at 0.95, 0.0003 elsewhere.
; Li, Liao, Dai and Xie (Materials 12, 2183, 2019) fit the gyroid sheet as
; rho = 0.675 c + 0.012 from meshed samples, 0.015 to 0.03 above this.
; Stiffness goes as about rho^n at low density: sheet gyroid 0.596
; rho^1.467 to rho 0.45, steeper above (Nguyen-Xuan, Tran, Thai and Lee,
; arXiv 2303.13838, 2023, homogenized cells).
(fn tpms-thickness (rho period most a1 a2 a3)
  (let ((r (clamp rho 0 most)))
    (/ (* period r (+ a1 (* r (+ a2 (* r a3))))) (* pi (sqrt 3)))))
(fn gyroid-thickness (rho period) (tpms-thickness rho period 0.95 1.5469 0.0361 -0.1249))

; ---------------------------------------------------------------- printing
;
; Overhangs of a part built up along +z on a bed at z0.
; (overhangs f z0 a h) is the air under faces leaning more than a degrees
; past vertical (0 < a < 90), read over a height h: points with material
; h above and none h tan a away at their height, in any of 8 directions
; 45 degrees apart.
; (empty (overhangs part 0 45 1)) holds where no face leans more than
; asin(2 e cos a / h) past the limit, e the claim's accuracy (how far off
; it may be): at a 45, h 1, e 0.05 that is 4.1 degrees (49 may hold); at
; h 4, 1.0. A face up to 2.3 degrees under the limit at 45 may fail,
; between the 8 directions; nearer the limit on either side it may be
; unknown. Only signs are read, so it holds for any field of the shape,
; and it decides along edges and creases. A ceiling narrower than
; 2 h tan a, or a flat overhang reaching out less than h tan a, counts as
; held. A bridge counts as an overhang; its span is not checked.
(fn overhangs (f z0 a h)
  (let ((r (* h (tan (deg->rad a))))
        (beside (reduce min (map (lambda (i) (let ((u (deg->rad (* 45 i)))) (move f (* r (cos u)) (* r (sin u)) 0)))
                                 (range 8)))))
    (max (neg f) (ahead f 0 0 1 h) (neg beside) (- z0 z))))

; ---------------------------------------------------------------- sheet metal
;
; A sheet-metal part designed flat, as cut, and folded from that one source.
; (sheet outline holes th k bends): outline a 2D field (std's polygon),
; holes a vec of 2D fields cut through ((vec) for none), thickness th > 0,
; and K-factor k: where the neutral surface (the layer a bend does not
; stretch) lies, a fraction 0..1 of th from the inside face. Geomiq
; tabulates 0.33 to 0.5 by radius, process and material.
; (sheet-bend a b angle r) folds what lies left of its line a -> b by angle
; degrees to inside radius r > 0: positive toward +z (the top face, z =
; th), at most 180 either way. a and b are (vec u v) points on the flat
; where the bent material ends, at the outline or a relief.
; Its zone, the strip that bends, is a -> b long and BA wide about the line.
; BA, the bend allowance, is the neutral surface's arc: angle (r + k th),
; angle in radians (SOLIDWORKS Help, K-Factor). A flange, all of the flat
; beyond its bend's zone, may be wider than the bend. A bend whose line
; lies beyond another's zone is on its flange and folds with it.
; sheet refuses a = b, crossing bends (a bend's line through another's
; zone), two bends each beyond the other, and a bend's line beyond two
; bends neither beyond the other; sheet-frame refuses such a point.
; Closed corners and tabs bent out of a plate's interior are not built;
; sheet-unrelieved fails where a design needs them.
;   (sheet-folded S)  the folded part: the flat at z 0..th where no bend
;     moved it, flanges turned rigidly about their bends' axes, zones rolled
;     round theirs, r inside to r + th outside. A point keeps its place
;     along the neutral surface; plane sections stay plane. Flanges reach
;     1/100 th (at most BA/4) into their zones, sticking out of the
;     bend's outside face there by under (th/100)^2 / 2(r + th).
;   (sheet-flat S)  the flat pattern, z 0..th. `ink slice --out FILE.dxf`
;     writes its cut file: outline and holes, and the bend lines its body
;     marks, `(mark bend (sheet-bend-lines S))`.
;   (sheet-bend-lines S)  each bend's centre line, a to b within the
;     outline: marks (vec g r) for a body's (mark NAME ...).
;   (sheet-frame S u v side)  the frame at the flat's point (u v) on its top
;     (side 1) or bottom face (-1), carried to the folded part: z out of
;     that face, x the flat's x.
;   Rules: regions a design leaves empty, `(empty RULE)`. A failing
;   claim's witness point shows where.
;   (sheet-short-flanges S 4)  flanges under 4 th from the bend's mold line
;     (where the outside faces meet, extended) anywhere along it, either
;     side (Protolabs; mold line: American Steel & Aluminum), past 90
;     degrees as at 90. It is air outside the outline, or another bend's
;     zone, within 4 th - (r + th) tan(min(angle, 90) / 2) of each zone,
;     and none from 90 degrees at r >= 3 th.
;   (sheet-holes-near-bends S 2.5)  holes within 2.5 th of a bend's zone,
;     beside it (American Steel & Aluminum, Geomiq: 2.5 th + r from the
;     adjacent flange's inside face, at 90 degrees 2.5 th past the zone).
;   (sheet-unrelieved S)  material:
;     - within th beside each end of a bend, from th + r behind its zone to
;       its far edge. A relief must be at least th wide and th + r deep
;       past the bend (American Steel & Aluminum, Geomiq), judged as that
;       rectangle;
;     - past a bend's line beside the bend, on the flange it is on, which
;       would fold with it;
;     - two bends both bend: a corner's zones meeting, a flange beyond
;       another bend's line, bends closer than their zones (a zone off its
;       flange), a lip's flat joining its parent's zone (a tear).
;   (sheet-collisions S)  where folded pieces no bend joins meet (walls
;     through each other, a curl's end on its start). The folded body joins
;     them.
;   (sheet-min-radius S)  the least inside radius, against the material's:
;     `(value (sheet-min-radius S) >= th)` (Geomiq; tighter for soft alloys).
(fn sheet-bend (a b angle r) (vec a b angle r))
(fn sheet-bend-allowance (angle r th k) (* (deg->rad (abs angle)) (+ r (* k th))))
(fn sheet (outline holes th k bends)
  (if (<= th 0) (error "sheet: the thickness th must be above 0")
    (if (+ (< k 0) (> k 1))
      (error "sheet: the K-factor k is where the neutral surface lies, a fraction of th from the bend's inside face: 0 to 1 (Geomiq tabulates 0.33 to 0.5)")
      (let ((gs (map (lambda (bd) (sheet-geom bd th k)) bends))
            (n (len gs))
            ; row i, column j: 1 where bend j's line lies beyond bend i's zone
            (beyond (map (lambda (gi) (map (lambda (gj) (sheet-is-beyond gi (sheet-mid gj))) gs)) gs))
            (crossing (fold + 0 (map (lambda (i) (fold + 0 (map (lambda (j)
              (if (== i j) 0 (sheet-crosses (nth gs i) (nth gs j)))) (range n)))) (range n))))
            (mutual (fold + 0 (map (lambda (i) (fold + 0 (map (lambda (j)
              (if (< i j) (* (nth (nth beyond i) j) (nth (nth beyond j) i)) 0)) (range n)))) (range n))))
            (depths (map (lambda (j) (fold + 0 (map (lambda (row) (nth row j)) beyond))) (range n))))
        (if (> crossing 0)
          (error "sheet: a bend's line crosses another bend's zone: crossing bends are not built (end the one at a relief beside the other, or make it two bends)")
          (if (> mutual 0)
            (error "sheet: two bends each lie beyond the other: a bend folds what lies left of a -> b, so swap a and b of the one that should fold away from the other")
            (vec outline holes th k gs
                 (map (lambda (j) (sheet-order (fold (lambda (acc i) (if (nth (nth beyond i) j) (concat acc (vec i)) acc)) (vec) (range n))
                                               depths "sheet: a bend's line")) (range n))
                 (if (> (len holes) 0) (apply subtract (concat (vec outline) holes)) outline)
                 depths)))))))
(fn sheet-flat (S) (extrude-z (nth S 6) 0 (nth S 2)))
(fn sheet-bend-lines (S)
  (map (lambda (g) (vec (sheet-across g here2) (max (sheet-lat g here2) (nth S 0)))) (nth S 4)))
(fn sheet-folded (S) (let ((ps (sheet-pieces S))) (reduce min (concat (vec (nth ps 0)) (nth ps 1) (nth ps 2)))))
; Tests only pairs of pieces no bend joins (a zone joins its flange and the
; one it is on) whose proven boxes overlap, boxes taken without the holes.
(fn sheet-collisions (S)
  (let ((ps (sheet-pieces S))
        (ws (sheet-pieces (vec (nth S 0) (vec) (nth S 2) (nth S 3) (nth S 4) (nth S 5) (nth S 0) (nth S 7))))
        (chains (nth S 5)) (n (len chains))
        (box (lambda (f) (bounds f (vec (vec 1 1 1) (vec -1 -1 -1)))))
        ; (vec field node-or-bend is-zone box); node -1 is the base
        (pieces (concat
          (map (lambda (a) (vec (if (< a 0) (nth ps 0) (nth (nth ps 1) a)) a 0
                                (box (if (< a 0) (nth ws 0) (nth (nth ws 1) a))))) (concat (vec -1) (range n)))
          (map (lambda (i) (vec (nth (nth ps 2) i) i 1 (box (nth (nth ws 2) i)))) (range n))))
        (m (len pieces))
        (joined (lambda (p q) (* (== (nth p 2) 0) (== (nth q 2) 1)
                                 (+ (== (nth p 1) (nth q 1)) (== (nth p 1) (sheet-parent chains (nth q 1)))))))
        (meet (lambda (p q) (let ((a (nth p 3)) (b (nth q 3)))
                (fold * 1 (map (lambda (k) (* (< (nth (nth a 0) k) (nth (nth b 1) k)) (< (nth (nth b 0) k) (nth (nth a 1) k))))
                               (range 3)))))))
    (reduce min (concat (vec 1) (flatten (map (lambda (i) (map (lambda (j)
      (let ((p (nth pieces i)) (q (nth pieces j)))
        (if (* (< i j) (- 1 (joined p q)) (meet p q)) (max (nth p 0) (nth q 0)) 1)))
      (range m))) (range m)))))))
; The folded pieces: (vec base flanges zones), each placed.
(fn sheet-pieces (S)
  (let ((th (nth S 2)) (gs (nth S 4)) (chains (nth S 5))
        (body (extrude-z (nth S 6) 0 th))
        (cut (lambda (node) (map (lambda (c) (neg (sheet-takes (nth gs c) here2 (nth (nth gs c) 11))))
                                 (sheet-children chains node))))
        (flange (lambda (i)
          (let ((chain (concat (nth chains i) (vec i))))
            (sheet-xform (reduce max (concat (vec body)
                           (map (lambda (a) (sheet-beyond (nth gs a) here2 (nth (nth gs a) 11))) chain)
                           (cut i)))
                         (sheet-xf gs chain)))))
        (zone (lambda (i) (sheet-xform (sheet-bend-zone S i) (sheet-xf gs (nth chains i))))))
    (vec (tag flange (reduce max (concat (vec body) (cut -1))))
         (map (lambda (i) (tag flange (flange i))) (range (len gs)))
         (map (lambda (i) (tag bend (zone i))) (range (len gs))))))
(fn sheet-frame (S u v side)
  (let ((gs (nth S 4)) (p (vec u v)) (xf (sheet-xf gs (sheet-order (sheet-beyond-of gs p) (nth S 7) "sheet-frame: the point"))) (m (nth xf 0))
        (o (+ (sheet-mul m (vec u v (if (> side 0) (nth S 2) 0))) (nth xf 1)))
        (ex (nth m 0)) (ez (* (nth m 2) (if (> side 0) 1 -1))))
    (if (> (fold + 0 (map (lambda (g) (sheet-in-zone g p)) gs)) 0)
      (error "sheet-frame: the point is in a bend's zone: put a frame on a flange")
      (axes (origin (nth o 0) (nth o 1) (nth o 2)) (x (nth ex 0) (nth ex 1) (nth ex 2))
            (z (nth ez 0) (nth ez 1) (nth ez 2))))))
(fn sheet-min-radius (S) (reduce min (map (lambda (g) (nth g 3)) (nth S 4))))
(fn sheet-short-flanges (S kf)
  (let ((th (nth S 2)) (gs (nth S 4)) (n (len gs)))
    (max (slab-z 0 th)
      (reduce min (concat (vec 1) (map (lambda (i)
        (let ((g (nth gs i))
              (w (- (* kf th) (* (+ (nth g 3) th) (tan (deg->rad (* 0.5 (min 90 (abs (nth g 2)))))))))
              (s (abs (sheet-across g here2))) (h (* 0.5 (nth g 7))))
          (if (> w 0)
            (max (sheet-lat g here2) (- h s) (- s h w)
                 (reduce min (concat (vec (neg (nth S 0)))
                   (map (lambda (j) (if (== i j) 1 (sheet-bend-strip (nth gs j) here2))) (range n)))))
            1)))
        (range n)))))))
(fn sheet-holes-near-bends (S kh)
  (let ((th (nth S 2)) (holes (nth S 1)))
    (if (== (len holes) 0) 1
      (max (reduce min holes) (slab-z 0 th)
        (reduce min (concat (vec 1) (map (lambda (g)
          (max (sheet-lat g here2) (- (abs (sheet-across g here2)) (* 0.5 (nth g 7)) (* kh th)))) (nth S 4))))))))
(fn sheet-unrelieved (S)
  (let ((th (nth S 2)) (gs (nth S 4)) (chains (nth S 5)) (n (len gs))
        (ends (map (lambda (g)
          (let ((l (sheet-along g here2)) (s (sheet-across g here2)) (bl (nth g 6)) (h (* 0.5 (nth g 7))))
            (max (min (max (- (neg th) l) l) (max (- bl l) (- l bl th)))
                 (- (neg (+ h th (nth g 3))) s) (- s h))))
          gs))
        (beside (map (lambda (i)
          (let ((g (nth gs i)) (s (sheet-across g here2)))
            (reduce max (concat (vec (neg (sheet-lat g here2)) (neg s) (- s (* 0.5 (nth g 7))))
                                (sheet-on-flange gs (nth chains i) here2)))))
          (range n)))
        (pairs (flatten (map (lambda (i) (map (lambda (j)
          (if (sheet-member i (nth chains j))
            ; j's zone off i's flange, and what j takes in i's zone
            (min (max (sheet-bend-strip (nth gs j) here2) (- (sheet-across (nth gs i) here2) (* 0.5 (nth (nth gs i) 7))))
                 (max (sheet-takes (nth gs j) here2 0) (sheet-bend-strip (nth gs i) here2)))
            (if (* (< i j) (- 1 (sheet-member j (nth chains i))))
              (reduce max (concat (vec (sheet-takes (nth gs i) here2 0) (sheet-takes (nth gs j) here2 0))
                                  (sheet-on-flange gs (nth chains i) here2) (sheet-on-flange gs (nth chains j) here2)))
              1)))
          (range n))) (range n)))))
    (max (nth S 6) (slab-z 0 th) (reduce min (concat (vec 1) ends beside pairs)))))

; A bend's numbers (sheet-geom) by index: 0 a, 1 b, 2 angle, 3 r, 4 u
; (unit a -> b), 5 n (left of u, toward what folds), 6 length, 7 BA, 8 the
; angle's sign, 9 a point on its axis (on the zone's start line, r beyond
; the inside face), 10 its fold (the transform of the flat beyond its zone:
; back BA, turned about the axis), 11 how far flanges reach into its zone.
; A transform (vec m t), m its columns, takes p to m p + t. A lineage is
; the bends a point lies beyond, from the base out. A node is a flange, -1
; the unmoved flat. S is (vec outline holes th k bends lineages flat
; depths).
(fn sheet-cs (deg)
  (let ((a (abs deg)) (sg (if (< deg 0) -1 1)))
    (if (== a 90) (vec 0 sg) (if (== a 180) (vec -1 0) (vec (cos (deg->rad a)) (* sg (sin (deg->rad a))))))))
(fn sheet-rodrigues (u deg)
  (let ((cs (sheet-cs deg)) (c (nth cs 0)) (s (nth cs 1)))
    (map (lambda (e) (+ (* e c) (* (cross u e) s) (* u (* (dot u e) (- 1 c)))))
         (vec (vec 1 0 0) (vec 0 1 0) (vec 0 0 1)))))
(fn sheet-mul (m v) (+ (* (nth m 0) (nth v 0)) (* (nth m 1) (nth v 1)) (* (nth m 2) (nth v 2))))
(fn sheet-after (x1 x2)
  (vec (map (lambda (col) (sheet-mul (nth x1 0) col)) (nth x2 0))
       (+ (sheet-mul (nth x1 0) (nth x2 1)) (nth x1 1))))
(fn sheet-xf (gs chain)
  (fold (lambda (acc i) (sheet-after acc (nth (nth gs i) 10)))
        (vec (vec (vec 1 0 0) (vec 0 1 0) (vec 0 0 1)) (vec 0 0 0)) chain))
(fn sheet-xform (f xf)
  (let ((m (nth xf 0)) (d (- here (nth xf 1))))
    (warp f (dot (nth m 0) d) (dot (nth m 1) d) (dot (nth m 2) d))))
(fn sheet-geom (bd th k)
  (if (== (norm (- (nth bd 1) (nth bd 0))) 0) (error "sheet-bend: a bend's a and b must differ: its line runs from a to b")
  (let ((a (nth bd 0)) (b (nth bd 1)) (ang (nth bd 2)) (r (nth bd 3))
        (u (unit (- b a))) (n (vec (neg (nth u 1)) (nth u 0)))
        (ba (sheet-bend-allowance ang r th k))
        (foot (- a (* n (* 0.5 ba))))
        (sg (if (< ang 0) -1 1))
        (p (vec (nth foot 0) (nth foot 1) (+ (if (< ang 0) 0 th) (* sg r))))
        (q (sheet-rodrigues (vec (nth u 0) (nth u 1) 0) ang)))
    (if (<= r 0) (error "sheet-bend: the inside radius must be above 0")
      (if (> (abs ang) 180) (error "sheet-bend: a bend folds at most 180 degrees either way")
        (if (== ang 0) (error "sheet-bend: a bend of 0 degrees folds nothing")
          (vec a b ang r u n (norm (- b a)) ba sg p
               (vec q (- p (sheet-mul q (+ p (* (vec (nth n 0) (nth n 1) 0) ba)))))
               (min (* 0.01 th) (* 0.25 ba)))))))))
(fn sheet-mid (g) (* 0.5 (+ (nth g 0) (nth g 1))))
(fn sheet-across (g p) (dot (nth g 5) (- p (nth g 0))))
(fn sheet-along (g p) (dot (nth g 4) (- p (nth g 0))))
(fn sheet-lat (g p) (let ((l (sheet-along g p))) (max (neg l) (- l (nth g 6)))))
(fn sheet-bend-strip (g p) (max (sheet-lat g p) (- (abs (sheet-across g p)) (* 0.5 (nth g 7)))))
(fn sheet-is-beyond (g p) (> (sheet-across g p) (* 0.5 (nth g 7))))
; Negative on the flange a lineage ends at: beyond each of its bends' zones.
(fn sheet-on-flange (gs chain p) (map (lambda (a) (- (* 0.5 (nth (nth gs a) 7)) (sheet-across (nth gs a) p))) chain))
; 1 where bend j's line a -> b passes through bend i's zone (touching is
; not crossing).
(fn sheet-crosses (gi gj)
  (let ((a (nth gj 0)) (d (- (nth gj 1) (nth gj 0))) (h (* 0.5 (nth gi 7)))
        (s (sheet-t-range (sheet-across gi a) (dot (nth gi 5) d) (neg h) h))
        (l (sheet-t-range (sheet-along gi a) (dot (nth gi 4) d) 0 (nth gi 6))))
    (> (- (min 1 (nth s 1) (nth l 1)) (max 0 (nth s 0) (nth l 0))) 1e-9)))
; The t range where f0 + t df lies strictly between lo and hi.
(fn sheet-t-range (f0 df lo hi)
  (if (< (abs df) 1e-12)
    (if (* (> f0 lo) (< f0 hi)) (vec -1 2) (vec 1 0))
    (let ((t0 (/ (- lo f0) df)) (t1 (/ (- hi f0) df))) (vec (min t0 t1) (max t0 t1)))))
(fn sheet-in-zone (g p) (<= (sheet-bend-strip g p) 0))
(fn sheet-member (i l) (fold (lambda (acc a) (max acc (== a i))) 0 l))
; Negative beyond the bend's zone less e, and on what the bend takes from
; the flange it is on: its zone from e past its start, and all beyond.
(fn sheet-beyond (g p e) (- (- (* 0.5 (nth g 7)) e) (sheet-across g p)))
(fn sheet-takes (g p e)
  (let ((s (sheet-across g p)) (h (* 0.5 (nth g 7))))
    (min (max (sheet-lat g p) (- (+ (neg h) e) s)) (- h s))))
(fn sheet-beyond-of (gs p)
  (fold (lambda (acc i) (if (> (sheet-across (nth gs i) p) (* 0.5 (nth (nth gs i) 7))) (concat acc (vec i)) acc))
        (vec) (range (len gs))))
; A point's lineage: the bends as it lies beyond, by depth (how many bends
; each one's line lies beyond).
(fn sheet-order (as depths what)
  (map (lambda (d)
      (let ((hit (fold (lambda (acc i) (if (== (nth depths i) d) (concat acc (vec i)) acc))
                       (vec) as)))
        (if (== (len hit) 1) (nth hit 0)
          (error (concat what " lies beyond two bends neither of which lies beyond the other: their flanges would overlap on the flat; end one flange before the other's line (a notch), or make one bend fold the other's flange")))))
      (range (len as))))
(fn sheet-parent (chains i) (let ((c (nth chains i))) (if (== (len c) 0) -1 (nth c -1))))
(fn sheet-children (chains node)
  (fold (lambda (acc j) (if (== (sheet-parent chains j) node) (concat acc (vec j)) acc)) (vec) (range (len chains))))
; Bend i's zone in the coordinates of the flange it is on: the point at
; angle psi round the axis from the zone's start, radius rho, reads the
; outline at psi (r + k th) across the zone, rho - r under the inside face.
(fn sheet-bend-zone (S i)
  (let ((th (nth S 2)) (k (nth S 3)) (g (nth (nth S 4) i))
        (r (nth g 3)) (u (nth g 4)) (n (nth g 5)) (ba (nth g 7)) (p (nth g 9))
        (cs (sheet-cs (abs (nth g 2))))
        (v (- here p))
        (vn (+ (* (nth n 0) (nth v 0)) (* (nth n 1) (nth v 1))))
        (ve (* (nth g 8) (nth v 2)))
        (va (+ (* (nth u 0) (nth v 0)) (* (nth u 1) (nth v 1))))
        (rho (len2 vn ve))
        (s (- (* (atan2 vn (neg ve)) (+ r (* k th))) (* 0.5 ba)))
        (m (+ (nth g 0) (* u va) (* n s))))
    (max (warp (nth S 6) (nth m 0) (nth m 1) z)
         (- r rho) (- rho r th) (neg vn) (+ (* vn (nth cs 0)) (* ve (nth cs 1)))
         (max (neg va) (- va (nth g 6))))))

; ---------------------------------------------------------------- load ratings

; A part's load over its record's "load-ratings": 1 at the rating, so a
; static safety factor s holds at a ratio of at most (/ 1 s). w is a
; support's share of the load, (vec fx fy fz mx my mz) in N and N mm in the
; part's own frame, as ext's dynamics gives it: a radial bearing's force
; across its z, a thrust bearing's along it, a guide carriage's across its
; x (its travel) and its moments. spec is the part's record (a catalog
; part's b.spec, imported as b). A missing rating is an error naming it.
(fn load-rating (spec key) (get (get spec "load-ratings") key))
; A rolling bearing's static equivalent load over its static rating C0:
; the share's force over C0 (ISO 76: P0 = Fr on a radial bearing, Fa on a
; thrust bearing).
(fn bearing-static-ratio (w spec) (/ (len3 (nth w 0) (nth w 1) (nth w 2)) (load-rating spec "C0-N")))
; Its basic rating life L10, millions of revolutions: (C / P)^p (ISO 281),
; P = Fr or Fa as above, at least C / 1e6. p is the record's life-exponent:
; 3 for ball, 10/3 for roller bearings. A record without one (a plain
; bearing's, a linear bushing's) is outside ISO 281.
(fn bearing-l10 (w spec)
  (let ((c (load-rating spec "C-N")))
    (pow (/ c (max (len3 (nth w 0) (nth w 1) (nth w 2)) (* 1e-6 c))) (load-rating spec "life-exponent"))))
; A guide carriage's load over its static ratings: the largest of each
; force across its travel over C0 and its moments about x, y and z over
; MR, MP and MY, each on its own (HIWIN Eq. 1.1).
(fn carriage-static-ratio (w spec)
  (let ((c0 (load-rating spec "C0-N")) (m (load-rating spec "moment-rating-Nm")))
    (max (/ (abs (nth w 1)) c0) (/ (abs (nth w 2)) c0) (/ (abs (nth w 3)) (* 1000 (get m "MR")))
         (/ (abs (nth w 4)) (* 1000 (get m "MP"))) (/ (abs (nth w 5)) (* 1000 (get m "MY"))))))
; A belt's tight-side tension in N, installed at t0 and locked (Gates
; Formulas 5 and 6): t0 + |f| / 2, its two equal spans sharing the force f
; it passes, and |f| once the slack span goes slack. f is the torque on a
; pulley it turns over the pitch radius, or its pull on what it moves. A
; constant-force tensioner holding the slack span at t0 gives |f| + t0.
(fn belt-tension (f t0) (max (+ t0 (* 0.5 (abs f))) (abs f)))
; A belt's tension over its working tension, from its record's belt
; section: working-tension-N, else working-tension-N-per-mm times width.
(fn belt-tension-ratio (tension spec)
  (let ((b (get spec "belt")))
    (/ tension (get b "working-tension-N" (* (get b "working-tension-N-per-mm") (get (get spec "geometry") "width"))))))

; ---------------------------------------------------------------- power transmission
;
; (involute-gear-2d z mo pa backlash): an involute spur gear, 2D (extrude
; it), of z teeth, module mo (mm) and pressure angle pa, one tooth centred
; on +x. backlash is its share of the play along the pitch circle (mm).
; Pitch radius mo z / 2. Mating gears sit mo (z1 + z2) / 2 apart, the
; second turned so a gap faces the first's tooth: 180 + 180/z2 deg when
; it sits on the first's +x (mesh-turn gives it in general).
;
; The flank is the true involute, radial below the base circle, and the
; root is a plain circle. The field is the arc length to the flank at each
; radius: exact on the flank, distance-like (gradient near 1) near it.

; inv(a) = tan a - a, the involute function, of tn = tan a. It is a
; series of positive terms in u = tan^2(a/2), which reads tn once, so its
; interval bounds stay tight where tn - atan(tn)'s would not. Its twelve
; terms are within 1e-10 rad for tn <= 1 (a tooth's tip, at 6 teeth or
; more and 20 degrees) and 1.3e-7 at tn = 1.5, one expression for every tn.
(fn involute-of-tan (tn)
  (let ((s (sqrt (+ 1 (square tn)))) (u (- 1 (/ 2 (+ 1 s)))) (w (/ tn (+ 1 s))))
    (* 2 w u (fold (lambda (acc j) (+ (* acc u) 1 (/ (if (== (mod j 2) 0) 1 -1) (+ (* 2 j) 3))))
                   0 (range 11 -1 -1)))))

; (involute-gear-shifted-2d z mo pa backlash shift): the same gear cut
; with its rack moved out shift modules (profile shift, or addendum
; modification). A small pinion with shift > 0 loses its undercut and
; runs at a wider centre distance. Tip and root circles move out shift
; mo. The tooth thickens 2 shift mo tan(pa) at the pitch circle. Where
; the vendor tops the tip, clip it with a circle.
(fn involute-gear-shifted-2d (z mo pa backlash shift)
  (let ((rp (* 0.5 mo z))
        (rb (* rp (cos (deg->rad pa))))
        (ra (+ rp (* mo (+ 1 shift))))
        (rf (- rp (* mo (- 1.25 shift))))
        (r (len2 x y))
        (pitch (/ (* 2 pi) z))
        ; Angle from the nearest tooth's centre line, in [0, pitch/2],
        ; folded inside the gear's material so interval bounds stay tight
        ; in the gaps, where a mating tooth sits.
        (th (- (* 0.5 pitch) (abs (- (mod (atan2 y x) pitch) (* 0.5 pitch)))))
        (tn (/ (sqrt (- (square r) (square rb))) rb))
        (psi (- (+ (/ pi (* 2 z)) (involute-of-tan (tan (deg->rad pa)))
                   (/ (* 2 shift (tan (deg->rad pa))) z))
                (+ (involute-of-tan tn) (/ (* 0.5 backlash) rp))))
        (tooth (max (* r (- th psi)) (- r ra))))
    (min tooth (- r rf))))
(fn involute-gear-2d (z mo pa backlash) (involute-gear-shifted-2d z mo pa backlash 0))
; The rack these gears run on, 2D (extrude it): along x, teeth facing +y,
; pitch line y = 0, a tooth centred on x = 0. It is ISO 53's basic rack:
; straight flanks at pa, pitch pi mo, addendum mo, dedendum 1.25 mo.
; backlash is its share (mm at the pitch line). The body fills all below
; the roots; intersect it with a slab for length and depth. The field is
; the distance to a flank near it.
(fn involute-rack-2d (mo pa backlash)
  (let ((p (* pi mo))
        (u (- (* 0.5 p) (abs (- (mod x p) (* 0.5 p)))))
        (flank (* (- u (- (* 0.25 p) (* 0.5 backlash) (* y (tan (deg->rad pa))))) (cos (deg->rad pa)))))
    (min (max flank (- y mo)) (+ y (* 1.25 mo)))))

; A straight bevel gear's toothed cone, 3D: n teeth, module mo and
; backlash (mm at the pitch circle) at the heel, pressure angle pa, pitch
; angle delta. A pair at 90 degrees has atan(z1 / z2) and its complement.
; The teeth run between cone distances toe and heel from the apex, along
; the pitch cone. The apex is at the origin, the axis z. The gear lies
; below the apex, teeth facing +z, a tooth centred on +x. The pitch radius
; at cone distance a is a sin(delta). Mating gears share the apex, their
; axes delta1 + delta2 apart.
;
; The teeth are Tredgold's form, the one bevel gear generators draw: std's
; involute gear of n / cos(delta) teeth, module mo, on the heel's back
; cone, scaled toward the apex so every cone of the tooth (tip, root,
; flanks) runs to it. The field is distance-like near the teeth. Clip it
; with the gear's blank (a revolve of its section) and cut its bore.
(fn straight-bevel-gear (n mo pa backlash delta toe heel)
  (let ((d (deg->rad delta)) (a (/ (* mo n) (* 2 (sin d))))
        (r (len2 x y)) (rho (len3 x y z))
        (phi (atan2 r (neg z))) (th (* (atan2 y x) (cos d)))
        (rv (* a (+ (tan d) (- phi d))))
        (tooth (warp (involute-gear-2d (/ n (cos d)) mo pa backlash) (* rv (cos th)) (* rv (sin th)) 0))
        (s (- (* r (sin d)) (* z (cos d)))))
    (max (* tooth (/ rho a)) (- s heel) (- toe s))))

; A gear in a mesh: its teeth, or (vec teeth centre-teeth) if
; profile-shifted to run at another count's centre distance. WCP-1008, 10
; teeth at a 12-tooth gear's centres, is (vec 10 12). Its teeth set the
; ratio and phase, its centre teeth where it sits.
(fn gear-teeth (g) (nth (flatten (vec g)) 0))
(fn gear-centre-teeth (g) (let ((v (flatten (vec g)))) (nth v (- (len v) 1))))
; Spur gears in mesh, module mo, za and zb teeth, sit mo (za + zb) / 2
; apart, by their centre teeth.
(fn mesh-distance (za zb mo) (* 0.5 mo (+ (gear-centre-teeth za) (gear-centre-teeth zb))))
; The second gear's turn when the first turns a, the second's centre at
; bearing from the first's, all from one x direction. It turns the other
; way, za / zb as fast, a gap facing the first's tooth on the line of
; centres (each gear has a tooth on its +x).
(fn mesh-turn (za zb bearing a)
  (let ((na (gear-teeth za)) (nb (gear-teeth zb)))
    (+ bearing 180 (/ 180 nb) (/ (* na (- bearing a)) nb))))

; Gears in mesh are relations (relate A B KIND), A the input: gear-mesh
; (spur), bevel-mesh (straight bevel at right angles) and rack-mesh (a
; pinion on a rack). Each reads its gears' axes and places from A's and
; B's frames. An end (PATH FRAME) reads a gear whose axis is not its z at
; FRAME. After teeth and module each optionally takes
; pa [backlash [efficiency]], involute-gear-2d's order:
; - pa: pressure angle, default 20, a number known when the model loads.
;   Under 10, no gear's, is an error: a backlash written where pa goes.
; - backlash: the pair's, both gears' shares (mm at the pitch circle).
; - efficiency: the share of A's power that B receives, the same either
;   way. Not given, unknown.
;
; `ink check` counts a mesh's gears as connected, whatever their backlash.
; A mesh has no zone (a region where overlap is meant), so interfere
; checks the teeth with the rest of the pair: over the pinion's angle,
; where the wheel turns by mesh-turn of its drive (one input, its whole
; span). A wrong ratio or phase is a clash at some angle.
; The claims: at the pitch point the pitch surfaces (circles, cones, the
; rack's pitch line) part by at most the room. They close by at most the
; room or a quarter module, whichever is less. A quarter module is std's
; gears' bottom clearance (ISO 53's basic rack, dedendum 1.25 module);
; past it, tips reach roots. Given a backlash, the room is
; backlash / (2 tan pa). Not given, the backlash is unknown and the room
; is a quarter module both ways: closing, std's clearance; opening, this
; relation's own bound, no standard's. The teeth still mesh there, with
; less contact: 15:45 at 20 degrees has a contact ratio of 1.61 in place,
; 1.37 opened a quarter module. Gears cut to another system get the same
; room (AGMA fine pitch, 20 DP and finer: clearance 0.2 module + 0.05 mm).
; The axes lie as the mesh needs, within an angle set by the room and
; the faces' overlap (mesh-tilt). The teeth's spans along the axis
; overlap (faces-meet). Each span is read from its shape where the mate's
; teeth work (half a module inside the pitch circle to the tips) by
; bounds, a box up to about a 256th of the gear's size too big. Failing,
; the faces are apart; holding, they overlap but for that slack. Its
; margin is the overlap (its upper bound, twice A's face, never binds).
; Each claim is decided to 1e-6 mm (1e-9 for an angle's sine), so a pair
; with no backlash placed at its numbers holds despite rounding.
; Its numbers may read params, but its parts' shapes are read when the
; model loads. If the shapes read a ranged param (one given a range of
; values), its claims over the range are refused. Unless the teeth clash,
; it does not check that the bodies have the teeth it says: a 20:40 pair
; related as 15:45 at their centres passes.
;
; Its coupling, data for packages built on the kernel, is the tooth force
; along either line of action. That runs through the pitch point (mid
; faces' overlap, or a bevel's mean cone distance) at pa to the pitch
; surfaces' common tangent, pushing B from A. Its tangential part is the
; torque over the pitch radius, its separating part that times tan pa
; (AGMA's and ISO 6336's Ft and Fr; a bevel's at its mean radius, AGMA
; 2003's).
;
; A mesh's numbers after teeth and module:
; (vec pa backlash room efficiency clearance). room is how far the pitch
; surfaces may part (mm). efficiency is a vec as a coupling takes it,
; empty where unknown. clearance is the bottom clearance, a quarter module.
(fn mesh-numbers (mo more)
  (let ((n (if (> (len more) 3) (error "a gear mesh takes, after its teeth and module, pa [backlash [efficiency]]") (len more)))
        (given (if (> n 0) (nth more 0) 20))
        (pa (if (< given 10) (error "a pressure angle under 10 degrees, no gear's: a gear mesh's numbers after its module are pa [backlash [efficiency]], involute-gear-2d's order, not the backlash first") given))
        (c (* 0.25 mo)))
    (vec pa (if (> n 1) (nth more 1) "unknown")
         (if (> n 1) (/ (nth more 1) (* 2 (tan (deg->rad pa)))) c)
         (if (> n 2) (vec (nth more 2)) (vec))
         c)))
; mesh-close: the least parting of the pitch surfaces, -min(room,
; clearance) (mm). mesh-tilt: the largest sine of the axes' error,
; sin(room tan pa / face), the angle (the backlash over twice the face)
; capped at a right angle. mesh-face: the faces' overlap, or A's own face
; where they do not overlap (faces-meet fails then).
(fn mesh-close (m) (neg (min (nth m 2) (nth m 4))))
(fn mesh-tilt (m face) (sin (min (* 0.5 pi) (/ (* (nth m 2) (tan (deg->rad (nth m 0)))) face))))
(fn mesh-face (ov span) (let ((f (- (nth ov 1) (nth ov 0)))) (select (neg f) f (- (nth span 1) (nth span 0)))))
; The two lines of action through pitch point p, A's axis ka through oa,
; n the pitch surfaces' normal there from A toward B. Each is at pa to
; their common tangent, the way A's surface moves turning about ka (for
; a rack, its pitch line).
(fn mesh-lines (oa ka p n pa)
  (let ((d (- p oa)) (along (cross ka (unit (- d (* (dot d ka) ka)))))
        (c (cos (deg->rad pa))) (s (sin (deg->rad pa))))
    (vec (vec p (+ (* c along) (* s n))) (vec p (- (* s n) (* c along))))))
; The outer box (bounds) of gear s's teeth, z and mo its teeth and module:
; s within the ring about its z where its mate's teeth work, from half a
; module inside its pitch circle (by centre teeth) to its tip circle,
; which leaves out a shaft or hub. Nothing there is an error: no gear of
; these teeth.
(fn mesh-teeth (s z mo)
  (let ((r (len2 x y)) (rp (* 0.5 mo (gear-centre-teeth z)))
        (ring (max (- r (+ rp mo)) (- (- rp (* 0.5 mo)) r))))
    (bounds (max s ring)
            (error "a gear of a mesh has no material where its mate's teeth work, from half a module inside its pitch circle to its tip circle about its z: is it a gear of that count and module, its axis its z (or read at a frame on its axis, (PATH FRAME))?"))))
; Where two spans overlap along a gear's axis ka through oa, (vec lo hi)
; along ka: its own a, (vec a0 a1) along ka, and another's b, (vec b0 b1)
; along kb from ob.
(fn mesh-overlap (oa ka a ob kb b)
  (let ((t0 (dot (- ob oa) ka)) (u (dot kb ka))
        (e0 (+ t0 (* u (nth b 0)))) (e1 (+ t0 (* u (nth b 1)))))
    (vec (max (nth a 0) (min e0 e1)) (min (nth a 1) (max e0 e1)))))
(fn z-span (b) (vec (nth (nth b 0) 2) (nth (nth b 1) 2)))
(fn gear-mesh (za zb mo more ...)
  (lambda (fa fb sa sb)
    (let ((m (mesh-numbers mo more)) (na (gear-teeth za)) (nb (gear-teeth zb))
          (oa (get fa "origin"))
          (axis (get fa "z"))
          (span (z-span (mesh-teeth sa za mo)))
          (ov (mesh-overlap oa axis span (get fb "origin") (get fb "z") (z-span (mesh-teeth sb zb mo))))
          (face (- (nth ov 1) (nth ov 0)))
          (d (- (get fb "origin") oa))
          (across (- d (* (dot d axis) axis)))
          (apart (norm across))
          (a0 (mesh-distance za zb mo))
          (e (unit across))
          (p (+ oa (* 0.5 (+ (nth ov 0) (nth ov 1)) axis) (* (/ (* na apart) (+ na nb)) e))))
      (relation gear-mesh
        (claim centre-distance apart (+ a0 (mesh-close m)) (+ a0 (nth m 2)) (accuracy 1e-6))
        (claim parallel-axes (norm (cross axis (get fb "z"))) 0 (mesh-tilt m (mesh-face ov span)) (accuracy 1e-9))
        (claim faces-meet face 0 (* 2 (- (nth span 1) (nth span 0))) (accuracy 1e-6))
        (optional (couple (mesh-lines oa axis p e (nth m 0)) 0 (nth m 3)))
        (data teeth (vec na nb) module mo ratio (/ nb na) pressure-angle (nth m 0) backlash (nth m 1))))))
; Straight bevel gears meshing at right angles (std's straight-bevel-gear):
; za and zb plain tooth counts, mo the module at the heel, face (mm) the
; teeth's length along the pitch cone. Each gear's z is its axis, toward
; its apex at its origin. The pitch angles are atan(za / zb) and its
; complement. The pitch point is on A's pitch cone at the mean cone
; distance, toward B. Claims:
; pitch-cones: B's pitch cone's distance from it (the parting). Moving a
;   pinion out along its axis by d parts them by d sin(its pitch angle).
; square-axes: the cosine of the angle between the axes.
; apexes: their distance apart, within the room.
; Its shapes are not read: face stands for their faces.
(fn bevel-mesh (za zb mo face more ...)
  (lambda (fa fb sa sb)
    (let ((plain (if (> (len (flatten (vec za zb))) 2) (error "bevel-mesh takes plain tooth counts, not (vec teeth centre-teeth)") 0))
          (m (mesh-numbers mo more))
          (oa (get fa "origin")) (ka (get fa "z")) (ob (get fb "origin")) (kb (get fb "z"))
          (da (atan2 za zb)) (db (atan2 zb za))
          (rm (- (/ (* 0.5 mo za) (sin da)) (* 0.5 face)))
          (kc (unit (- kb (* (dot kb ka) ka))))
          (p (- oa (* rm (+ (* (cos da) ka) (* (sin da) kc)))))
          (n (- (* (sin da) ka) (* (cos da) kc)))
          (q (- p ob))
          (h (neg (dot q kb)))
          (parting (- (* (norm (- q (* (dot q kb) kb))) (cos db)) (* h (sin db)))))
      (relation bevel-mesh
        (claim pitch-cones parting (mesh-close m) (nth m 2) (accuracy 1e-6))
        (claim square-axes (abs (dot ka kb)) 0 (mesh-tilt m face) (accuracy 1e-9))
        (claim apexes (norm (- ob oa)) 0 (nth m 2) (accuracy 1e-6))
        (optional (couple (mesh-lines oa ka p n (nth m 0)) 0 (nth m 3)))
        (data teeth (vec za zb) module mo face face ratio (/ zb za) pitch-angles (vec (rad->deg da) (rad->deg db))
              pressure-angle (nth m 0) backlash (nth m 1))))))
; A pinion A, z teeth (or (vec z centre-teeth)) of module mo, on a rack B
; along its x, its teeth along its z facing +y, its pitch line its x axis
; (std's involute-rack-2d, extruded along z). Claims:
; pitch-line: the height of the pinion's axis above the rack's pitch
;   plane, its pitch radius within the parting.
; parallel-axes: the sine of the angle between the pinion's axis and the
;   rack's z.
; faces-meet: the least of how far the pinion's teeth overlap the rack's
;   across it, and how far inside the rack's teeth along it the pitch
;   point lies.
; Its lead, the rack's travel per turn of the pinion, is pi mo z.
(fn rack-mesh (z mo more ...)
  (lambda (fa fb sa sb)
    (let ((m (mesh-numbers mo more))
          (oa (get fa "origin")) (axis (get fa "z")) (ob (get fb "origin")) (xb (get fb "x")) (yb (get fb "y"))
          (span (z-span (mesh-teeth sa z mo)))
          (rack (bounds (max sb (- y mo) (- (* -0.5 mo) y))
                        (error "a rack of a mesh has no material where its pinion's teeth work, from half a module below its pitch line to its tip line, along its x with its teeth facing +y")))
          (ov (mesh-overlap oa axis span ob (get fb "z") (z-span rack)))
          (face (- (nth ov 1) (nth ov 0)))
          (r0 (* 0.5 mo (gear-centre-teeth z)))
          (e (unit (- (neg yb) (* (dot (neg yb) axis) axis))))
          (p (+ oa (* 0.5 (+ (nth ov 0) (nth ov 1)) axis) (* 0.5 mo (gear-teeth z) e)))
          (along (dot (- p ob) xb)))
      (relation rack-mesh
        (claim pitch-line (dot (- oa ob) yb) (+ r0 (mesh-close m)) (+ r0 (nth m 2)) (accuracy 1e-6))
        (claim parallel-axes (norm (cross axis (get fb "z"))) 0 (mesh-tilt m (mesh-face ov span)) (accuracy 1e-9))
        (claim faces-meet (min face (- along (nth (nth rack 0) 0)) (- (nth (nth rack 1) 0) along)) 0 (* 2 (- (nth span 1) (nth span 0)))
               (accuracy 1e-6))
        (optional (couple (mesh-lines oa axis p e (nth m 0)) 0 (nth m 3)))
        (data teeth (gear-teeth z) module mo lead (* pi mo (gear-teeth z)) pressure-angle (nth m 0) backlash (nth m 1))))))

; Couplings by a belt, a chain or a screw: relations (relate A B KIND), A
; the input (the motor's side), giving what A and B exert on each other
; as data for packages built on the kernel. Each optionally takes its
; efficiency last: the share of the power A gives that B receives, then
; of what B gives A, by default the same (a screw's: below). Not given, it
; is unknown, as is any effort it moves power for.
; A belt's pretension t0 (N) is installed and locked, as in Gates'
; synchronous drives (Light Power & Precision manual, V.D: an idler locked
; after installation). Its spans' tensions less t0, weighted by length,
; sum to 0: two equal spans carry t0 +- the torque over twice the pitch
; radius (Formulas 5 and 6). A span that would push goes slack. (held t0)
; is a constant-force tensioner (a spring idler) holding the slack span
; at t0.
(fn held (t0) (vec "held" t0))
; Among a belt's or chain's optional arguments, (tensile-stiffness EA)
; makes each span a spring of EA / its length, for packages simulating the
; drive: (belt-span span t0 0.95 (tensile-stiffness 16903)). EA (N) would
; stretch it by its own length: a maker's tension for 0.1% elongation,
; times 1000. Not given, spans are rigid.
(fn tensile-stiffness (ea) (vec "tensile-stiffness" ea))
;
; An open (uncrossed) belt's or chain's two spans on pitch circles ra
; about A's z through its origin and rb about B's: their outer tangents,
; each a coupling line from its point on B's circle toward A's, with its
; length.
(fn open-belt-lines (fa fb ra rb)
  (let ((axis (get fa "z"))
        (ca (get fa "origin"))
        (cb (get fb "origin"))
        (d (- cb ca))
        (e (unit (- d (* (dot d axis) axis))))
        (n (cross axis e))
        (s (/ (- ra rb) (norm (- d (* (dot d axis) axis)))))
        (c (sqrt (- 1 (square s)))))
    (map (lambda (side) (let ((m (+ (* s e) (* side c n))) (q (+ cb (* rb m))) (p (+ ca (* ra m)))) (vec q (unit (- p q)) 0 (norm (- p q))))) (vec 1 -1))))
; A belt drive: an open belt at t0 (N) a span on pulleys of pitch radii ra
; and rb, on parallel axes, each frame's origin on the belt's centre plane
; (as a catalog pulley's is). The spans pull the pulleys together; Gates'
; belt pull is their vector sum (VI).
(fn belt-drive (ra rb t0 more ...)
  (lambda (fa fb sa sb)
    (relation belt-drive (couple (open-belt-lines fa fb ra rb) t0 more) (data pitch-radii (vec ra rb)))))
; A chain drive: sprockets of pitch radii ra and rb (p / (2 sin (180 /
; N)), p the chain pitch, N teeth). The tight span pulls with the torque
; over the pitch radius, the slack span not at all (iwis' chain pull,
; centrifugal pull ignored).
(fn chain-drive (ra rb more ...)
  (lambda (fa fb sa sb)
    (relation chain-drive (couple (open-belt-lines fa fb ra rb) 0 more) (data pitch-radii (vec ra rb)))))
; One span of a belt at t0: span (vec p q) runs from p, where it leaves A,
; to q, where it meets B (the relating assembly's coordinates). A belt
; round several pulleys and idlers is a relation per span; the parts'
; balance tightens them together, keeping its length. Give each span the
; whole drive's efficiency.
(fn belt-span (span t0 more ...)
  (lambda (fa fb sa sb)
    (relation belt-span (couple (vec (vec (nth span 1) (unit (- (nth span 0) (nth span 1))) 0 (norm (- (nth span 0) (nth span 1))))) t0 more))))
; A screw A driving a nut B along A's z: the nut advances lead (mm) per
; turn of A, signed as the nut's slide reads A's turn (a right-hand thread
; turned about +z draws the nut toward -z). Thrust F comes with torque F
; lead / (2 pi), plus the efficiency's loss: drive torque F lead / (2 pi
; eta), back-driven F lead eta' / (2 pi) (Rexroth, screw assemblies, eq. 9
; and 10). eta' defaults to 2 - 1 / eta, the thread's friction to first
; order. At or below 0 the screw self-locks (Ewellix's "not
; backdriving"): lowering a load, the motor drives it down.
(fn lead-screw (lead more ...)
  (lambda (fa fb sa sb)
    (let ((o (get fa "origin")) (axis (get fa "z")) (h (/ lead (* 2 pi)))
          (eff (if (> (len more) 1) more (if (> (len more) 0) (vec (nth more 0) (- 2 (/ 1 (nth more 0)))) more))))
      (relation lead-screw (couple (vec (vec o axis h) (vec o (neg axis) h)) 0 eff) (data lead lead)))))

; ISO metric external screw thread (60 degree profile) along +z, unbounded:
; intersect with a slab for length. Major diameter d, pitch p (M6: d 6 p
; 1). Crest flat p/8 at d, root flat p/4 at the minor diameter d - 1.0825
; p. The field is distance-like, measured in the axial plane (the helix's
; lead angle is ignored). So offsets give clearances, provably: a tapped
; (internal) thread is the void (- (thread-iso d p) clearance).
(fn thread-iso (d p)
  ; Rotating (cos, sin) pairs puts atan2's one jump at the thread roots.
  (let ((r (len2 x y))
        (depth (* 0.5412659 p))
        (r-major (* 0.5 d))
        (r-minor (- (* 0.5 d) (* 0.5412659 p)))
        (k (/ (* 2 pi) p))
        (cx (/ x r))
        (sx (/ y r))
        (cz (cos (* k z)))
        (sz (sin (* k z)))
        ; Distance along z from the nearest crest centre, 0 .. p/2.
        (w (/ (abs (atan2 (- (* cz sx) (* sz cx)) (neg (+ (* cz cx) (* sz sx))))) k))
        ; The field's Lipschitz bound (its greatest slope). A pitch too
        ; coarse for the diameter (r-flank <= 0) makes the field jump at the
        ; axis, so it has none (0 promises nothing).
        (r-flank (- r-minor (* 0.65 p)))
        (slope (select (neg (min r-flank p))
                 (sqrt (+ 1 (* 0.75 (square (/ p (* 2 pi r-flank))))))
                 0)))
    (lipschitz slope
      (min (- r r-minor)
           (max (- r r-major)
                (+ (* 0.5 (- r r-minor)) (* 0.8660254 (- w (* 0.375 p)))))))))

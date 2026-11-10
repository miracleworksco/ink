; Two nameplates lettered in std's sans: one with its text engraved into
; the top face, one with it standing proud. Each plate is sized from its
; text's box (lettering-box): change a string or the letter height and the
; plate resizes.

; The two nameplates side by side.

(import "std")

(param letter-h 8)       ; cap height
(param margin 3)         ; letters' box to plate edge
(param plate-t 3)
(param depth 0.6)        ; engraving depth
(param relief 0.8)       ; embossing height

(def face (font "sans" (wght 500)))
(def engraved-text "PART No. 042-B")
(def embossed-text "Hello, World! 25°C")

; A plate round its text's box, the box's centre on the origin.
(fn plate-2d (L)
  (let ((b (lettering-box L)))
    (rounded-rect-2d (+ (- (nth b 2) (nth b 0)) (* 2 margin)) (+ (- (nth b 3) (nth b 1)) (* 2 margin)) 2)))
(fn label-2d (L)
  (let ((b (lettering-box L)))
    (move (lettering-2d L) (* -0.5 (+ (nth b 0) (nth b 2))) (* -0.5 (+ (nth b 1) (nth b 3))) 0)))

(def layout-a (lettering engraved-text face letter-h))
(def layout-b (lettering embossed-text face letter-h))
(def plate-a (extrude-z (plate-2d layout-a) 0 plate-t))
(def text-a (label-2d layout-a))
(def text-b (label-2d layout-b))

(body engraved
  (shape (move (subtract plate-a (extrude-z text-a (- plate-t depth) (+ plate-t 1))) 0 15 0))
  (attr color "#B8C2CC"))
(body embossed
  (shape (move (min (extrude-z (plate-2d layout-b) 0 plate-t)
                    (extrude-z text-b (- plate-t 0.5) (+ plate-t relief)))
               0 -15 0))
  (attr color "#C8A060"))

; The claims read the bodies in their own frame (moved back from the layout).
(def a (move (shape-of engraved) 0 -15 0))
(def b (move (shape-of embossed) 0 15 0))
(fn slab (z0 z1 g) (extrude-z g z0 z1))

; The text keeps the margin from each plate's edge, less 0.1.
(assert engraved-text-on-face
  (subset (slab 0.2 0.8 text-a) (slab 0 1 (- (plate-2d layout-a) (- 0.1 margin)))))
(assert embossed-text-on-face
  (subset (slab 0.2 0.8 text-b) (slab 0 1 (- (plate-2d layout-b) (- 0.1 margin)))))

; Engraving: every letter is cut to the full depth, to 0.05: nothing of
; the plate lies within the letters shrunk by 0.05, at any height above
; the groove floor. Under the groove the plate is solid.
(fn groove (g) (slab (+ (- plate-t depth) 0.05) plate-t g))
(assert engraving-cut (empty (max a (groove (- text-a -0.05)))))
(assert engraving-floor
  (subset (slab 0.1 (- plate-t depth 0.1) (- (plate-2d layout-a) -0.1)) a))

; Embossing: every letter, in full to 0.05, stands `relief` proud of the
; plate, and nothing else does.
(assert embossed-letters
  (subset (slab plate-t (- (+ plate-t relief) 0.05) (- text-b -0.05)) b))
(assert embossed-only-letters
  (empty (max b (slab (+ plate-t 0.05) (+ plate-t relief 1) (neg (- text-b 0.05))))))

; Two nameplates from std's single-stroke text (text-2d): one with its text
; engraved into the top face, one with it standing proud. It shows a part
; sized from its text (text-width): change a string or the letter height
; and the plate resizes.

; The two nameplates side by side.

(import "std")

(param letter-h 8)       ; cap height
(param stroke 0.8)       ; stroke width: the engraving cutter, or a printed line
(param margin 3)         ; text to plate edge
(param plate-t 3)
(param depth 0.6)        ; engraving depth
(param relief 0.8)       ; embossing height

(def engraved-text "PART No. 042-B")
(def embossed-text "Hello, World! 25°C")

; A plate sized to its text: the line's advance width (letters and the gaps
; the font puts round them) plus the margin, and the cap height plus room
; for descenders (0.35 of the cap height in this font) and the strokes' ends.
(fn plate-2d (s)
  (rounded-rect-2d (+ (text-width s letter-h stroke 0) (* 2 margin))
                   (+ (* 1.35 letter-h) stroke (* 2 margin))
                   2))
; The text, centred: its cap band on y = 0.175 h so that, with descenders,
; it sits centred on the plate.
(fn label-2d (s) (move (text-2d-styled s letter-h stroke 0 0.5 0.5) 0 (* 0.175 letter-h) 0))

(def plate-a (extrude-z (plate-2d engraved-text) 0 plate-t))
(def text-a (label-2d engraved-text))
(def text-b (label-2d embossed-text))

(body engraved
  (shape (move (subtract plate-a (extrude-z text-a (- plate-t depth) (+ plate-t 1))) 0 15 0))
  (attr color "#B8C2CC"))
(body embossed
  (shape (move (min (extrude-z (plate-2d embossed-text) 0 plate-t)
                    (extrude-z text-b (- plate-t 0.5) (+ plate-t relief)))
               0 -15 0))
  (attr color "#C8A060"))

; The claims read the bodies in their own frame (moved back from the layout).
(def a (move (shape-of engraved) 0 -15 0))
(def b (move (shape-of embossed) 0 15 0))
(fn slab (z0 z1 g) (extrude-z g z0 z1))

; The text keeps the margin from each plate's edge, less 0.1 (at the top
; it is exactly the margin: the plate is sized to the caps' strokes).
(assert engraved-text-on-face
  (subset (slab 0.2 0.8 text-a) (slab 0 1 (- (plate-2d engraved-text) (- 0.1 margin)))))
(assert embossed-text-on-face
  (subset (slab 0.2 0.8 text-b) (slab 0 1 (- (plate-2d embossed-text) (- 0.1 margin)))))

; Engraving. Every stroke is cut to the full depth and width, to 0.05:
; nothing of the plate lies within the text shrunk by 0.05, at any height
; above the groove floor. Each stroke is a segment swept by a disk of
; diameter `stroke`, so this also proves the groove at least stroke - 0.1
; wide all along (a disk of that diameter, clear, round every point of
; every stroke's centre line). Under the groove the plate is solid.
(fn groove (g) (slab (+ (- plate-t depth) 0.05) plate-t g))
(assert engraving-cut (empty (max a (groove (- text-a -0.05)))))
(assert engraving-floor
  (subset (slab 0.1 (- plate-t depth 0.1) (- (plate-2d engraved-text) -0.1)) a))

; Embossing: every stroke, in full to 0.05 (so at least stroke - 0.1 wide,
; as above), stands `relief` proud of the plate, and nothing else does.
(assert embossed-letters
  (subset (slab plate-t (- (+ plate-t relief) 0.05) (- text-b -0.05)) b))
(assert embossed-only-letters
  (empty (max b (slab (+ plate-t 0.05) (+ plate-t relief 1) (neg (- text-b 0.05))))))

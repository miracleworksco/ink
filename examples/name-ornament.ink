; A name ornament made to order: the customer's name, raised in a joined
; script on a disc with a hole to hang it by. The name is a text param
; (`ink assert name-ornament.ink --param name=Ada`); the letters are
; sized to fill the space the disc leaves them.
;
; Its claims hold for this name (the letters fit the disc's face, clear
; the hole, and print: no stroke under 0.8 mm), and one for every name the
; ornament allows (up to 8 letters, a hyphen, an apostrophe, a space): its
; letters are never under 3 mm tall.

(import "std")

(def letters "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz -'")
(param name "Grace" (chars letters) (len 1 8))

(param r 40)            ; the disc's radius
(param thick 3)         ; its thickness
(param relief 1.5)      ; how far the letters stand proud
(def font-name (font "script"))
(def text-w 64)         ; the box the letters fill
(def text-h 26)
(def hole-r 2.5)
(def hole-y (- r 6))

; The cap height that fits a name whose letters, at cap height 1, span a
; box w by h: at most 14 mm.
(fn fit (w h) (min 14 (/ text-w w) (/ text-h h)))
(def at-1 (lettering-box (lettering name font-name 1)))
(def cap (fit (- (nth at-1 2) (nth at-1 0)) (- (nth at-1 3) (nth at-1 1))))
(def L (lettering name font-name cap))
(def box (lettering-box L))
; The letters, their box centred on the disc.
(def letters-2d (move (lettering-2d L) (* -0.5 (+ (nth box 0) (nth box 2))) (* -0.5 (+ (nth box 1) (nth box 3))) 0))

(def disc (subtract (extrude-z (circle-2d r) 0 thick) (move (cylinder-z hole-r (* 3 thick)) 0 hole-y 0)))
(body ornament
  (shape (min disc (extrude-z-round letters-2d (- thick 0.5) (+ thick relief) 0.4)))
  (attr material "pla")
  (attr color "#B03A2E"))

; This name: its letters lie on the disc's face, 3 mm in from its rim, and
; clear the hanging hole by 2 mm (claims of 2D profiles: asked in the plane).
(assert on-face (subset letters-2d (circle-2d (- r 3))))
(assert clear-of-hole (clearance letters-2d (move (circle-2d hole-r) 0 hole-y 0) 2))

; Every name allowed: the widest and tallest any 8 letters can be (a sound
; bound over every glyph the font has, its swashes and joining forms
; among them, and all its kerning) still fit at 3 mm.
(def worst (lettering-max letters 8 font-name 1 0))
(def smallest-cap (fit (nth worst 0) (- (nth worst 2) (nth worst 1))))
(assert legible-for-any-name (value smallest-cap >= 3))
; At that size not every letter of this font prints: a specimen of them
; all, (lettering-styled letters font-name smallest-cap smallest-cap 1
; 0.5), has a Q whose tail is 0.34 mm thick. So each order is checked at
; its own name and size: Ada's last a, the form the font ends a word with,
; tapers to 0.78 mm, and that order fails.
(assert strokes (min-thickness letters-2d > 0.8))

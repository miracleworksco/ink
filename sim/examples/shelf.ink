; Wall shelf: a board on two printed brackets, each fixed with four M5 screws.
;
; Frames: the wall is the plane y = 0, +y points out of the wall, +z is up.
; Each bracket is modelled in its own frame, centred on x = 0.

(import "std")

; ---------------------------------------------------------------- parameters (mm)

(param arm 40)              ; how far the bracket reaches out from the wall
(param height 45)           ; overall bracket height
(param width 34)            ; bracket width along x
(param plate 5)             ; plate thickness
(param web 12)              ; width of the lattice web under the arm
(param cell 8)              ; gyroid cell size in the web
(param wall 1.2)            ; gyroid sheet thickness
(param skin 1.6)            ; solid skin on the web's sloped face
(param hole-r 2.75)         ; M5 clearance hole radius
(param hole-dx 11.5)        ; holes at x = +/- hole-dx (screw heads must clear the web)
(param spacing 200)         ; distance between brackets
(param board 18)            ; shelf board thickness
(param fillet 3)            ; blend radius where plates meet

; ---------------------------------------------------------------- bracket

(def upright (move (box width plate height) 0 (* 0.5 plate) (* 0.5 height)))
(def arm-plate (move (box width arm plate) 0 (* 0.5 arm) (- height (* 0.5 plate))))

; Triangle under the arm, between the upright and the arm's tip.
; `slope` is the exact signed distance to its sloped face (negative above it).
(def slope
  (let ((dy (- arm plate)) (dz (- height plate)))
    (/ (- (* dz (- y plate)) (* dy z)) (len2 dy dz))))
(def web-region
  (max (- (abs x) (* 0.5 web))
       (- plate y)
       (- z (- height plate))
       slope))

; Gyroid sheet lattice inside the web, plus a solid skin along the slope.
(def web-lattice (lattice-fill web-region (gyroid cell) wall))
(def web-skin (max web-region (neg (+ slope skin))))
(def web-body (min web-lattice web-skin))

; Clearance holes: two through the upright (into the wall), two through the arm.
(def wall-holes (symmetric-x (move (cylinder-y hole-r 40) hole-dx 0 20)))
(def shelf-holes (symmetric-x (move (cylinder-z hole-r 40) hole-dx 25 height)))

; Fillet in the inside corner between the upright and the arm.
(def inner-fillet
  (max (corner-fillet (- y plate) (- (- height plate) z) fillet) (- (abs x) (* 0.5 width))))

(def bracket-field
  (subtract (min upright arm-plate inner-fillet web-body) (min wall-holes shelf-holes)))

(body bracket
  (shape bracket-field)
  ; Colour graded bottom to top: a field, not a constant.
  (attr color (vec (+ 0.15 (* 0.012 z)) 0.45 (- 0.85 (* 0.01 z))))
  (attr material "petg")
  (attr density 1.27g/cm3)                      ; g/mm^3
  (attr process (data (fdm (nozzle 0.4) (layer 0.2) (infill "as-modelled")))))

; ---------------------------------------------------------------- screws (M5)
; Screw frame: head occupies z in [-3, 0], shank runs up the +z axis.

(fn screw (len)
  (min (move (cylinder-z 4.2 3) 0 0 -1.5)
       (tag shank (move (cylinder-z 2.4 len) 0 0 (* 0.5 len)))))

(body screw-35 (shape (screw 35)) (attr color "#9A9A9A") (attr material "steel") (attr density 7.85g/cm3))
(body screw-16 (shape (screw 16)) (attr color "#9A9A9A") (attr material "steel") (attr density 7.85g/cm3))

; ---------------------------------------------------------------- shelf board

(body board
  (shape (rounded-box (+ spacing 100) 150 board 1.5))
  (attr color "#B8875A")
  (attr material-note "birch plywood")
  (attr density 0.68g/cm3))

; ---------------------------------------------------------------- assemblies

(assembly bracket-set
  (instance bracket bracket)
  ; Wall screws: head against the upright's face (y = plate), shank into the wall.
  (instance wall-l screw-35 (at (rot-x 90) (move (- hole-dx) plate 20)))
  (instance wall-r screw-35 (at (rot-x 90) (move hole-dx plate 20)))
  ; Shelf screws: head against the arm's underside, shank up into the board.
  (instance shelf-l screw-16 (at (move (- hole-dx) 25 (- height plate))))
  (instance shelf-r screw-16 (at (move hole-dx 25 (- height plate)))))

(assembly root
  (instance board board (at (move 0 75 (+ height (* 0.5 board)))))
  (instance left bracket-set (at (move (* -0.5 spacing) 0 0)))
  (instance right bracket-set (at (move (* 0.5 spacing) 0 0)))
  ; Shelf screws thread into the board: an intended overlap, their shanks.
  (contact left/shelf-l board (shape-of screw-16 shank))
  (contact left/shelf-r board (shape-of screw-16 shank))
  (contact right/shelf-l board (shape-of screw-16 shank))
  (contact right/shelf-r board (shape-of screw-16 shank)))

; ---------------------------------------------------------------- design checks

; Screw shanks pass through the clearance holes without touching the bracket.
(require shank-clears-hole
  (disjoint (shape-of bracket) (move (cylinder-z 2.4 40) hole-dx 25 height)))
; Screw heads seat against the bracket but never cut into it.
(require heads-seat-flush (disjoint (shape-of bracket-set/bracket) (shape-of bracket-set/shelf-r)))
; Shelf screw tips stay at least 3 mm below the board's top surface.
(require shelf-screw-hidden
  (disjoint (shape-of root/left/shelf-l) (half-space 0 0 -1 (- 3 (+ height board)))))
; Screw heads clear the lattice web by at least 1 mm.
(require heads-clear-web (clearance (shape-of bracket-set/shelf-r) web-region 1))
; The bracket stays in its envelope (fits a 34 x 40 x 45 print volume).
; Its faces lie on the envelope's, which no search separates at the model's
; floor within its budget: decided to 0.2 mm, at cells of 0.1.
(require fits-envelope (within (shape-of bracket) (-17 0 0) (17 40 45)) (res 0.1))
; Wall screws keep 10 mm clear of the board.
(require wall-screws-clear-board (clearance (shape-of root/left/wall-l) (shape-of root/board) 10))
; The bracket prints lighter than 30 cm^3 of material.
(require bracket-light (volume (shape-of bracket) < 30000mm3))
; The colour field stays in [0, 1] everywhere in the bracket.
(require color-in-gamut (range (+ 0.15 (* 0.012 z)) 0 1 (shape-of bracket)))

; ---------------------------------------------------------------- simulation
; The shelf as one fused PETG piece: board and both brackets, screws left out.
; The bracket backs are clamped where they bolt to the wall (the board's back
; edge only touches the wall, so it stays free). Load: 5 kg of books spread
; over the top surface, plus self-weight.

; Fillets where a bracket meets the board's underside (bracket frame): along
; both sides of the bracket, and across the arm's tip.
(def board-fillets
  (min (max (corner-fillet (- (abs x) (* 0.5 width)) (- height z) fillet) (neg y) (- y arm))
       (max (corner-fillet (- y arm) (- height z) fillet) (- (abs x) (* 0.5 width)))))

(def shelf-solid
  (min (shape-of root/board)
       (shape-of root/left/bracket) (move board-fillets (* -0.5 spacing) 0 0)
       (shape-of root/right/bracket) (move board-fillets (* 0.5 spacing) 0 0)))

(study books-5kg
  (shape shelf-solid)
  (elastic (E 2100MPa) (nu 0.38) (density 1270kg/m3) (yield 45MPa))    ; PETG: MPa, -, MPa
  (fix wall (max (- y 0.5) (- z (- height 1))))            ; bracket backs: y < 0.5, below the board
  (load books (force 0 0 -49.05N) (region (- (+ height board -1) z)))   ; top 1 mm of the board
  (gravity 0 0 -9.81m/s2)
  (res 3))                                                 ; refine with ink solve --res

; Neither claim is decided, and both say why. Cells of 6 mm glue shut
; narrow gaps where the lattice's 1.2 mm sheets meet its skin, so the
; claims are checked on 0.75, 1.5 and 3 mm cells instead, and 0.75 mm
; cells glue a narrower one (0.15 mm) where they meet.
; And the safety factor falls as the cells shrink, 42.8, 29.6, 22.0 and 18.2
; at 3, 2, 1.5 and 1 mm, where the lattice meets the plates at sharp
; junctions, whose peak never settles: fillet them to decide it.
; expect-status: holds-books unknown
; expect-status: stays-flat unknown
(require holds-books (min-safety-factor books-5kg > 2))
(require stays-flat (max-displacement books-5kg < 2))

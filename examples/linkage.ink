; A four-bar linkage, its crank turning once per unit time. The coupler and
; rocker angles are solved at every moment so the loop closes.
(import "std")
(t-span 0 1)

(param crank 20)      ; ground pivot A at the origin, crank A-B
(param coupler 60)    ; coupler B-C
(param rocker 40)     ; rocker D-C, ground pivot D
(param ground 50)     ; A to D along x
(def crank-angle (* 2 pi t))

(solve loop
  (unknowns (coupler-angle 0.8) (rocker-angle 1.6))
  (== (+ (* crank (cos crank-angle)) (* coupler (cos coupler-angle)))
      (+ ground (* rocker (cos rocker-angle))))
  (== (+ (* crank (sin crank-angle)) (* coupler (sin coupler-angle)))
      (* rocker (sin rocker-angle))))

; A bar of length l from the origin along +x, 6 wide, 3 thick. The coupler
; lies on the crank and the rocker.
(fn bar (l) (move (rounded-rect-2d-bar l) 0 0 0))
(fn rounded-rect-2d-bar (l) (extrude-z (move (rounded-rect-2d (+ l 6) 6 3) (* 0.5 l) 0 0) -1.5 1.5))

(body crank-link (shape (bar crank)) (attr color "#D79B00"))
(body coupler-link (shape (bar coupler)) (attr color "#6C8EBF"))
(body rocker-link (shape (bar rocker)) (attr color "#82B366"))

(assembly root
  (instance crank crank-link (at (rot-z (* (/ 180 pi) crank-angle))))
  (instance coupler coupler-link
    (at (rot-z (* (/ 180 pi) coupler-angle)) (move (* crank (cos crank-angle)) (* crank (sin crank-angle)) 3)))
  (instance rocker rocker-link (at (rot-z (* (/ 180 pi) rocker-angle)) (move ground 0 0))))

; The coupler's far end meets the rocker's tip, closing the loop over the
; whole turn: the solved angles are enclosed over each piece of it.
(def coupler-end (vec (+ (* crank (cos crank-angle)) (* coupler (cos coupler-angle)))
                      (+ (* crank (sin crank-angle)) (* coupler (sin coupler-angle)))))
(def rocker-tip (vec (+ ground (* rocker (cos rocker-angle))) (* rocker (sin rocker-angle))))
(assert loop-closes (value (distance coupler-end rocker-tip) < 0.001))

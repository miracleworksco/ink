; Validation P3's four-bar with the coupler placed on the crank's pin, so
; the solve's unknown th3 is its angle relative to the crank and turns a
; whole turn while the crank does. Each moment's root is followed along t
; from the load-time root, so the rocker stays on the open assembly (87 to
; 126 deg). Solved afresh from the load-time root at each step, it landed
; on the crossed assembly at some (th4 = 234.8 deg at t = 0.245) and the
; claim failed.
(import "std")
(t-span 0 0.6283185307)
(body crank (shape (move (box 100 4 4) 50 0 0)) (frame pin (at (move 100 0 0))))
(body coupler (shape (move (box 400 4 4) 200 0 0)) (frame pin (at (move 400 0 0))))
(body rocker (shape (move (box 300 4 4) 150 0 0)) (frame pin (at (move 300 0 0))))
(param th3 0.85) (param th4 1.52)
(assembly (fourbar (p 0))
  (instance crank crank (at) (turn z p))
  (instance coupler coupler (at (rot-z (rad->deg th3)) (frame-of crank pin)))
  (instance rocker rocker (at (rot-z (rad->deg th4)) (move 350 0 0))))
(assembly root (instance f (fourbar (* (rad->deg 10) t))))
(def gap (- (get (frame-of root/f/coupler pin) "origin") (get (frame-of root/f/rocker pin) "origin")))
(solve loop (unknowns th3 th4) (== (nth gap 0) 0) (== (nth gap 1) 0))
(assert open (value (sin th4) > 0.5))
; expect-status: open holds

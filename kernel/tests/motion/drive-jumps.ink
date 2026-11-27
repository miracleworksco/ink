; Claims written on a def's diff. A def no joint has is not checked for
; jumps, its units not known: where it may jump (a select or floor
; switching on t), the claim is unchecked there and unknown. So a def
; stepping by 30 is unknown at the step, through a diff of its diff and a
; warp of t too, and so is a select chain whose branches meet in reals
; (t1 = 1/3, its constants folded to floats a few ulps apart), its folded
; error amplified, or not. None of them holds or fails on the jump.
(import "std")
(t-span 0 0.75s)
(def e (select (- t 0.5s) (* 10 (/ t 1s)) (+ 30 (* 10 (/ t 1s)))))
(require speed (value (diff e t) < 100))
(require accel (value (diff (diff e t) t) < 1000))
(require warped (value (warp (diff e t) x y z (* 2 (/ t 1s))) < 100))
(def f (select (- (/ t 1s) 0.5) (* 4 (/ t 1s)) (+ 3 (* 4 (/ t 1s)))))
(require energy (value (* 0.5 2 (square (diff f t))) < 100) (accuracy 5))
(def d (select (- (/ t 1s) 0.5) 0 (* 100 (- (/ t 1s) 0.5))))
(require kick (value (diff (diff d t) t) < 1000))
(def acc 1500mm/s2)
(def v 500mm/s)
(def t1 (/ v acc))
(def p1 (* 0.5 acc t1 t1))
(def belt (select (- t t1) (* 0.5 acc t t) (+ p1 (* v (- t t1)))))
(require belt-speed (value (diff belt t) < 510mm/s))
; The 30 step with a branch never taken whose bound overflows.
(def a (+ (select (- (/ t 1s) 0.5) 0 30) (select (- (/ t 1s) 2) (* 10 (/ t 1s)) (exp (* 1000 (/ t 1s))))))
(require inf-branch (value (diff a t) < 20))
; ln t reads -inf over the first step, the step in it.
(def g (+ (select (- (/ t 1s) 0.03) 0 30) (select (- (/ t 1s) 0.01) 0 (* 0.001 (ln (/ t 1s))))))
(require ln-start (value (diff g t) < 20))
; A cap no value reaches.
(def b (min (+ (* 10 (/ t 1s)) (select (- (/ t 1s) 0.5) 0 30)) 1e20))
(require big-cap (value (diff b t) < 20))
; A huge value cancelling to a small one, and a 0.25 step.
(def c (+ (select (- (/ t 1s) 0.5) 0 0.25) (* 1e-15 (- (+ (/ t 1s) 1e15) 1e15))))
(require cancel (value (diff c t) < 1))
; Two foldings of 1/3 meeting at the switch, divided by 1e-6.
(def u1 (/ 1 3))
(def u2 (- 1 (/ 2 3)))
(def h (select (- (/ t 1s) u1) (/ (- (/ t 1s) u1) 1e-6) (/ (- (/ t 1s) u2) 1e-6)))
(require amplified (value (diff h t) < 2e6))
; expect-status: speed unknown
; expect-status: accel unknown
; expect-status: warped unknown
; expect-status: energy unknown
; expect-status: kick unknown
; expect-status: belt-speed unknown
; expect-status: inf-branch unknown
; expect-status: ln-start unknown
; expect-status: big-cap unknown
; expect-status: cancel unknown
; expect-status: amplified unknown

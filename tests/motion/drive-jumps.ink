; Claims written on a def's diff. A def no joint has is not checked for
; jumps, its units not known: where it may jump (a select or floor
; switching on t), the claim is unchecked there and unknown. So a def
; stepping by 30 is unknown at the step, through a diff of its diff and a
; warp of t too, and so is a select chain whose branches meet in reals
; (t1 = 1/3, its constants folded to floats a few ulps apart), its folded
; error amplified, or not. None of them holds or fails on the jump.
(import "std")
(t-span 0 0.75)
(def e (select (- t 0.5) (* 10 t) (+ 30 (* 10 t))))
(assert speed (value (diff e t) < 100))
(assert accel (value (diff (diff e t) t) < 1000))
(assert warped (value (warp (diff e t) x y z (* 2 t)) < 100))
(def f (select (- t 0.5) (* 4 t) (+ 3 (* 4 t))))
(assert energy (value (* 0.5 2 (square (diff f t))) < 100) (accuracy 5))
(def d (select (- t 0.5) 0 (* 100 (- t 0.5))))
(assert kick (value (diff (diff d t) t) < 1000))
(def acc 1500)
(def v 500)
(def t1 (/ v acc))
(def p1 (* 0.5 acc t1 t1))
(def belt (select (- t t1) (* 0.5 acc t t) (+ p1 (* v (- t t1)))))
(assert belt-speed (value (diff belt t) < 510))
; The 30 step with a branch never taken whose bound overflows.
(def a (+ (select (- t 0.5) 0 30) (select (- t 2) (* 10 t) (exp (* 1000 t)))))
(assert inf-branch (value (diff a t) < 20))
; ln t reads -inf over the first step, the step in it.
(def g (+ (select (- t 0.03) 0 30) (select (- t 0.01) 0 (* 0.001 (ln t)))))
(assert ln-start (value (diff g t) < 20))
; A cap no value reaches.
(def b (min (+ (* 10 t) (select (- t 0.5) 0 30)) 1e20))
(assert big-cap (value (diff b t) < 20))
; A huge value cancelling to a small one, and a 0.25 step.
(def c (+ (select (- t 0.5) 0 0.25) (* 1e-15 (- (+ t 1e15) 1e15))))
(assert cancel (value (diff c t) < 1))
; Two foldings of 1/3 meeting at the switch, divided by 1e-6.
(def u1 (/ 1 3))
(def u2 (- 1 (/ 2 3)))
(def h (select (- t u1) (/ (- t u1) 1e-6) (/ (- t u2) 1e-6)))
(assert amplified (value (diff h t) < 2e6))
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

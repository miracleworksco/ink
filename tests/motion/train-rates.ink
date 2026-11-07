; Lines turned at exact rates of one drive, as a roller floor's are
; written: a at 1/3 of it, b and c at 2/3. A pair of arms on two of them
; reads the drive only through their angles, so it repeats once both have
; turned whole turns: 1080 degrees of the drive, a second. b's arm meets
; a's only near 700 degrees of the drive in each 1080 (a at 233, b at 467),
; found at a moment there. c's arm sweeps across a's but never meets it
; (2.9 mm apart near 896): proven over one second, and so over the ten.
(import "std")
(t-span 0 10)
(def drive (* 1080 t))
(body arm (shape (box 30 2 2)))
(assembly line (instance arm arm (at (move 25 0 0))))
(assembly root
  (instance a line (at) (turn z (/ drive 3)))
  (instance b line (at (move 0 -35 0)) (turn z (/ (* 2 drive) 3)))
  (instance c line (at (move 45 0 0)) (turn z (/ (* 2 drive) 3))))
(assert meets (disjoint (shape-of root/a/arm) (shape-of root/b/arm)))
(assert misses (disjoint (shape-of root/a/arm) (shape-of root/c/arm)))
; expect-status: meets fails
; expect-status: misses holds

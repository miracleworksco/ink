# Changes

## 0.2.0

### Upgrading a 0.1 file

Each line is a change that stops a 0.1 file loading or changes what an
answer means, the message 0.2 gives, and the fix. The messages are quoted
as `ink require FILE --plain` prints them.

**`assert` is `require`.** The form is the same; only the word changed.

    a.ink:3: `assert` is now `require`
      write (require big (volume (shape-of b) > 500)): the form is the same, only the word changed

**`ink assert` is `ink require`.** One command asks every requirement,
the ones a joint's limits and a relation make included. Separate
`ink-require` programs, where a build had them, are this command too.

    unknown command 'assert'
    `ink assert` is now `ink require`: the same file, with (assert ...) written (require ...)

**A number that is not a length says its unit.** A bare number is a
length in mm (an area, a volume), a pure number, or an angle in degrees
where a pose is taken. A mass, a force, a torque, a time, a speed, a
pressure, a temperature and the rest say theirs: `(attr mass 12g)`,
`(gravity 0 0 -9.81m/s2)`, `(t-span 0 2s)`, `(speed-limits -10mm/s
10mm/s)`, `(elastic (E 200GPa))`. A bare 0 stands. A value worked out
says it by a product: `(* 1g (get spec "mass-g"))`.

    a.ink:3: attr mass takes a mass, and 12 is a bare number: say its unit (12g, if it is in g)
      a bare number is a length in mm or a pure number; a mass says its unit: g mg kg lb oz
    a.ink:2: gravity takes an acceleration, and -9.81 is a bare number: say its unit (-9.81m/s2, if it is in m/s2)
    a.ink:2: t-span takes a time, and 2 is a bare number: say its unit (2s, if it is in s)
    a.ink:3: speed-limits takes a speed, and -10 is a bare number: say its unit (-10mm/s, if it is in mm/s)
    a.ink:3: E takes a pressure or stress, and 1000000 is a bare number: say its unit (1000000MPa, if it is in MPa)

**What the model works out says its unit.** `(reaction PATH)` gives
forces and torques, `(tension A B)` forces and
`(effort-rms PATH)` a torque or a force, where 0.1 gave bare numbers a
design typed the unit of (`(* 1N (nth r 2))`). That product is now a
force times a force; drop it. What a model number is compared with says
its unit, and so does a load handed to std's bearing, carriage and belt
helpers, which read their records' ratings in newtons themselves.

    a.ink:13: '-' takes values of one quantity, but a torque worked out at a.ink:13 is a torque and 330 is a bare number: say its unit (330N*mm, if it is in N*mm)

A bare number beside such a quantity in a sum is refused the same way
(`(- t 300)` beside a time, `(+ p 0.2)` beside a pressure), and so is a
bare value given to a param whose default is one (`--param m=50` where
`(param m 5g)`: write `m=50g`).

    a.ink:3: '-' takes values of one quantity, but a time worked out at a.ink:3 is a time and 300 is a bare number: say its unit (300s, if it is in s)
    param m is written 5g; give m=50 a unit

The message offers the unit 0.1 read the number in. If the number was
meant in another, write that one: 0.1 read it wrong.

**`t` is seconds, and a motion is written with units.** A joint's drive
is an angle (a turn) or a length (a slide). A number times `t` is a time,
so say the rate: `(* 360deg/s t)`, `(* 20mm/s t)`. A time is taken from
`t` with its unit, `(- t 0.3s)`; a share of a time is a pure number, `(/ t
0.5s)`. A formula that works in bare numbers takes `t` out as its number
of seconds: `(/ t 1s)`.

    a.ink:4: a turn's drive takes an angle (deg), and (* 360 t) is a time (its t is written at a.ink:4)
      t is in seconds, so a number times t is a time: say the rate's unit, (* 360deg/s t), or the travel's, (* 90deg (/ t 0.5s))

**A sine takes an angle, not a time.** `(sin (* pi t))` was radians a
second: write `(sin (* 180deg/s t))`.

    a.ink:4: 'sin' takes an angle, and (* pi t) is a time: t is in seconds, so say how fast the angle turns
      write the rate with its unit, (sin (* 180deg/s t)) or (sin (* 2rad/s t)); a bare number here was radians a second

**Angles are `deg` and `rad`.** `30deg` and `1.2rad` are numbers with a
unit, which 0.1 refused. A sine converts one: `(sin 30deg)` is 0.5. A
bare number is read as it was: degrees where a pose is taken (`rotate-z`,
a turn), radians in `sin`, `cos` and `tan`, so `(sin 30)` is still -0.988.
The warning that a sine's number "looks like degrees" is gone. A turning
rate is `deg/s`, `rpm` or `rad/s`, a turning acceleration `deg/s2` or
`rad/s2`, and `std`'s profile takes them: `(trapezoid-position t
9000deg/s2 1800deg/s 4500deg)`.

**`(solve ...)`, `(optimize ...)` and `ink optimize` are gone**, until
the constraint solver and optimiser return as a package. Give the value
as a `param` and state what it must meet as a `require`. (`ink solve`,
which runs a study, stays.)

    a.ink:3: unknown top-level form 'solve'
      (solve ...) left in 0.2: the constraint solver and optimiser return as a package. Until then give the value as a param and state what it must meet as a (require ...)
    unknown command 'optimize'
    `ink optimize` left in 0.2 with (solve ...) and (optimize ...): the constraint solver and optimiser return as a package

**A yes says what it covers.** `ink require` exits 0 only where every
requirement is yes over everything it reads; 3 where one is no, unknown,
or yes over part of it; 1 on an error. A requirement that reads `t` in a
file with no `(t-span)` is asked at one moment and says so:

    clear: YES, proven, at one value of t  (a.ink:4)
      note: asked at t = 0 only: it reads t, and no (t-span) covers it

and the file is not all yes (exit 3). Add the span the motion runs over.
`(require NAME (within-ratings PATH))` says the parts under PATH it found
no rating for (`NOT COVERED`) and what it assumed (`ASSUMING`), and
every answer over a design's parts says first what it checked of where
they sit (`PLACEMENT`: no clash at the start, and what `ink interfere`
and a picture are still owed): below.

**A yes over a motion covers what lies between the moments read.** An
answer read at moments of the motion (`YES, estimated, at 17 values of
t`) was a yes wherever every moment passed, whatever happened between two
of them. Now each step between two moments is asked too, with `t` ranged
over that step alone; a clash found there is a no, and a step nothing
shows is not settled, saying where (exit 3). The same holds of
`ink interfere`'s moving pairs, of a study whose load the model works
out (read at the model's own steps, which follow how fast its loads
change, not at even ones), and of an estimate's condition
(`estimate-when`), which is shown over every moment its value is read at.
A file whose parts pass through each other between two moments, or whose
load swings faster than it was sampled, changes its answer; an
`overlaps` that holds through a motion (a pin that stays in its hole) is
now proven, by a point that stays in both. A `touches` over a motion
whose shapes move is not settled: claim a clearance or an overlap.

    clear: NO, proven, at every value of t  (a.ink:7)

    on: not settled: a source could not decide  (a.ink:5)
      note: not settled between the moments read, over t 0 to 10: what it reads may move further within a step than the room the claim has at the moments. A yes at the moments says nothing of what lies between them there: ask a shorter span (--t-span A:B), or slow what changes. A contact found at one moment says nothing of the next, so no step between them is shown: claim what can hold at every moment (a clearance, an overlap, disjoint, a value), or ask it at rest with why the motion can't change it

**A polygon is built within an allowance of the one written.** Every
`polygon` is one construction, which may be off the profile written by
its allowance (a ten-millionth of its points' box). A claim decided
within that of it is not settled, never read the wrong way: a point
exactly on a drawn corner or edge, a fit closer than the allowance at
an accuracy finer than it. Claim a point a little inside, or at an
accuracy coarser than the allowance.

    corner-in: not settled: a source could not decide  (a.ink:3)
      note: it reads a shape built within an allowance of the one written (polygon allowance 0 0.000000 within 0.000004), a bound on how far the shape built may be from it. It neither holds at every value in that band nor fails at every one, or it reads the band through a solve, which carries none

**A drafted polygon whose section changes within its span is refused.**
Walls that lean in until they meet, an edge drafted until it vanishes,
or a round that shrinks to nothing within z0..z1 are refused, saying
where. Draft it over the span where its section keeps its edges, and
write what lies beyond as its own shape (a polygon over the next span,
or the closed end as a cap).

    b.ink:2: the polygon's edge from (10, 0) (point 2) to (10, 1) (point 3): its corners' rounds overlap, or its walls cross, at z = 2.8356409098088546

**The least time is held back.** 0.2 leaves it out until it decides
claims: a `(least-time ...)` value and `ink effort --least-time` are
refused, saying so.

    a.ink:3: (least-time ...) is not in this release of ink: the least time to run a path within its actuators' limits is held back until it is ready to rely on

**Flow is held back.** 0.2 leaves flow (lattice Boltzmann on the GPU,
with conjugate heat) out until it is ready to rely on: a `(flow ...)` form
and `ink flow` are refused, saying so.

    a.ink:3: (flow ...) is not in this release of ink: flow (lattice Boltzmann on the GPU, with conjugate heat) is held back until it is ready to rely on

### What is new

**A fix holds its region's material, and no more.** A fix held each cut
cell its region touched whole, so a hold reached up to a cell past where it
was written, and three grids, rounding its edge alike, agreed on the wrong
answer. Now a fix holds the material inside its region, traced at its edge
(Nitsche's method), and a part held there reads as written. A hold thinner
than two cells with free material on both sides, or one just inside a face
with a sliver of free material beyond it, is refused, saying what to write
(make it two cells thick, or extend it to the face). Some deflections and
temperatures near a hold change; a hold slab inside a bar is often not
settled where it was answered.

**A peak displacement or temperature is read on the part's face.** It was
read at the last node inside the part, up to a cell short of its surface,
and every grid short alike. A bar's end now reads its stretch, and a
heated face its temperature.

**Three grids are trusted only where they converge, or where nothing could
run away.** Differences no larger than what the solves resolve are noise,
not a trend. Three grids converging faster than their elements can are
luck: their grid error is three times their whole movement, not their last
difference. Three grids moving apart, or swinging wider, are given three
times their larger change, unless the solve finds something that could make
the value grow without bound as the cells shrink (a load or hold under a
cell across, a peak stress at a hold's edge or a sharp corner, a mode that
shrinks with the cells): then the claim is not settled, saying which. What
the solves resolve is taken of what a displacement or temperature varies by
over the part, not of an offset it carries whole: a small temperature rise
on a hot part, or a small deflection on a part a hold moves far, now asks
for finer cells where it is still growing. Some claims that held on three
grids now ask for finer cells.

**A held temperature may vary over its face and in time.** `(temperature
NAME REGION VALUE)` takes a value of x, y, z and params (a cylinder's side
at 80 cos(theta) degF), as a fix's move does, and in a transient of t: each
step holds it as it is at the step's end, and the steps end where it may
jump. A held temperature that reads t outside a transient is refused.

**A mass is proven, not estimated.** A body's mass, centre and inertia, and
an integral, mean or volume claim, are decided by their proven bounds, and
`ink mass` states them. An actuator's effort carries how far those bounds
move it, so an effort within them of its limit is not settled.

**A value claim is decided in real numbers.** A number as written is the
exact fraction its digits write (`0.1 + 0.2 <= 0.3` is yes; `1 ft / 12 in
> 1` is no), and so is a record's number and a param written in digits; an
arithmetic step that rounds is the interval its rounding cannot leave. A
claim that was decided by a float's rounding either way now reads as the
numbers do. A band through an estimate's condition is carried where the
condition holds over it.

**Ratings.** A part can carry what its maker rates it for (a motor's
torque curve, a gearbox's permissible torque and shaft loads, a bearing's
load ratings and life), and `(require NAME (within-ratings PATH))` asks
all of them of everything under PATH at once, against what the model
works out: the efforts, speeds and reactions of the motion written.

**The model's own numbers are facts.** A design no longer types the
torque its motor sees. Simulation works out each joint's effort, speed
and reaction from the bodies, gravity, loads and motion, and the ratings
read those. `(reaction PATH)`, `(effort-rms PATH)` and the rest are
values a requirement can use.

**Couplings.** `(relate A B KIND)` joins two parts' motion and passes
force between them: `gear-mesh`, `belt-drive`, `chain-drive`,
`lead-screw`, `cord`. A coupling that states no efficiency loses nothing,
and the answer says it assumed so.

**Simulation ships.** `sim/` is public for the first time: studies
(static, modal, buckling, steady and transient heat; a modal study of a
part held in part or not at all finds its elastic modes), topology
optimisation (`ink topology`) and dynamics (`ink effort`), as packages
built on the kernel the way a third party's would be. Flow and the least
time are held back (above).

**Units throughout,** as above: a quantity carries its unit through
arithmetic, and one read where another belongs is an error, not a number
a thousand times off.

A hoist, whole: a geared stepper, a 20 mm drum on its shaft, and 12 kg
hanging from the drum by a cord. `MOTOR.ink` stands for a rated motor
with a `drive`, from a part library (here a NEMA 17 stepper with a 5:1
planetary gearbox, 17HS19-1684S-PG5).

```lisp
(import "std")
(import "MOTOR.ink" as m)
(gravity 0 -9.81m/s2 0)
(t-span 0 1s)
(body drum (shape (cylinder-z 10 60)) (attr material "aluminum-6061"))
(body weight (shape (box 100 100 100)) (attr mass 12kg))
(assembly (hoist (angle 0deg))
  (instance motor (m.drive angle))
  (instance drum drum (frame-of motor shaft))
  (instance weight weight (at (move 10 -200 30)) (slide y (* 10mm (/ angle 1rad))))
  (relate drum weight (cord 10mm)))
(assembly root (instance h (hoist (* 10rad/s t))))
(require ok (within-ratings root))
```

`ink require hoist.ink --plain` answers (its notes and some of its lines left out here):

```
PLACEMENT: no two parts overlap at the start of the motion. Not checked here: the rest of the motion (`ink interfere FILE` proves it), and whether each part touches what it mounts on (`ink render FILE --out look.png` shows it)
ok: NO, estimated, at 18 values of t from 0 to 1: for the 4 ratings checked; 2 things not covered; 1 assumed; no clash at the start pose; not all yes  (hoist.ink:13)
  NOT COVERED: 2 of the parts the model loads under it carry no rating, so nothing here checks them: root/h/drum (a coupling passes force between it and root/h/weight); root/h/weight (a coupling passes force between it and root/h/drum)
  ASSUMING no loss in root/h/drum~root/h/weight (cord): it states no efficiency, so the efforts through it are the least the motion can need. To count its loss, write the share of power it passes last in its kind, after what is there: (cord ... 0.95), with the maker's number in place of 0.95
  the 4 ratings its parts carry, each checked:
  side-load on root/h/motor/rotor: NO, estimated, at 18 values of t  (17HS19-1684S-PG5 gearbox shaft-radial-load-N-max)
    the size of radial-load on root/h/motor/rotor shaft: 118.304 N, between 115.937 and 120.67 N, estimated  (ink effort on root)
    at most
    100 N
  thrust on root/h/motor/rotor: YES, estimated, at 18 values of t  (17HS19-1684S-PG5 gearbox shaft-axial-load-N-max)
  permissible-torque on root/h/motor/rotor: YES, estimated, at 18 values of t  (17HS19-1684S-PG5 gearing rated-torque-Nm: the maker's Max. Permissible Torque)
    the size of torque on root/h/motor/rotor shaft: 1177.2 N*mm, between 1177.186 and 1224.302 N*mm, estimated  (ink effort on root)
    at most
    2000 N*mm
  pull-out-torque on root/h/motor/rotor: YES, estimated, at 18 values of t  (17HS19-1684S-PG5 motor curves 0: pull-out, 24 V, 1.69 A, 2000 steps per revolution, DM542T; ...)
```

The torque is the hand sum, 12 kg x 9.81 m/s² x 10 mm = 1177.2 N mm, and
the motor has it. The no is the one thing nobody typed: the weight hangs
on the gearbox's shaft, 118 N across it where the maker allows 100. Two
bearings and a `(supports ...)` clause take that load off the shaft.

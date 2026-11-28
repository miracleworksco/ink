# Your first hour: a machine that moves

What to read first, in order, to build a machine with ink and trust what
it answers. Each step points into the long references for more:
[DYNAMICS.md](DYNAMICS.md) for motion, efforts and ratings,
[SIM.md](SIM.md) for a part's strength and stiffness, and
[CLAIMS.md](../../kernel/docs/CLAIMS.md) for how a requirement is answered.

## 1. ink, and a part to build with

Build `ink` as [../README.md](../README.md) says. `ink spec` is the language
and `ink std` the standard library; both are what an agent reads too.

A machine is made of parts somebody sells: a motor, a pulley, a rail.
Where a parts catalog is installed, `ink parts WORDS` finds them, each with
the import path a model reads it by:

<!-- not run: it needs a parts catalog, which this repository does not hold -->
```sh
ink parts 17HS15-1504S
```

The answer's `import` (and the frames each part has, where things attach)
is what goes in `(import "..." as m)`. Below, `MOTOR.ink` stands for that
path: replace it with the one the search gives.

## 2. The smallest machine

A belt axis: a stepper turns a 20-tooth GT2 pulley, the belt runs round an
idler of the same size and is clamped to a carriage, which slides.

```lisp
(import "std")
(import "MOTOR.ink" as m)            ; a rated stepper with a `drive`: a part library's NEMA 17, its shaft 23.5 mm
(gravity 0 0 -9.81m/s2)
(t-span 0 4s)
(def r 6.366)                        ; a 20-tooth GT2 pulley's pitch radius, mm
(def L 400)                          ; motor to idler, mm
(def bz 18.5)                        ; the belt's height: the pulley's middle, 5 mm under the shaft's end
(def t0 20N)                         ; the belt's pretension
(def travel (* 150 (- 1 (cos (* 90deg/s t)))))   ; the carriage: 300 mm out and back in 4 s, again and again
(def turned (* (/ travel r) 1rad))        ; the motor's angle that moves it there

(body pulley (shape (subtract (move (cylinder-z r 6) 0 0 -5) (cylinder-z 2.5 30))) (attr material "aluminum-6061"))
(body idler (shape (subtract (cylinder-z r 6) (cylinder-z 1.5 10))) (attr material "aluminum-6061"))
(body carriage (shape (box 60 40 30)) (attr mass 500g))
(assembly (axis (a 0deg))
  (instance motor (m.drive a))
  (instance p pulley (frame-of motor shaft))                        ; on the shaft: the motor turns it
  (instance i idler (at (move L 0 bz)) (turn z a))                  ; the same size: it turns as the motor does
  (instance c carriage (at (move 50 (- r) bz)) (slide x (* r (/ a 1rad))))   ; clamped to the run a +z turn moves +x
  ; the belt, a span at a time: pulley to carriage, carriage to idler, and back from the idler to the pulley
  (relate p c (belt-span (vec (vec 0 (- r) bz) (vec (+ 50 (* r (/ a 1rad))) (- r) bz)) t0 0.95))
  (relate c i (belt-span (vec (vec (+ 50 (* r (/ a 1rad))) (- r) bz) (vec L (- r) bz)) t0 0.95))
  (relate i p (belt-span (vec (vec L r bz) (vec 0 r bz)) t0 0.95)))
(assembly root (instance x (axis turned)))
(require axis-ok (within-ratings root/x))
```

Read it top down:

- `(gravity ...)` and `(t-span ...)`: the world's weight, and the motion's
  time. Every number with a unit says it (`9.81m/s2`, `2s`, `500g`); a
  length with none is millimetres.
- `travel` is the carriage's position over time, written as a formula of
  `t`; `turned` is the motor angle that puts it there. A motion is always
  written this way: what drives, as a function of time.
- A part's look and mass: `(body NAME (shape ...) (attr material ...))`, or
  `(attr mass ...)` where the shape is only its envelope.
- `(assembly (axis (a 0deg)) ...)` places the parts, the motor's angle its
  argument. What the motor turns rides its shaft: `(frame-of motor shaft)`.
  The idler turns on a joint of its own (`turn z a`), the carriage slides
  (`slide x ...`), each by a formula of the same angle.
- The belt is a `belt-span` a span at a time: the line it runs along, from
  where it leaves one part to where it meets the next, its pretension and
  its efficiency.
- `(require axis-ok (within-ratings root/x))` asks whether every rated part
  under the axis is within what its maker rates.

## 3. Ask, and read the answer

<!-- not run: it needs MOTOR.ink replaced by a catalog's path -->
```sh
ink require axis.ink --plain
```

```text
axis-ok: YES, estimated, at 52 values of t from 0 to 4: for the 1 rating checked; 4 things not covered; no clash at the start pose; not all yes
  NOT COVERED: 3 of the parts the model loads under it carry no rating, so nothing here checks them: root/x/p (a coupling passes force between it and root/x/c); root/x/c (a coupling passes force between it and root/x/p); root/x/i (a coupling passes force between it and root/x/c)
  NOT COVERED: a load on a rated part that none of its ratings reads, so nothing here checks it (its record states no limit for a load across or along its shaft or travel): root/x/motor/rotor (a coupling pulls on root/x/p, which it carries, across or along its joint)
  pull-out-torque on root/x/motor/rotor: YES, estimated, at 52 values of t
```

- The verdict is one of four: **YES** and **NO** (as stated, shown, or
  estimated: read at the moments named, each a value with its band), **not settled** (it
  says why, and what would decide it), and **not available** (nothing could
  answer it). A requirement over a motion is read at many moments; the
  answer says how many and where it has least room.
- **NOT COVERED** is what the yes is not about. Here the pulley, idler and
  carriage carry no rating, and the motor's ratings say nothing of the
  belt's pull across its shaft. A yes is only as wide as what it covers:
  read these lines before you trust it.
- Each row below says where its number came from: a maker's table, the
  motion's inverse dynamics, a formula. By hand the carriage needs 0.5 kg at
  0.37 m/s², 1.18 N mm on the pulley's 6.37 mm radius; ink says 1.58 with
  the pulleys' turning and the belt's loss. Check one number like this on
  every new machine.

## 4. What you still check by hand

ink answers what the model says. It does not know what you meant:

- **Units.** A missing unit is refused, saying which (`attr mass takes a
  mass, and 1200 is a bare number`). A wrong one that is still a unit is
  not: `1.2g` for a 1.2 kg camera, or gravity of `9.81mm/s2`, gives a
  confident yes. Read every mass and every unit once.
- **What holds what.** ink says when parts overlap where they start (the
  answer's first line, and `ink interfere` over the motion), when a part
  in a small assembly (up to 60 parts) is proven more than 1 mm from all
  others with no relation joining it (`LOOSE`), and when a part moves by
  a motor's angle with nothing between them (a belt left out): NOT COVERED
  names clashes and missing couplings, and the answer is not all yes. It
  still does not see a deleted bolt when the part rests against its mount:
  `ink render FILE --out look.png` shows it.
- **The motion's end.** A motion that stops by `clamp` or a `select` on
  `t` is read to its end. One that goes on past the span is said (`the
  written motion goes on after the span's end`): give a span that covers
  it, or a motion that stops or repeats within it.
- **What is not rated.** Everything NOT COVERED lists: a pulley's bore, a
  belt's tension, a bracket's strength. A bracket's strength is a study
  (SIM.md's quick start); a part with no rating needs its own number.

## 5. Where next

- A rail and carriage under the axis, its load and moments held against
  the carriage's ratings: DYNAMICS.md, Supports and Ratings.
- A geared motor, a lead screw, gears or a chain: DYNAMICS.md, Couplings.
- The force on each part, to size a bracket: DYNAMICS.md, Reactions; then
  SIM.md for its stress and deflection under them.
- How fast it can go, and what a motor's curve allows: DYNAMICS.md, How
  fast; `ink effort` for every actuator's torque and speed at every moment.
- What each verdict rests on, and how to make a not settled settle:
  CLAIMS.md.

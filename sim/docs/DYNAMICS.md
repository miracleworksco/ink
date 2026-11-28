## Dynamics

The effort a mechanism's actuators supply along its written motion, by
inverse dynamics: what each motor must give for the motion as written,
and whether it can; the load each joint and the model's support carry,
and whether it stays standing. Every answer here is an estimate unless it
says otherwise (Scope and limits). How each is checked against hand
calculations: docs/DYNAMICS-VALIDATION.md.

### Quick start

A motor swinging an arm, with its hand calculation in its comments:

<!-- example examples/arm.ink -->
```lisp
; A motor swings an aluminium arm 45 degrees each way about a level shaft.
; What torque must it give? By hand: the arm (200 x 20 x 10 mm, 108 g) held
; level takes m g L / 2 = 0.108 x 9.81 x 100 = 105.9 N mm.
; ink effort examples/arm.ink: 105.9 as it passes level, and its peak,
; 106.9, just below level, where it also slows the arm's swing. ink require
; checks the motor's limits (+-150 N mm) along the whole motion: they hold.
(import "std")
(gravity 0 0 -9.81m/s2)                      ; m/s^2, z up
(t-span 0 2s)                             ; s

(body arm (shape (move (box 200 20 10) 100 0 0)) (attr material "aluminum-6061"))

(assembly (joint (angle 0))              ; angle: deg, an argument the motion passes
  (instance arm arm (at) (turn y angle) (effort-limits -150N*mm 150N*mm)))   ; the motor: N mm
(assembly root
  (instance j (joint (* 45deg (sin (* 180deg/s t))))))
```

Run from sim/ (well under a second):

<!-- run -->
```sh
ink effort examples/arm.ink
```

```json
{"actuators": [{"path": "root/j/arm",
                "peak": {"effort": -106.9, "t": 0.055},
                "rms": 91.2,
                "margin": {"margin": 43.1, "lo": -150.0, "hi": 150.0}}]}
```

- `actuators`: each joint with `(effort-limits ...)`, here the motor at
  `root/j/arm`. Its `peak` effort over the motion is -106.9 N mm (by hand,
  105.9 to hold the arm level, and 1 more to slow its swing), at
  0.055 s (as much at 0.945 s, its mirror: the earlier is shown);
  negative, as it turns the arm about -y, holding it up. `rms`
  is what heats a motor; `margin` is how near the effort comes to its
  limits, at its nearest.
- `--t 0` gives one moment instead: -105.9 N mm as the arm passes level.
- An effort is an estimate (`accuracy: null`): measured at moments
  through the motion, closely enough to find its peak (Steps).

Each actuator's limits are a claim, which `ink require` checks:

<!-- run -->
```sh
ink require examples/arm.ink
```

```json
{"counts": {"yes": 1},
 "requirements": [{"name": "effort-limits:root/j/arm", "verdict": "yes", "value": 43.1}]}
```

Its value is the margin at its nearest, 43.1 N mm to spare, less the
0.015 the arm's measured mass, centre and inertia can move the effort
within their proven bounds: 43.08. A motor too weak,
`(effort-limits -100N*mm 100N*mm)`, fails, with a moment it falls short (`t`).

### Forms

```lisp
(gravity 0 0 -9.81m/s2)                      ; m/s2: the world's, z up
(instance NAME CHILD WHERE [CFRAME] (turn z a) (effort-limits LO HI [(accuracy x)]))
(instance NAME CHILD WHERE [CFRAME] (turn z a) (actuator [LO HI]))   ; driven, and what it can give, claimed by nothing here: its part's rating is what is asked
(instance NAME CHILD WHERE (slide x s) (stiffness K))   ; a joint held by a spring (Modes, Simulation)
(attr spin-inertia 5400g*mm2)                 ; g mm2: a body's published moment about its own z
(load NAME PATH (force fx fy fz) [(torque tx ty tz)] [(at FRAME)])   ; N and N mm, the world's axes; in an assembly, its own
(reaction PATH [FRAME])                  ; a value: (vec fx fy fz mx my mz), forces and moments
(reaction PATH (by S ...))               ; a value: each part S's share of it, (vec W ...)
(supports PATH S ...)                    ; the parts S hold what joint PATH carries: each is given its share (Ratings); in an assembly, its instances' paths
(zmp R)                                  ; std: where a support's reaction R acts on its plane
(slip R)                                 ; std: its force across the plane over its push
(tension A B)                            ; a value: each line of A's and B's coupling's tension, a force
(effort-rms PATH)                        ; a value: the RMS over the span of actuator PATH's effort, a torque or a force
(simulate X0 [V0])                       ; a coordinate's motion, simulated forward from its forces
(command PATH EFFORT)                    ; actuator PATH's effort in a simulation, clamped to its limits
(springs)                                ; the file's stated stiffnesses simulated as springs
(following-error PATH)                   ; a value: a sprung joint's deviation in the simulation
```

Units: lengths mm, angles deg, masses g, forces N, torques N mm, time s,
gravity m/s2; an effort is N mm on a turn and N on a slide, its rate
deg/s or mm/s.

- Its clauses are `(effort-limits ...)`, `(actuator)` and `(stiffness K)` on an
  instance with a motion, after its limits, and `(couple ...)` in a
  relation (Couplings); an ink without this package refuses them, naming
  them (IR.md, Data for other packages).
- `(gravity gx gy gz)`, in m/s², is the world's (`(gravity 0 0 0)` too).
  Only the file run's counts; twice is an error. There is no default:
  `ink effort` without one is an error. Studies keep their own.
- `(effort-limits LO HI)` on a joint makes it an actuator, states what
  it can give and claims it. `(actuator LO HI)` states the same and
  claims nothing: the solver knows what the actuator can give (the
  effort several actuators on one coordinate share, a command's clamp,
  `ink effort`'s margin), and what is asked of it is its
  part's rating. `(actuator)` states no limit: its effort is worked out,
  with no margin. A catalog motor whose record carries its torque as a
  rating (a stepper's pull-out curve, the Kraken's torque line) writes
  `(actuator LO HI)` from the same record, and the rating is asked with
  `(require NAME (within-ratings PATH))`: one check of the curve, and a
  design placing the motor and asking nothing is asked nothing. Every
  other joint supplies no effort. LO and HI are fields of params and t, so they
  may read the drive's rate: a torque-speed line is an expression (std's
  `dc-torque-max` and `dc-torque-min`, a DC motor's at full voltage, capped
  at a current; `curve-at`, a published curve through its points). Limits
  that cross (LO > HI) hold no effort: where an actuator's data ends. A
  file with gravity claims each actuator under its root,
  `effort-limits:PATH` (its instance path, `root/arm/s/output`):
  `min(effort - LO, HI - effort) >= 0`, with `terms` effort, rate, lo and
  hi. The run file owns the claim and supplies its motion, even when the
  joint or root is imported. Without gravity
  `ink check` names each unchecked limit under the root in `not_checked`.
  `ink require` includes these placements, though their parts are imported. An unchecked
  limit in the requested scope clears both hold summaries and gives exit 3. See Steps for how a
  step decides a limit when gravity is given.
- `(attr spin-inertia J)` is a body's published moment about its own z
  axis, and states it a rotor balanced about that axis (a datasheet's rotor
  inertia, where the catalog models the rotor as its keyed shaft). It
  replaces the measured moment about that axis; the mass and other moments
  stay measured. Spinning about that axis, the body is read as round
  (Steps).
- A catalog drive is asked what it can give once. A motor whose record
  carries its torque as a rating is asked by that rating, `(require NAME
  (within-ratings PATH))`: the catalog's generated steppers (their
  published pull-out curve, read at the shaft's speed at each moment;
  below the chart's first point as its torque there, where the record's
  holding torque is no less, and otherwise not available there, saying
  so; not available past its last speed), the FRC motors and the servos
  (each one's line from stall torque to free speed, at its record's
  voltage) and the Pololu gearmotors and DC motors (that line,
  and what their gearbox takes: the root mean square of the torque over
  the duty at most the torque Pololu recommends for a continuously
  applied load, `continuous-torque`, and the torque at every moment at
  most its momentary one, where its record gives it). Its joint is `(actuator LO HI)`, the
  same curve from the same record for the solver, claimed by nothing. A
  design that sets the drive's `volts` (a stepper's supply curve) or
  `amps` (a DC motor's or servo's current cap) passes that argument to
  the rotor's ratings (`pull-out-torque` reads the matching supply curve,
  and `current-limit` caps torque by `amps`), so the rating itself holds
  the drive's setting. The others state
  `effort-limits` always: a rod actuator's line through its
  no-load and full-load points, capped at its rated force; and the five
  hand-written NEMA steppers' pull-out curve at their supply `volts`. Its rotor's inertia is its
  `spin-inertia`. Where the record gives no mass, the moving body is
  steel, an estimate its record lists.

**Coordinates.** A mechanism moves by the numbers the target's instances
pass as arguments (root's by default) that some joint's drive reads, at
any depth: its coordinates, named `INSTANCE.ARG` (`arm.shoulder`). An
expression of t is a motion, its speed and acceleration its exact
derivatives; a constant is a pose held still. An argument no joint reads
(a length), or a constant a body's shape reads too (a lead read by a
thread and a drive), is a size: it does not move. Write the motion as
arguments, `(assembly (arm (shoulder 0) (elbow 0)) ...)` placed by
`(instance arm (arm (swing 50 3) -30))`; every body whose placement reads
a coordinate moves with it, however the expressions couple them (a rotor
turning 48.8 times an arm, a CoreXY's carriage and two motors on X and
Y).

**The effort.** Each moving body's weight and inertia, from its mass
properties as `ink mass` measures them, are projected onto the
coordinates through exact derivatives of its placement (Kane's method):
the effort along q is the sum over bodies of ∂x/∂q · m (a − g) + ∂θ/∂q ·
(I α + ω × I ω). The actuators' efforts τ solve Bᵀ τ = f, B the rate of
each actuator's drive along each coordinate (radians for a turn), by
least squares of least sum of (τ / span)², span an actuator's limits'
width there: actuators alike share alike, a larger one takes more (as its
span squared). This split is a fixed policy, not what the drives' own
stiffness would give. Drives that share no coordinate are solved
independently. The solve scales each coordinate to its own drive
coefficients and decides independence on the unweighted directions.
A direction too close to separate is shared only if its dependence can
be established. Otherwise the group is undecided. A fixed-limit drive
takes its stated effort when the other drives can supply the balance.
Where that is impossible, the required departure from its fixed effort
is reported, shared by least squared departure among fixed drives. A
nonzero difference within the rounding of that drive’s load leaves its
equality claim undecided. A tension or reaction that uses this balance is
refused for a fixed-effort discrepancy only beyond the same rounding.
A moving coordinate or one held by a declared drive reports any effort
left unsupplied as `unpowered`, and every actuator depending on that
group of coordinates is undecided. A passive
joint fixed in place carries its load.

**Carried bodies.** A body whose `carried` says its material moves (its
path or shift reads t or a moving coordinate, its profile neither: Paths)
is counted as that material: station m at the path's station m + SHIFT,
its mass per mm the profile's (measured once, spread evenly along the
pitch line), its effort integrated along the material (Gauss–Legendre on
each piece). So a belt counts where its material goes: a CoreXY belt adds
its mass to each axis and couples x and y. Not counted: the section's
offset from the path and its turning round an arc; the change of speed
where an arc advancing other than its radius per radian (a chain on a
sprocket) meets the next piece; stretch, sag and tension; material
entering or leaving (a length that changes).

**Loads.** `(load NAME PATH (force fx fy fz) (torque tx ty tz) (at
FRAME))`, the file run's (an imported file's top-level load is on that
file's root, which an importer does not place: it is not applied, and
where the importer places what it is on, `ink check` lists it apart,
`loads_not_applied`, and `ink require` says it is not covered and is
not all yes; a load that goes wherever its assembly is placed is
written inside the assembly, below), is a force in N and a torque in N mm, components in the world's
axes, on what PATH places (a body, or an assembly nothing moves inside,
as one piece; one whose parts move inside is refused, naming them),
acting at FRAME's origin, read at each moment, by default PATH's own.
`(at (frame-of root/arm/fore tip))` is a point moving with the part;
`(at (move 0 0 150))`, written as an instance's place is, a point the
force stays at while the part moves under it. A load's name is its own:
the same name stated twice is refused (it would be summed), and `ink
check` lists every load the model applies with its force at the moment asked. Components are fields of params and t and may read
frames: a spring between two parts is a load on each along the line
between their frames. Loads, placements, motions, limits and attributes
may not read a reaction (a friction load as μ times a normal force is
refused at load; whether a support slides is std's `(slip R)`: Tipping).

In an assembly, `(load NAME PATH ...)` beside its instances is a load its
parts exert, PATH an instance it places, components in its axes, FRAME
in its coordinates, all fields of its arguments, params and t. Each
instance of the assembly carries it, turned and moved as placed: a catalog
strut or a sprung hinge states its force once (`ink effort`, claims,
reactions and simulations alike). A load supplies effort along each
coordinate through its point's velocity (Kane's generalized force), so
the actuators need that much less, and it is in every reaction carrying
its path.

### Reactions

What the model works out carries its unit: a reaction's forces are forces
and its moments torques, a tension a force, an RMS
effort a torque (a turn's) or a force (a slide's). A design never types
the unit of a number the model produced; what it compares one with says
its own: `(value (nth (reaction root/x) 2) < 200N)`, `(value (effort-rms
root/m/rotor) < (* 1N*m rated))`.

`(reaction PATH [FRAME])` is a value, `(vec fx fy fz mx my mz)`: what the
rest of the model exerts on the bodies PATH carries, N and N mm about
FRAME's origin in its axes, at each moment. PATH (root, or an instance
path under it) carries the bodies under it and those placed on what it
carries (their WHERE reading its frames: in `(instance upper upper
(frame-of s output) hub)` the upper arm rides the servo's output). Every
body it carries counts, moving or still: the sum of each one's m (a − g)
and I α + ω × I ω, less the loads on them. FRAME is any frame, read at each
moment (`(at)` the world's); without one, PATH's joint's frame (its WHERE
moved by its motion). The joint keeps its declared local axis. A slide on x
uses force component 0 for its axial load, a slide on y uses 1, and a
slide on z uses 2. The corresponding moments are components 3, 4 and 5. A
jointless instance uses its placement. Root uses the world's axes. In a
joint's frame the moment about its axis is the torque the joint passes,
the force along it its axial load, and the force and moment across it the radial load and
tilting moment its bearings hold.

That is the joint's load where PATH meets the rest through that joint
alone. A coupling's lines and frictions are not in it (Couplings); a belt
not declared one, a second support or a contact reaching what it carries
is in the sum, unseparated (Supports divides it). A belt carried along its
path (std's `carried`) whose every piece lies on a coupling's span or
wraps a part a span meets is separated: its material's inertia is on the
parts it lies on (a wrap's on its wheel; a span's on its two ends, at
them, by where it is along it, as a string's load is; an open belt's tail
past a clamp on the clamp), so on what carries them. `ink effort` names the bodies
moving with its coordinates that it does not carry (place a rider by
`(frame-of ...)` to have it carried) and those placed partly on what it
carries (`moves_with`): the coordinates its joint's drive reads, or, on
no joint, those its own placement reads (a linkage's rocker placed by a
formula of the crank's angle). Where any is listed, the reaction sums
their push on what it carries where no joint is declared (a linkage's
pins), which the model does not determine: it is refused there, and a
claim reading it is `unknown`, saying so: a linkage's pin loads are not
found (Scope and limits). One carrying a body placed by a formula of t
through no coordinate or joint (in `driven_by_nothing`: a linkage's rocker
placed by a formula of t; a body riding a joint, a fan spun on its own
joint by t or a disc on a joint whose frame t moves, is held by that
joint) is refused too: nothing declared moves or holds that body.
Where any such body is in the model, every claim on an effort or a
reaction names it in its `notes`: no actuator is charged for its motion,
and where it touches what the claim reads its push is not in it.

A claim reading a reaction is measured at steps of its own: the efforts'
rules (Steps), what only the reaction carries, and each component in the
claim's frame within 4% of its peak of the line through its neighbours.
Bracketed minima and maxima near each component's extremes are located in
that same frame and measured in addition to the adaptive grid.
A group spinning with one link that it carries is read with its spin at
rest, round, where its bound on the reaction is within 1% of the
reaction's peaks, or where its spin is past the step budget; the claim
then holds where it holds for every value within that bound, fails where
it fails at the moment's real value, and is undecided between, saying so.
`--t` reads the real value. A claim reading several reactions, or a
reaction and an effort, is measured at every step any of them chooses.

It is `unknown` where a body it carries has no mass (named: no density is
guessed) or is not proven bounded (a ground half-space with a density),
deforms, or moves along a coordinate that is unpowered there. It needs
`(gravity ...)`. One an imported file asks for is answered in that file,
with its gravity and loads: the importer's check is `unknown`, saying so.

### Supports

A part on a joint is held by that joint: its load is the joint's
reaction, and nothing more is written for it (a carriage on its slide, a
motor's rotor). `supports` is for a body held by parts that are not its
joint (a shaft in two bearings, a table on four carriages). A part is
never its own support (an error), and a support is at what it holds: one
off a turning joint's axis by more than it is across, or clear of what a
slide carries by more than they are across together, is refused.

`(supports PATH S ...)` says what holds a joint's load, where the design
is written: at the top level with paths from root, or in an assembly
with its own instances' (`(supports motor/rotor ba bb)`). Each S is then
given its share of PATH's reaction as a fact its ratings read (Ratings).
The same division is a value where a claim reads it:

`(reaction PATH (by S ...))` divides PATH's reaction among the parts that
carry it, each S an instance path (a bearing, a carriage): `(vec W ...)`,
one per S, W `(vec fx fy fz mx my mz)` what S exerts on what PATH
carries, N and N mm about the origin of S's `load` frame in its axes, or
of FRAME where S is written `(S FRAME)`. Each rated catalog part declares
its `load` frame, `(frame-of S load)`, where its rating acts: a bearing's
centre with z along its bore (an angular-contact or tapered bearing's
effective load centre), a bushing's centre, a flange or pillow block's
insert centre, a carriage on its rail under its centre with x along it.
A part declaring none is an error naming it (an assembly placed around
one such part stands for it).

What it carries is its record's `load-ratings`: its body's `(attr spec
RECORD)` (examples/slider-loads.ink), else the `spec` its part's file
reads (every rated catalog part's): a guide carriage's force across its x
and moment about each axis; a thrust bearing's or screw nut's (`direction`
axial) force along its z, a push along +z alone unless its `thrust` is
double-direction (a pull is refused, never checked as a push); a radial
bearing's or bushing's (`direction` radial) force across its z.

PATH's joint is frictionless: no part carries force along a slide or
moment about a turn's axis, and its drive passes that; where an axial part
(a screw's nut) lies along the slide, it passes the drive along its own
line, and the carriages take that line's moment. On a slide each part is
set to its motion; one at a slant is refused. Where no nut is named, a
claim on a carriage's share says in its `notes` where the drive is taken
to pull (a belt clamped h above moves h times its pull of the moments; a
belt related as a coupling pulls along its own line).

The shares are the reaction's six equations solved: two bearings share a
shaft's radial load by statics over their spacing, one carriage takes the
force and three moments. Guide carriages of one record, two or more on
parallel rails, share as HIWIN's and THK's catalogs size them (HIWIN
G99TE24-2410, Table 1-3): a rigid table on identical blocks of equal
stiffness, their moments carrying only what their forces cannot, their
forces of least sum of squares (P = W/4 ± W a / 2c ± W b / 2d down each of
four). That is the catalogs' assumption, which a claim's `notes` state,
not the model's stiffness. A claim reading a share is `unknown`, saying
which, where statics does not determine them (three bearings on a shaft,
carriages of different records, two nuts on one slide), and where the
parts named cannot carry all of the reaction (an axial force on two radial
bearings: name the thrust bearing), however large what they do carry
beside it: what none carries is weighed against the arithmetic's own
rounding of what is summed (4096 times a double's 2.2e-16), not against
the reaction's size, so 1e12 N across two bearings hides no 20 N along
them. A share is never guessed. A support
whose record has no load-ratings, or neither a direction nor a guide's
moments, is an error naming it. A share is an estimate as the reaction
is, and carries its bound where the steps read spinning groups round. A
coupling's lines are out of the reaction, so the shares carry them (a
belt-driven shaft's bearings its belt's pull); an axial support that is
one of a coupling's parts (a screw's nut) carries the coupling's force
along its axis besides its share.

### Couplings

```lisp
(relate A B (belt-drive ra rb t0 [e [eb]] [(tensile-stiffness ea)]))   ; std: an open belt, pitch radii ra, rb (mm), t0 (N) a span
(relate A B (chain-drive ra rb [e]))     ; a chain: its slack span pulls with none
(relate A B (belt-span span t0 [e [eb]] [(tensile-stiffness ea)]))  ; one span of a belt, (vec p q): p leaving A, q meeting B
(relate A B (gear-mesh za zb mo [pa [backlash [e]]]))   ; spur gears; bevel-mesh, rack-mesh alike
(relate A B (lead-screw lead [e [eb]]))  ; A the screw, B the nut
(held t0)                                ; std: a pretension a constant-force tensioner holds
```

A coupling's two parts move by one coordinate: the driven part's motion
is written as a function of the driver's, never as a coordinate of its
own passed in step (that is refused: "only one of its two parts moves
with that coordinate"). A motor, a lead screw and its nut, the shortest
whole machine with one (MOTOR.ink stands for a rated gearmotor's import
path, the Pololu 4751's: `ink parts 4751` gives it where a parts catalog
is installed):
the screw sits on the motor's shaft, the nut's slide reads the motor's
angle, and the coupling says what passes between them.

```lisp
(import "std")
(import "MOTOR.ink" as m)
(gravity 0 0 -9.81m/s2)
(t-span 0 2.7s)
(param lead 8mm)
(def turned (trapezoid-position t 9000deg/s2 1800deg/s 4500deg))    ; 100 mm up at 40 mm/s, 0.2 s ramps
(body screw (shape (move (cylinder-z 4 200) 0 0 100)) (attr material "steel"))
(body nut (shape (subtract (cylinder-z 11 15) (cylinder-z 4.2 20))) (attr material "steel"))
(body load (shape (move (box 100 60 50) 0 45 0)) (attr mass 3kg))
(assembly (lift (angle 0))
  (instance motor (m.drive angle))
  (instance screw screw (frame-of motor shaft))                      ; on the motor's shaft: it turns with it
  (instance nut nut (at (move 0 0 60)) (slide z (* (/ angle 360deg) lead)))   ; the nut's slide reads the motor's angle
  (instance load load (frame-of nut))
  (relate screw nut (lead-screw (/ lead 1mm) 0.4)))                 ; 8 mm a turn, 40% efficient
(assembly root (instance z (lift turned)))
(require lift-ok (within-ratings root/z (excluding root/z/screw root/z/nut)))
```

The motor's torque is the load's weight and what speeding it up takes,
through the lead and the coupling's efficiency (3 kg on 8 mm a turn at
40%: 96.9 N mm). The screw and the nut here are plain bodies that carry
no rating, so the requirement leaves them out on purpose and says so.
It is yes for the motor and not all yes: nothing holds the screw but the
gearmotor's own shaft, whose record states no thrust limit, and the
answer says that load is checked by nothing. Bearings named as holding
the rotor, `(supports motor/rotor ba bb)`, take it.

What hangs from a drum is joined to it by std's `(cord R)`: one line from
where it is tied to where it leaves the drum, which pulls and never
pushes. A hoist is then a motor, a drum on its shaft, a weight whose
slide is the drum's turning, and `(relate drum weight (cord 10mm))`: the
motor's torque, its shaft's side load and the cord's pull all answer.

A relation's `(couple LINES T0 [E [EB]])` (IR.md: Assemblies, Relations)
is what A and B exert on each other: each line `(vec P U [H [L [K]]])`, P
a point on it and U its direction, pulls B along U and A the other way
with a tension T >= 0 and a moment H T about U (H, mm, its pitch: a
screw's lead / 2 pi), L its length, K its stiffness (N per mm its ends
draw apart; not given, rigid), in the relating assembly's coordinates.
std's `(tensile-stiffness EA)` (EA a force, its unit said, 16903N: a belt
maker's tension for 0.1% elongation, times 1000) makes each line a spring of K = EA / L (a line
with no length, or one stating its K, is then an error). A is the input
(the motor's side). std's kinds give a belt's or chain's two spans (outer
tangents of the pitch circles about each frame's z), one span of a belt
routed round several parts (std's `belt-route-span`), a gear mesh's two
lines of action through the pitch point at pa, and a screw's axis either
way at pitch lead / 2 pi. Each line must do no work along any coordinate
(its ends move alike along it, to 1e-3 of their motion's terms); else the
coupling does not match the joints' motion (a crossed belt written open, a
ratio's or lead's sign), and its tensions are refused, naming the line
and the coordinate.

**Tensions.** A passive joint carrying one end of a coupling passes
nothing along its motion: what its bodies need along it is the lines'
pull. Those balances, solved at each moment, give the tensions; a tension
they leave free (a belt's spans, a gear's two flanks) is set as installed.
A plain T0 is installed and locked (Gates' Light Power & Precision
manual, V.D): the spans' tensions less T0, weighted by length, sum to 0,
so two equal spans carry T0 ± half the torque over the pitch radius, and
a span that would push goes slack at 0. `(held T0)` is a constant-force
tensioner: the slack span keeps T0, the tight one T0 plus the torque over
the radius. A chain's slack span pulls with none. A gear's tooth force is
the torque over the pitch radius along the tangent, with Ft tan pa
separating them (a bevel's at its mean radius, AGMA 2003); a screw's
thrust is the nut's force. `(tension A B)` is `(vec T ...)` in its lines'
order; `ink effort` gives each coupling's `tensions`, at a moment or where
its tightest line peaks.

A claim reading them is `unknown`, saying which, where a line no balance
reads (two actuators driving one train), where the balances cannot all
hold (a part rigidly driving a joint through no declared coupling: place
it on its driver's frame), where a tension is left free beyond the
pretension (two belts in parallel), where a line would push, where a
balance's coordinate moves with no actuator supplying it, or where a
coupling of its train passes power and gives no efficiency (its loss is
in the tensions: Efficiency). A train is the couplings whose tensions rest
on one another's (a balance reads both, or they share a part): a balance
of one train refuses no other's, and machines sharing no coupling stand
apart. A fixed-effort drive is checked against the effort the prescribed
motion requires. If they differ, its effort claim reports the mismatch
and a tension or reaction depending on that balance is undecided; where
they differ within the rounding of the drive's load, its effort claim is
unknown, and a tension or reaction resting on that motion is read to the
shortfall's share of its own rounding either way: a value claim holds
where it holds over that band, and is undecided where its margin is
within it. A belt's own material lying on its spans and wraps is in the
balances (Reactions), a span's tension its mean along it. What else moves
with a balancing joint's coordinate that no actuated or balancing joint
carries (a coupler on a belt-driven crank, a deck sliding beside a
screw's nut) is left out of the tensions, and `unseparated` names it. A
loss resting on them can be off by at most its worst efficiency's share
of what those bodies need (each stage's worst, multiplied along the stages
moving it), and of the friction charged where several stages meet: an
effort claim is decided only where its margin passes that
(`loss_bound`). A group spinning under a balancing joint with
a departure from round is resolved, never read at rest: its bound on a
tension is not found.

**Reactions.** A reaction leaves out the lines and frictions of a
coupling on what its path carries: a belt-driven shaft's joint passes no
torque, its bearings carry the belt's pull; a carriage a belt is clamped
to carries the belt's line; a motor's shaft passes its effort, the loss
among it. What moves with a coupling's part through accounted joints is
separated, and so is it from a part on what such a joint carries that no
coupling reaches (a servo bolted to a belted carriage: the belt pulls the
carriage, so the servo's reaction is its own bodies'); the rest stays in
`moves_with`.

A line with a length (a belt's span) meets an end that turns within
that wheel's width along its axis: one drawn more than a millimetre past
its faces is refused, saying where (a belt beside its pulley).

**A study's load from the motion.** A study's `(force ...)` or
`(torque ...)` may read what the model works out, a reaction or a
tension, with params; nothing is typed, and the study follows the
design:

```lisp fragment
(def shaft (reaction root/x/m/rotor))             ; the belt's pull on the motor's shaft
(study hold ... (load belt (force (nth shaft 0) (nth shaft 1) 0N) (region screw-rings)))
```

It is a number when the study is asked, not when the file loads. A
claim on the study answers for the whole motion, and no better than it
would for any part of it. `ink require` reads the model's loads at the
span's ends and the middle of each step the model's own answers are read
at, which follow how fast its motion and loads change (above; 64 even
steps where nothing chooses any). A step the model does not follow (a
load that swings faster than its finest step, a jump) leaves the claim
not settled, saying where: a yes at the moments read says nothing of it.
A static elastic study's
deflection and stress fields are linear in its loads, and every moment's
loads are a weighted sum of a few moments': as many as the loads have
independent directions over the motion (one for a pull that only grows
and shrinks, six at most for one force and torque; one more where a
typed load or gravity acts too, which every moment holds once). The
study is run at those few moments, at once where memory leaves room,
and each moment's fields are that sum of those runs'. The claim is then
decided at every moment whose loads differ, each as a run of its own
would be: its value, the stress by the supports that it leaves out, its
read's error, a sharp corner, and its own comparison of three grids.
Where the model's loads are all that drive the study and every moment's
are one moment's times a number (a pull that only grows and shrinks), the
response is that moment's times the number and every size a claim reads
is largest where the number is: the claim is decided there, and where it
is most negative, not at each moment.
What the model gives only within a band (a reaction of a spinning body
read round within what its measure is known to; a value whose own
measure leaves an error) is carried into the answer: each load is that
value times a number, so the answer moves by the band's share of the
study's response to that value alone, itself a sum of the runs (one more
run where no moment's loads lie along it). The claim holds where it holds
over that band and is not settled where the band crosses its bound,
saying the band comes from what the model gives. No band is dropped for
its size.
The answer is the least favourable moment's: a no before a not settled
before a yes. A smaller load in a direction the part bears worse is so
counted exactly.

A buckling factor is not a sum, and is bounded instead: one over it
grows no faster than the loads add, so a moment's is at least one over
the sum of each run's share, the runs' loads reversed where a weight is
negative (those runs solve their buckling both ways). A moment whose
bound reaches 1 is not settled, saying the bound; the runs' own moments
keep their own factors. A part cut at a mirror has no bound, and is not
settled over a motion.

Only the model's loads may change over the motion. Where anything else
in the study follows `t` or the motion (a pressure, a material, a
support's move, a region, the shape), a moment is no sum of others, and
the claim is not settled, saying so. A load is taken as a sum of the
picked moments' where it differs from one by under a billionth of the
largest (the rounding of the model's numbers). Between the moments read
nothing is bounded, and the answer says so. A metric that is not linear
in the loads (a buckling factor, a frequency, a compliance) is given at
the one moment the loads are largest together and is not settled over a
motion, saying why. `ink solve` shows one run, at that moment
(`run_with`). What the model cannot give (no gravity, a body with no
mass) refuses the study with what is missing.

**The books.** `ink require` adds up its own accounting from the numbers
it has worked out, at every moment read, and an answer whose books do
not balance is not a yes (`books_do_not_balance`, saying which part, how
much and which way): the shares handed to a joint's supports, summed in
the world, against the joint's reaction, force and moment; each
support's offered load against every share any joint gave it (a bearing
named as holding two shafts carries both); and the power the actuators
supply against the rate of work the motion takes along its coordinates
(equal with no efficiency stated, no less with one). And each joint's
own balance, added up apart from how its reaction was found: the weight
of every body it carries, the loads written on them, what couplings pull
on them with from outside, and their masses times the accelerations of
their centres (second differences of where the centres are, a little
before and after), against the reaction's force, at up to 33 of the
moments; and the same of moments about the world's origin: the
reaction's moment, the weights at their centres, the loads at their
points and the couplings' lines against the rate the carried bodies'
angular momentum changes at, their inertias turning with them
(`books_balanced` counts the joint-moments so balanced, `forces` and
`moments`: where a rotor spins too fast for the differences, its
moments are not counted). Not added up: a joint that carries a belt
counted along its path, or a part a carried belt lies on, or a body whose
shape moves; a moment where the motion's acceleration jumps (a ramp's
start).

**Efficiency.** E is the share of the power A gives that B receives, EB
of what B gives that A receives (default E). A screw's EB is by default
2 − 1/E (its thread's friction to first order); at or below 0 it
self-locks, and lowering a load the motor drives it down. The loss is
friction the coupling exerts on the part power enters it by, whichever
end is written first (E where it enters by A, EB by B), along the joint
that part rides, so what leaves is that share of what enters: a screw's
drive torque is F lead / (2 pi E) and back-driven F lead EB / (2 pi)
(Rexroth, screw assemblies, eq. 9 and 10). A single coupling's loss is on
the power through each end's own joint. Where an end's axle rides a
joint that moves, that power is the parts' motion in the carrier's frame
(Willis's). A turning carrier (a planet's) is so read: an end it carries
turns there by its own joint, and another coaxial with it, or held, by its
own turn less the carrier's (a sun held still turns back at the carrier's
rate), and the loss is charged there: a 0.98 mesh on a carrier driven
round a held sun needs 1010.2 N mm where the ideal is 1000. Where the
carrier slides (a pulley on a carriage rolling on a fixed span), or the
other end turns about another axis, the loss is not found: such a
coupling with an efficiency under 1 is taken as one that gives none, its
efforts the least the motion needs and said as assumed, the reactions and
tensions resting on it not settled. A routed belt's (one
mode of `belt-span`s) on the power each part gives it along each joint
moving it, its own and those carrying that (a gantry's slide pulling its
idlers' axles), but one moving every part of it. The friction is the
actuators' to supply, and the parts': it is in the reactions carrying
the part it acts on (a motor's shaft passes its effort), and where that
part's joint balances couplings' tensions (an intermediate shaft, a
passive output driving back), in its balance, so the stage driving it
carries it, in its efforts and its tensions: each stage's loss is on the
power entering it, its tensions what passes through it, and stages in
series multiply, either way power flows (two 0.5 meshes lifting need the
ideal torque over 0.25, the first's tooth force twice its ideal one;
lowering, the motor brakes 0.25 of it). Where several stages meet in the
joint, or a routed belt's part meets another stage there (a head two
belts drive: they work beside each other, not in series), its friction is
in the joint's balance, but the stages driving it do not pass it on in
the efforts: what they would add is in `loss_bound` (above).

Newton's third law pairs the friction with the coupling's other part. A
screw's is its thread's, a moment about its axis on the nut as on the
screw, whichever end power enters by: raising, the nut's guide holds the
motor's F lead / (2 pi E), and lowering F lead EB / (2 pi), the thrust the
load either way. Otherwise it is the least force that leaves the other
part nothing along its own joint: a moment alone where one does (a rack's
pinion's on the rack; across a bevel pair, about the driver's axis on the
wheel), else a force through the other's axis level with the lines,
across the line between the axes: two wheels on parallel shafts C apart
hold the friction's moment as a couple of that over C on their bearings,
as the frame must (a belt drive's 1000 N mm at 0.9 puts 1.11 N across
100 mm). A mesh's sliding friction along its line of centres, which E does
not give, is not in them. A routed belt's friction has no one other part:
what reacts it through its other parts' supports is not found and not in
their reactions, which say so. Through a joint carrying several of its
parts (a gantry's idlers) it is shared among them by the size of each
one's pull along the joint, so none takes more than the joint's loss.
A coupling whose two ends lie on one line (two slides along it, two
wheels on one axle) moves its ends as one along it, with nothing turning
or bending in it, and loses nothing, whatever efficiency it states.

A coupling giving no efficiency that passes power leaves the old effort
claim (`effort-limits`) `unknown`, and the reactions and tensions of its
train, which carry its friction; a requirement of ratings takes it as
lossless and says so once ("ASSUMING no loss in ..."), its efforts then
the least the motion can need. A synchronous belt drive's published bound
is "greater than 95%" (SDP/SI Technical Section, T-22).

### Ratings

A rated part's limits are its own (its `ratings`: CLAIMS.md), asked with
`(require NAME (within-ratings PATH))`. A design types no load: it says
how the machine is put together, and each rating is answered from the
model. A hoist: a geared stepper turns a drum that winds a cord, and two
bearings hold the drum. (The two imports stand for rated parts' import
paths, the 17HS19-1684S-PG5 geared stepper's and the 608ZZ bearing's:
`ink parts 17HS19-1684S-PG5` and `ink parts 608ZZ` give them where a
parts catalog is installed. With them every rating below answers yes.)

```lisp
(import "std")
(import "MOTOR.ink" as m)            ; a rated motor with a `drive`: a part library's geared stepper
(import "BEARING.ink" as b)          ; a rated bearing: a part library's deep-groove ball bearing
(gravity 0 -9.81m/s2 0)
(t-span 0 1.3s)
(param r 10mm)                       ; drum radius
(param v 0.1m/s)                     ; lift speed
(def rate (* (/ v r) 1rad))                               ; the drum's turning rate: a speed over a radius is radians a second
(def turned (trapezoid-position t (/ rate 0.3s) rate (* (/ 100mm r) 1rad)))   ; its angle at t

(body drum (shape (min (move (cylinder-z 4 110) 0 0 60) (move (cylinder-z r 60) 0 0 60))) (attr material "aluminum-6061"))
; what the bearings sit in: two plates bored for them, on a back plate
(body housing (shape (subtract (min (move (box 60 60 7) 0 0 38) (move (box 60 60 7) 0 0 118) (move (box 60 8 87) 0 34 78))
                               (cylinder-z 11 300)))
  (attr material "aluminum-6061"))
(assembly (hoist (angle 0))
  (instance housing housing)
  (instance motor (m.drive angle) (at (move 0 0 -2)))
  (instance ba b.bearing (at (move 0 0 38)) (fact required-life 5000h))
  (instance bb b.bearing (at (move 0 0 118)) (fact required-life 5000h))
  (instance drum drum (frame-of motor shaft))   ; on the motor's shaft: what the motor turns rides its frame
  (supports motor/rotor ba bb))                 ; what holds it: the two bearings
(assembly root (instance h (hoist turned)))
(load cord root/h/drum (force 0 -121.7N 0) (at (move r 0 78)))   ; the cord's pull, tangent to the drum
(require hoist-ok (within-ratings root/h))
```

Three things carry the answer, and each is said where the machine is
written:

- **What a motor turns rides its shaft.** The drum is placed on the
  motor's frame, `(frame-of motor shaft)`, so the rotor's joint carries
  it and its torque is the motor's. A drum on a joint of its own that
  reads the same angle moves with the rotor but is not on it: nothing
  says what joins them, and the load between them is refused, saying so.
- **What holds a joint's load is named:** `(supports PATH S ...)`, in
  the assembly with its own instances' paths (or at the top level, from
  root: `(supports root/h/motor/rotor root/h/ba root/h/bb)`). Each S is
  then given its share of PATH's reaction by statics (Supports), and
  PATH's own part carries none of it. A bearing nothing names has no
  load from the model, and its rating says to name it. A part on a
  joint needs no line: the joint holds it.
- **What moves is driven, and what is loaded is held.** A body on a
  joint of its own, moved by its own formula, is joined to nothing: the
  motor that should move it is charged nothing for it (`ink effort`'s
  `unpowered`), and `ink require` says so and is not all yes. Write its
  motion as a formula of what moves it (`(slide y (* r (/ angle 1rad)))`
  for what a drum winds) and relate the two by what is between them
  (`(relate drum weight (cord 10mm))`), or place it on that part's frame.
  The same is said of a part on a joint that does not move and is no
  actuator's, loaded along it, and of a load on a part on no joint.
- **What passes the force is in the model.** A joint moved by an
  actuator's coordinate, not riding its moving frame, that no relation
  reaches (a belt, a cord, a mesh) has nothing between it and the motor:
  the motor is charged for what it moves, but the part that would carry
  that force (its mass, tension, stiffness, losses and its own rating) is
  not checked. `ink require` names the joint and is not all yes. A joint
  that reads the coordinate only through a loop's closure is not named:
  the links that close the loop are the parts between.
- **The loads applied are loads:** `(load ...)`, `(gravity ...)`.

A belt axis, the machine most often built first: a stepper turns a
20-tooth GT2 pulley, the belt runs round an idler of the same size and is
clamped to a carriage, which slides. The belt is a `belt-span` a span at a
time (std), each span the line it runs along, from where it leaves one
part to where it meets the next, in the assembly's coordinates; the
carriage's slide is the motor's angle times the pulley's pitch radius.
(MOTOR.ink stands for the 17HS15-1504S stepper's import path:
`ink parts 17HS15-1504S` gives it where a parts catalog is installed.)

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

It answers yes. The motor's peak torque is 1.58 N mm: by hand 0.5 kg at
the motion's peak acceleration, 150 mm x (pi/2 per s)^2 = 0.37 m/s², is
0.185 N, 1.18 N mm on the 6.37 mm radius, and the rest turns the pulleys
and pays the belt's loss. Its answer also says what a yes here does not
cover: the pulley, idler and carriage carry no rating, and the motor's
ratings do not read the belt's pull across its shaft.

What the solvers then offer as facts, in the order `ink require` gives
them (after anything the design states):

- On the part a joint moves (its path, or a frame the part has: a
  motor's `shaft`): `torque` and `speed` (`force` and `speed` on a slide)
  where the joint is an actuator; `radial-load` and `axial-load`, the
  joint's reaction across and along its axis (Reactions); `moment-x`,
  `moment-y` and `moment-z` (N mm), that reaction's moment about each
  axis of the instance the joint moves (a rail's carriage: its roll,
  pitch and yaw, by its part's frame).
- On a coupling's part (an end of one, or a part within an end):
  `tension` (N), the tightest of its couplings' lines there, as the
  balances give them (Tensions): a screw's thrust on its nut, a belt's
  tight span on its pulleys, a tooth's force on meshing gears.
- On each support: `radial-load` and `axial-load`, its share's force
  across and along its load frame's z (a guide's carriage: across its
  x, its travel, with `moment-x`, `moment-y`, `moment-z` about its
  axes); `speed`, PATH's own.
- `speed` on the part of any joint, an actuator or not: the rate of its
  drive as written, which the design's formulas fix (an idler a belt
  turns, a carriage it pulls). So a speed typed under it is held beside
  it, as any typed number is.
- Each from the least it needs. A `speed` is the written motion's: it
  is answered with no mass and no gravity in the design, and solves
  nothing. A torque, a force or a load is answered where the moving
  bodies it rests on have a mass and a gravity is declared. A moving
  body with no mass leaves open ("the design lacks something", naming
  it) only what it would have changed: the efforts of the actuators it
  moves with, and the reactions that carry it. The other actuators'
  efforts are the ones the design would get without that body.
- Thrust on radial bearings: where the supports are radial bearings and
  some of them carry thrust too (a deep-groove ball bearing: its
  record's `thrust`), their radial shares are given whatever the thrust,
  and the thrust is theirs between them. Statics does not say which one
  locates the shaft, so each one's `axial-load` is between 0 and the
  whole of it, less what the design states of the others: state it,
  `(fact axial-load 20N)` on the one that locates the shaft and `0N` on
  one free to float, and each is exact. A bearing whose record carries
  no thrust takes none; with no bearing that can, the split is refused,
  saying to name a thrust bearing.

Each is read at every moment of the motion (`ink effort`'s: the efforts'
steps), as it is there, and the rating's answer is the worst moment's. A
load is an estimate within what the masses it rests on can move it, each
body's mass, centre and inertia, and each carried body's mass per mm,
anywhere in its proven bounds. A reaction's peak between those moments is left out,
and each fact says so. A moment across a joint's axis is left out of its part's
`radial-load`, saying how much: a shaft's rating is of a force at a
point. Where a coupling (a belt, a mesh, a screw) pulls on what a joint
carries, the joint holds that too, and its reaction has it: a carriage
a belt is clamped to, a motor's shaft under its pulley. Where the reaction or its
split is refused (Reactions, Supports), so is the fact, for that reason.

A typed number never makes an answer better than the model's. A stated
fact is held beside what the model works out about the same thing (a
load, a torque, a speed; a study's stress or deflection), by size:

- stated at or past the model's most (its band's high end), it stands:
  an allowance is the designer's to state;
- stated short of the model's least, the two are sources that disagree:
  `not settled: sources disagree`, each shown with its number and where
  it is from;
- stated inside the model's band, it does not narrow it toward the
  better: the answer is worked with the model's band, from the worse of
  the two least to the model's most, and says what was typed and where.
  (A bearing rated 1370 N whose load the model works out as 1352 to
  1408 N is not settled, typed 1360 N or not.)

Leave the fact out to be answered from the model alone.

A bearing's life is of its duty, not of a moment: over a motion it is
worked from ISO 281's equivalent load of a varying duty (`duty-load`:
the mean of the equivalent load to the life exponent, weighed by the
turns made at each moment) and the mean speed (`mean-speed`), each read
at every moment and weighed by its step's time. Where load and speed are
both stated they are the duty. At a still pose there is no duty, and a
life is not available, saying so: ask over a motion that turns, or state
the speed the bearing will run at.

`ink check` lists what the design applies and how it says it is held:
`gravity`, `loads` (each with its force and torque at the moment asked)
and `supports` (each divided reaction and its parts).

**In an expression.** std's ratios hold a carriage's share or a belt's
tension against a record written in the model, each the load over the
rating, so a static safety factor s holds where the ratio is at most 1/s.
Each takes the part's whole record (the one holding `load-ratings`, or a
belt's `belt` section), never a section of it. From
examples/slider-loads.ink:

```lisp
(def records (read-json "slider-loads.json"))
(def c (nth (reaction root/s/c (by root/s/c)) 0))
(require carriage-static (value (carriage-static-ratio c (get records "carriage")) <= 0.5))
(require belt-tension (value (belt-tension-ratio (belt-tension (nth (reaction root/s/c) 0) 30N) (get records "belt")) <= 1))
```

A bearing's static load, life and limiting speed are carried as ratings
on the bearing body and checked with `(within-ratings ...)` (as in
examples/slider-loads.ink and Ratings above).

- `(carriage-static-ratio W RECORD)`: the largest of each force across its
  travel over C0 and each moment over MR, MP and MY (HIWIN's fSL and fSM,
  each on its own).
- `(belt-tension f t0)`: a belt's tight side, installed at t0 (a force), f
  the force it passes (`(nth (reaction root/d/f) 5)` over the pulley's pitch radius, `10mm`,
  or `(nth (reaction root/axis/c) 0)` along a slide): max(t0 + |f| / 2,
  |f|) installed and locked; |f| + t0 with a constant-force tensioner. A
  belt related as a coupling gives its spans' tensions themselves,
  `(tension A B)`.
- `(belt-tension-ratio T RECORD)`: T over the belt's working tension.

A motor's continuous rating is against its RMS effort over a duty cycle:
`(effort-rms PATH)` is the RMS over the file's `(t-span A B)` of the
effort actuator PATH supplies, `ink effort`'s `rms`: the efforts' steps'
middles (Steps) by the midpoint rule, an estimate; where groups are read
round, its band is the most that moves the effort at a step. It is
refused where a step is `unchecked` (an impulse).

```lisp
(param rated-torque 0.98N*m)   ; the gearmotor's continuous rated torque, from its record
(require cool (value (effort-rms root/arm/g/output) < rated-torque))
```

Each reads its rating through the record's ratings vocabulary; one the
record lacks is an error naming it (`no key 'C0-N'`, with the keys there
are). No typical value stands in. They are estimates: the shares rest on
masses measured to 2%, and the ratings are the maker's.

### Tipping

`(zmp R)` (std), R `(reaction root FRAME)` with FRAME on the support's
plane and its z up into the model, is where on that plane the support's
force acts, `(vec x y)` in FRAME (x = −my / fz, y = mx / fz); where the
support does not push (fz ≤ 0: the model lifts off) it is at infinity.
The model stays on its support while that point is inside its footprint,
a field of x and y negative inside, read with `warp`:

```lisp
(def foot (max (- (abs x) 70) (- (abs y) 40)))   ; 140 x 80, centred
(def p (zmp (reaction root (at))))                ; the floor: z = 0, up
(require stands (value (warp foot (nth p 0) (nth p 1) 0) < 0))
```

Its value is how far the point is outside the footprint in mm (negative:
inside); failing, its moment. examples/servo-arm-effort.ink: a slow
swing stands, a fast one tips its stand.

It takes the support as holding by friction; whether it does is std's
`(slip R)`, R in the same frame: |(fx, fy)| / fz, infinite where the
support does not push. The model slides where that passes its friction
coefficient, which the claim states:

```lisp
(require grips (value (slip (reaction root/cart/wheel (frame-of root/cart/wheel contact))) < 0.8))
```

It is one contact's: root's reaction sums every contact, whose ratio is
no single contact's (a rear-wheel-drive cart at a/g = 0.5 reads 0.5 there
while its driven wheels may slide). A base on several contacts needs each
one's.

### ink effort

`ink effort FILE [--target A] [--t T | --t-span A:B]` answers at a moment
(`--t`), or over the span (`--t-span`, else the file's `(t-span ...)`;
neither: at t = 0). `--modes` (Modes) gives its own answer.

- `actuators`: each one's `path`, `effort`, `rate`, `lo`, `hi`, `margin`
  (min(effort − lo, hi − effort)) and `reflected_inertia` (g mm² a turn, g
  a slide: the effort per unit of its own acceleration); what its effort
  `leaves_out`; the couplings whose loss it leaves out for want of an
  efficiency (`loss_unknown`); the most its couplings' loss can be off by
  (`loss_bound`); why its claim would be refused (`refused`, with no
  rating margin in that case); and its `basis`: `rigid` (inverse dynamics), `command` (it moves simulated
  coordinates alone: its command, clamped), or `spring` where the file's
  simulation releases its joint (Simulation). A `spring` actuator's effort
  is its spring's force, as its claims read it, its limits at the
  commanded speed; the rigid model's demand for the written motion is
  under `rigid`, `method` says so, and over a span its `band` is its
  stiffness times the following error's band (`band_kind` `measured`, no
  bound). Over a span it is read at the efforts' steps (Steps) and where
  the spring's force turns back, with its own `loose` and `partial`.
- `coordinates`: each one's `name`, `value`, `rate`, `acceleration` and
  `effort`, per unit of the argument as written (N mm per degree for an
  angle in degrees; N for a length in mm).
- `reactions`: root's (the world's axes about its origin: its support's)
  and each joint's under it (in its joint's frame): `path`, `force` (N),
  `moment` (N mm); the bodies it carries `without_mass` or `unbounded`, or
  whose shape moves (`leaves_out`); `moves_with`; where its claims would
  be refused, force and moment null and `refused` saying why.
- `couplings`: each relation's `tensions`, or why they are `refused`, and
  what they leave `unseparated`.
- `unpowered` (each coordinate's effort nothing supplies),
  `driven_by_nothing` (bodies placed by t other than through a
  coordinate: their motion is charged to no actuator), `carried` (with
  each one's mass), `deforming` (with mass, left out).
- `simulated`, where the file simulates (Simulation).
- `method`: how it was found; `accuracy` null (an estimate); `masses`:
  the bodies measured, those weighed from a published mass, those with a
  `spin-inertia`, the cells measured at, those `unsettled` (not measured
  to 2%).

Over a span, each actuator's `peak` (effort, moment, rate, reflected
inertia), `rms` and least `margin` (moment, effort, rate, limits); each
coordinate's `peak` (value, rate, acceleration); each unpowered
coordinate's `peak`; each reaction's `peak_force` and `peak_moment` (the
moment its force, or moment, is largest, with both there), or where it is
first refused, `refused` {t, why}; each coupling where its tightest line
peaks. Where groups are read with their spin at rest (Steps), each
actuator, coordinate and unpowered coordinate gives `bound`, and each
reaction `bound` {force, moment}: the most their departure from round
moves it over the span at any angle of their spins (absent where 0).
`motion` gives the `span`, the `steps`, how they were chosen (`basis`),
the intervals `unchecked` and `loose` (adjacent steps joined), and
`partial`. Every step counts but an unchecked one, where nothing is read
(a simulation that stops short leaves the span past its stop unchecked):
a loose one is read at its middle, coarser than the rule. The peaks and
margins are of the moments counted, so at least that much where some are
not; the RMS is of the same moments, each weighing its step, and null
where a step is unchecked (a mean needs the whole span). `partial` is
true where any step is unchecked or loose.

Refused: a moving body without mass (named: no density is guessed); an
actuator whose drive reads no coordinate. An actuator's effort rests on
the coordinates it moves and those coupled with them; where a deforming
body with mass moves with them (`leaves_out`), or one of them
is unpowered at a moment, its claim is `unknown` there, saying which. An
effort is found at one value of each param: over a toleranced param it
reads (a moving body's size or density, an actuator's limits), its claim
is `unknown`. An actuator held at a pose a body's shape follows leaves
that body out, and its claim is `unknown` saying so.

### Modes

`ink effort FILE --modes [--t T] [--target A]` gives a mechanism's natural
frequencies at the pose its coordinates have at t (default 0), lowest
first, and how it moves in each: small undamped oscillations about the
pose at rest, as a multibody model is linearised (Adams/Linear, Simscape's
linearize). It needs no gravity. A simulation steps the same springs
forward (Simulation).

Each coupling line stating a stiffness K (a belt's spans, given
`(tensile-stiffness EA)`) acts as a spring along the line while under
tension. With zero pretension, equal spans with opposite stretch per unit
of each coordinate contribute one spring between them. The stiffness and
each corresponding stretch coefficient may differ by floating-point rounding. An unmatched line
that can stretch at this pose leaves the linear modes undetermined.
A joint stating `(stiffness K)` has a spring on its value. K is a number
or a field of params, in N/mm for a slide and N mm per radian for a turn
(`2N*m` is 2000). Other units are an error. A catalog stepper states its
rotor's holding stiffness N_r T_h (its rotor's teeth, 90 / its step angle, times its record's holding
torque: the slope of T_h sin(N_r angle) at rest, two phases; Acarnley,
Stepping Motors). Each joint a stiff line's end rests on (and every one
supporting it, up to the target) is released: it gets a coordinate of its
own, carrying what it carries (a belt stretches between a pulley and its
carriage, round idlers turning on their own). So is each joint its own
stiffness holds that no coordinate drives. A joint a coordinate drives
keeps that coordinate, unless another joint stating a stiffness moves
along it too (each released, each its own spring). A line stating no
stiffness is rigid: one reaching a released joint releases the joint at
its other end and holds its stretch at 0 (a rack adds J / r² to its
carriage). The mass matrix is every rigid body's measured inertia along
the coordinates and released joints; the stiffness matrix each spring's
stretch per unit of each, squared, times its K. Directions moving no body
are dropped; those with stiffness and no inertia (an idler given
`(attr spin-inertia 0)`; an inertia within a billionth of what its
bodies' masses and inertias sum to there is none, as a simulation reads
it) are condensed out (Guyan); then K φ = ω² M φ.
A mode is held where some spring stretches along it, each spring judged
on its own (a stretch under a millionth of its stretch per unit of each
coordinate times the mode's motion along it is none), and its ω² is its
springs' stretch squared times their stiffness over its inertia: a stiff
spring elsewhere hides no soft one (a 1 N/mm slider beside a 1e15 N/mm one
keeps its 5.03 Hz). A pose held still answers as one a motion passes.

- `modes`: each one's `frequency` (Hz; 0 where nothing holds it) and
  `moves`, the joints whose bodies carry its kinetic energy, most first:
  `joint`, `energy` (its share, at least 0.001) and `motion` (deg a turn,
  mm a slide, the first scaled to 1); a body no joint moves gives `body`
  and `energy`.
- `pose`, `springs` (each line's and joint's `stiffness`, N/mm, or N mm
  `per` radian, and its `site`, the file and line stating it),
  `released`, `rigid` (lines held as constraints), `carried` and
  `deforming` (left out, with their masses), `method`, `masses`,
  `accuracy` null.

Refused: nothing under the target stating a stiffness; a stiffness below
0 or not finite; a coupling's numbers an effort refuses (a tensile
stiffness written as a number after the efficiency, where EB goes, is
named); a rigid line that separates its parts; a body a coordinate moves
otherwise than through the released joints, where that coordinate drives
one of them (tied by the drives' expressions alone: relate it to what
drives it), naming it.

It is an estimate, not a bound: what it leaves out (damping, the
motion's speed, gravity's and the tensions' own stiffness, the frame's
and parts' flex, tooth and mesh stiffness and backlash, the belts' own
mass, deforming bodies) mostly lowers the frequencies, but the stated
stiffnesses err either way. No catalog belt record carries an EA (the
makers of the catalog's GT2 belts publish none); the model states its own,
as each spring's `site` shows. Gates publishes its 2MGT fibreglass-cord
belting's (Light Power & Precision manual, Table 6: LL2MGT06, 3.8 lb for
0.1%, EA 16903 N at 6 mm); other GT2 belts are measured -12% to +43% from
it, a mode the belts hold -6% to +20%. A stepper's holding stiffness is at
the current and phases its holding torque was measured with, which the
records do not state: a microstepping driver at the rated current gives
1/√2 to √2 of it, a mode the steppers hold 0.84 to 1.19 times ink's
figure. Stiffness no record states (a steel-cord belt's, a chain's, a
mesh's, a servo loop's) is not found. An abrupt change of acceleration a
leaves a mode ringing by a / ω² about where it settles; a change spread
evenly over whole periods leaves none (an input shaper).

### Simulation

A coordinate whose motion is not written but follows from the forces on
the mechanism (an arm dropping when power is cut, an arm steered by a
controller, a carriage at full voltage) is passed as `(simulate X0
[V0])`: its position X0 and speed V0 (default 0; fields of params, in the
argument's units) at the start of the file's `(t-span A B)`, from which
the mechanism is stepped forward over the span. It is the argument itself
of an instance root places, `(instance arm (arm (simulate -90)))` or
through a def; every coordinate so passed is simulated together, the
others follow their written motion. A file that says `(springs)` in
its run simulates its springs too (below), for every claim it makes and
`ink effort` alike; one that does not simulates its coordinates alone,
its couplings and stiffnesses as an effort reads them, a stiffness no
spring. What a file reads does not change its physics: a
`(following-error PATH)` is read only where it says `(springs)`, and is
refused otherwise, saying so.

The state model measures each moving body's mass and inertia once.
Its density may vary across the body in its own frame. A density that
depends on time or a moving coordinate makes the simulation unavailable,
including claims on the simulated coordinate.

```lisp
(def swing (simulate 0))
(assembly root (instance a (arm swing)))
(command root/a/s/output (- (* 40 (- 45 swing)) (* 1.2 (diff swing t))))   ; a PD to 45 deg
(command root/m/rotor 1e9)                                                ; full effort
```

- **Forces**: gravity, the loads, and each actuator moving a simulated
  coordinate, which supplies its `(command PATH EFFORT)` (N mm a turn, N
  a slide) clamped to its effort limits at that state: a command past
  them is full effort. EFFORT is a field of params, t and the simulated
  coordinates. An actuator with no command supplies none (power cut). A
  command on no actuator under root, or on one moving no simulated
  coordinate, is an error naming it; so is an actuator moving both a
  simulated and a written coordinate. A command may not read a reaction,
  a tension, an RMS effort or a following error.
- **How**: at each state the model's inverse dynamics and generalized
  inertia give the simulated coordinates' accelerations: the equations of
  motion in the mechanism's own coordinates, constraints embedded (Kane's
  method). Where a coupling's loss turns with the power's direction, they
  are found by Newton's method. The integrator is Dormand and
  Prince's 5(4) pair, switching to Radau IIA of order 5 where its
  stiffness test finds the motion stiff and back once it is not; a
  simulation of springs uses Radau IIA throughout. `simulated.method`
  says which. Positions are held to
  1e-8 a step. A switch in t (a `select` or `floor` on t in a command, a
  limit, a load or a written motion) is found by bisection and stepped
  across, so a pulse is never stepped over. Between steps the motion is
  the quintic Hermite through each end's position, speed and
  acceleration. It is stepped only as far as what reads it asks, once for
  every load asking the same question. A simulated coordinate and its speed
  share their endpoint values between steps. Speed and effort claims may
  therefore use the integrated path as an estimate. They check the written
  motion and speeds used by its commands from the simulation's start,
  including for a question about a later single moment.
- **Springs**: a file saying `(springs)` releases the joints its
  springs hold (as Modes does, but an actuator that is commanded or
  states no stiffness): each one's deviation from what the coordinates
  command (a small turn, radians, or slide, mm) is a state beside the
  simulated coordinates. A stiff line pulls with its pretension plus K
  times its stretch, never below 0 (a span goes slack and pulls again once
  taut); a joint's stiffness holds its deviation at 0, at its drive's
  value (a stepper's rotor at its commanded angle). A joint carrying no
  inertia along its own motion is held where the springs leave it. The
  deviations start where the springs hold the mechanism at rest at the
  span's start; where the written motion jumps, each keeps its place and
  speed. A springs-only position claim can read the resulting jump in
  following error. A speed claim across that jump is unknown. Speeds
  used in the command itself must still be defined.
  `(following-error PATH)` is joint PATH's deviation at t (mm a
  slide, deg a turn; an error where PATH is no released joint or carries
  no inertia): `(require lag (value (abs (following-error root/xy/gantry))
  < 0.11))`. It is an estimate known either way to within its band: 100
  of its tolerances (each deviation held to 1e-6 of the largest so far at
  the parts it moves, at least a nanometre), and 0.3 more for each period
  of the springs' fastest mode stepped. The band is measured, no bound
  (DYNAMICS-VALIDATION.md); a claim within it, and its masses' band
  (below), of its bound is undecided.
  An actuator whose joint is released supplies its spring's force, its
  stiffness times its following error, which its effort-limits and
  effort-rms claims and `ink effort` read (`basis` `spring`), within its
  stiffness times the band at the span's end and its masses' band (at
  the moment; over the span, its largest). An actuator whose
  coordinates move what the springs release has its effort claims
  `unknown`, as have a reaction carrying a part they move and a sprung
  coupling's tension, each saying why.
- **What reads it** reads an ordinary motion, as an estimate
  (`solved_inputs` naming the simulated coordinates, `notes` its method):
  claims over the span are measured at the middles of 32 even steps and
  where each simulated coordinate they read turns back or its speed reaches
  an extreme of the interpolated motion. The minimum and maximum of each
coordinate’s acceleration, jerk and snap are also measured. These include
interior polynomial extrema and both sides of interpolation knots. At
launch its jerk and snap come
  from the first accepted step. An actuator moving only simulated
  coordinates supplies its clamped command, within its limits by
  construction, so its effort-limits claim nothing. What no actuator
  supplies along a simulated coordinate is never `unpowered`.
- **Masses**: each body's mass, centre and inertia are measured, each of
  its ten moments known to its proven bounds (Scope and limits), and the
  motion is known to within what those move it. Each moment alone at the
  far side of its bounds moves the motion by its first-order sensitivity:
  the variational equations stepped beside the motion (the change it
  makes to the accelerations at each state, by inverse dynamics with that
  moment moved, and the motion's Jacobian; the exponential midpoint rule
  on its steps; a released joint's rest moved by what the springs give
  for it), read between steps by its own quintic Hermite. A simulated
  coordinate, a following error and their rates to the fourth are each
  known either way to within the sum of those sizes at the moment read
  (moments moving it alike are carried as one): a claim within it of its
  bound is not settled, one past it decided. At an extreme a moment that
  only shifts the motion in time moves it by nothing. A spinning body
  read round within what its measure is known to carries that departure
  at any angle: its centre off its axis each way across it, and its
  second moment not round about it four ways. A change within the
  rounding of what it moves is none (a block falling free moves by none
  of its mass). A body's own published mass (`(attr mass M)` on it, read
  at the configuration it is published for) is its mass exactly: its
  moments are the measured ones scaled to it, the scale known as far as
  its volume is (a curved body so published moves as its mass does). A
  share of an assembly's published mass, or of a definition's a body is
  built from, or one read at another configuration, is its density times
  its measured volume, each known to its bounds: the density's over the
  volumes it is spread by, and the body's own. A carried body's mass per
  mm is a source too, at what its profile's measure (and its density) is
  known to. Against reruns with each moment at either side of its bounds,
  on a double pendulum and a ringing spring, the band is within 13% where it
  is a twentieth of its largest (validation.rs's
  `the_masses_band_is_what_each_measured_moment_moves_the_motion_by`).
  First order, no bound: a rerun at one side reads its second-order part
  too, up to 1.8 times the band where a moment's change turns at an
  extreme of a speed.
  A claim on an effort, a reaction or a tension in a file whose simulated
  coordinate has a band is not settled, the band read through its own
  solve, which carries none (every such solve reads every coordinate).
- `ink effort`'s `simulated`: the simulated `coordinates`; how far it was
  stepped (`end`), why it stopped short (`stopped` {t, why}, null where it
  ran its span), its `steps`, `rejected` steps and `evaluations`, its
  `method`, `accuracy` null; and, with springs, `springs`: `joints`, each
  released joint's largest following error (`peak`, `at`), its value at
  the `end`, its spring's largest force (`effort`) and, once the written
  motion is at rest, its `ringing` (`from` then, `amplitude` half its
  swing, the `frequency` of the mode carrying most of it and that mode's
  `share`); the joints `held` without inertia; `spans`, each spring line's
  `least` and `most` tension with when (0 is slack); `still`, when the
  written motion comes to rest for good. Where it does not start,
  `simulated.refused` says why, the efforts standing.
- Only root's mechanism is simulated: `ink effort --target A` on another
  assembly reads its coordinates as written.

It stops, the motion read up to that moment and refused past it, saying
when and why, where the mechanism cannot go on: a dead point of a loop
(the loop does not close near its last root, or the simulated
coordinates' generalized inertia has grown 10⁸-fold); a drive past its
curve's last speed; no inertia along a simulated coordinate pushed along
it (one within a billionth of what its bodies' masses and inertias sum to
there is none, as the modes read it: an idler given `(attr spin-inertia
0)` turned by it alone); accelerations not found; a command or state not finite; past 20000
stages evaluated; friction holding the mechanism still (at rest, no
acceleration agrees with its couplings' loss: a PD-held lift on a lossy
screw, a self-locking screw holding its load); a self-locking coupling
driven back while moving (a negative effective inertia: write the motion
and read its effort instead); or where inverse dynamics refuses the
forces (a coupling with no efficiency passing power along a simulated
coordinate, tensions not found, a body whose shape moves). Refused at the
start: no rest where the springs can hold the mechanism (naming the
joint); a motion of the released joints that nothing resists (naming its
joints); a rigid line between joints the springs release (naming it); a
body tied to a released joint by the drives' expressions alone; a
simulated coordinate moving nothing but released joints.

### Steps

Over a span an effort is measured at the span's two ends and at the
middles of steps chosen from the motion, and its claims at the same
moments. Each moment an actuator's limit may jump (a select or floor on
t: a drive disabled for 0.2 ms) is found by bisection to the resolution of
t (adjacent floats) and given a step of its own, so two jumps apart in t
are found apart and the moments between them measured. A step a limit may
jump in (one of those, or one holding jumps past the 512 the search
locates) is judged by the limit's range over the whole step, bounded by
intervals: where that can cross the effort read at its middle, the step
is `loose`. A smooth limit given directly as a field of time is also
checked over the step. A range that can cross the measured effort leads
to finer rigid-effort steps and leaves the step `loose` if unresolved.
Steps are halved until, inside each, no coordinate's speed or
acceleration, no moving body's rotation, translation's acceleration or
partial velocity, no carried body's circulation and no load's force or
torque can stray from its value at the middle by more than 1% of its peak
(max|f″| h²/8, from bounds of rates over the step). A partial velocity
that nearly cancels can use higher time derivatives to enclose its
curvature more tightly. Each component still uses its own measured peak.
Steps are then halved until each
effort lies within 4% of its peak of the line through its neighbours, a
difference under its rounding (1e-9 of the terms summed into it) being
none. Then, where an actuator's effort or its margin to either limit
turns back between the moments measured, within 8% of its peak of its
extreme there (a peak the line test reads up to a few percent low), the
moment it turns back is found by bisection on the sign of its rate (by
differences), as a simulated coordinate's extremes are, and given a step
of its own: the peak and least margin are read there. A step where a
coordinate's position or speed may jump is `unchecked`, and its answer
`unknown`. When a motion switches on a linear clock, the position and
speed can also be shown continuous by comparing their values exactly at
the switch, including when several terms switch together or a cycloid
stops at an exact full turn. Where only acceleration jumps (a trapezoid's
corner), the step is halved to the finest and measured either side.
At most 4096 adaptive steps: those left coarser are `loose`.
A claim whose steps include a loose or unchecked one is
`unknown` unless a moment measured fails it.

A spinning body's effort reads the spin's rate and acceleration, and its
angle only through its departure from round about the spin's axis. A body
spins with a link where the lowest link of its placement that moves turns
it: a joint; a WHERE riding one, placed on a frame of what the joint turns
(a pulley on a motor's `(frame-of m shaft)`); or a WHERE that is a
rotation about the child's own axis. A body whose shape and density show a
turn symmetry about the axis (`turn_symmetry`: a roller, a pulley, a gear,
a hex shaft), or whose axis is its own z and that has a `spin-inertia`, is
round: its steps follow its placement with the spin at rest. The other
bodies spinning with one link (a D-cut rotor and the pulley and grub
riding it) are one group, whose departure from round is their measured
mass off the axis and second moment a turn moves. Each body's measure
knows its moments only so far (twice what its last refinement moved
them): a part of its departure within that is its measure's error, not
an imbalance, and the body is measured round there, at every moment and
angle. A body whose shape repeats about a line of its own along the spin
(a wheel moved a hair off its shaft; a pulley but for its set screw's
hole and its bore's flat) is not left to its measure: the part that
repeats has its centre on its own line exactly, and no product of
inertia with it, so what it adds is exact (a 308 kg wheel 2e-7 mm off
its axis, measured 0.5 um off, reads its real 0.024 N at 100 rev/s, and
nothing of it is unknown), and the small parts that break the symmetry
(what is cut from it, what is added) are measured, each known only so
far: the body's unknown is theirs alone. What an imbalance is under what
is known is not
known, and is not read as none: it is carried as a departure that may
lie at any angle, by the same bounds a real one read at rest has, at
every moment, whether the group's spin is resolved or not. An effort,
a reaction and the ratings that read the reaction are each known only to
what an imbalance that large would move them at the body's spin (a body
with no line it repeats about, whose measure is known to 31,498 g mm,
leaves 12,435 N on what carries it at 100 rev/s), so a claim inside that
is not settled and one past it holds. Where the
motion has no bound over a step (simulated, or through a solve) the
bound is the moment's own. `motion.basis` names each such body with
what its moments are known to. A tension is not so bounded: a group
spinning under a balancing joint has its band left out of the tensions.
Its effect on each
effort, B, is bounded over each step; a B no bound is found for (a group
moving through a solve) is infinite. A group is read with its spin at rest
where the B of the groups so read, summed, stay within 1% of each effort's
peak; or where resolving its spin is past the step budget and the B
summed stay under each claim's margin at every moment measured (a fan
held at speed, far from its limits). Else its rotation is resolved as any
body's. A group read at rest is read round, so the effort read is the one
it would need at every angle, and the real one within B of it.
`motion.basis` names each group with its departure and the most it moves
each effort. A spinning body riding a WHERE's own turn, rather than a
joint, is resolved on its own.

**An effort claim's value** at a step is `min(effort - LO, HI -
effort)`, read with the groups read at rest made round (outside the
file's span, and in a file without one, every body as measured). With B
the most those groups move that effort there, the step holds where the
value is at least B, fails where it is below −B or where the effort at
that moment's own angles is past a limit (a witness), and is undecided
between: the claim is `unknown` there, saying the effort, the limits and
B. An unpowered coordinate a claim rests on is likewise one where the
effort nothing supplies is past B, undecided under it. The claim's
`motion.basis` gives the most B over its steps. A claim failing gives the
first failing step's moment and terms; holding, its value nearest its
bound over the steps. Between steps it is an estimate, as every effort
is (`"accuracy": null`, its `solved_inputs` naming its effort), unless
it is bounded there (below).

**Between its moments** an effort claim on a mechanism the bound covers
is bounded, not estimated. Over each step its effort is inverse dynamics
in interval arithmetic rounded outward, over every moment of the step at
once. Each moving body's share is the virtual work of its mass moving
with its placement x = R u + p: along coordinate k, the integral over its
mass of (dR_k u + dp_k) . (R″ u + p″ − g), linear in its mass, first
moment and second moment, so it reads only its placement's second rates
and its partials along the coordinates; each load acts through its
point's partial velocities, as at a moment. Every number is read over the
step as its bound there, its bound at the step's middle and its rate's
bound over the step, and each sum is the meet of its plain bound and its
mean value form, which shrinks with the step's square. Each body's
measured mass, centre and inertia lie anywhere in their proven bounds;
the drives' constant ratios are inverted in intervals; the limits are
bounded over the step, so a torque-speed line is read over every speed
the step passes. A step whose bound lies inside the limits throughout is
proven. One whose bound crosses a limit is halved while that may settle
it, at most 4096 halvings in all, and where the claim fails at a halved
step's middle, past its masses' band, that moment is its witness and the
halved steps are read (a 5 N bump 0.2 ms wide between 16 even moments, a
dip 2 µs wide in a rating that follows the speed). A step not settled
(its bound still crossing a limit at the budget, as a body spinning too
fast for a step's bound of its rotation leaves it, or crossing at its
middle within what the masses or the rounding are known to) keeps its
moments' reading, as above: an estimate there. `motion.basis` says how
many steps were proven, how many halvings it took, and how many were
left to their moments. The bound covers a rigid mechanism whose motion is
written as a field of t, its drives' ratios constant and as many drives
as coordinates in each set they join. It is not made, and the steps are
read as above, where the motion is simulated or sprung, moves through a
solve made again at each moment (a loop's closure), has couplings or a
carried body's material running round its path, where drives share a
coordinate (their split rests on their limits) or a drive's ratio moves,
where a spinning body's departure from round is not known within its
measure, or where a rating is one number; `motion.basis` says which.

### Scope and limits

What these answers are:

- Efforts, reactions, shares, tensions and RMS are measured at moments
  (Steps): estimates, not bounds over the span between them, but for an
  effort claim the bound covers (Steps, between its moments), which is
  bounded over every moment of each step it settles. The masses
  are measured to 2%, and their centres and inertias to the same share
  of their own scale, each proven. An effort, and an RMS effort, carries
  how far its bodies' measured mass, centre and inertia, and its carried
  bodies' measured mass per mm, can move it within those bounds: a claim
  within that of its limit is not settled. A
  simulated motion carries how far they move it, to first order
  (Simulation, Masses).
- Modes are a linearisation about a pose; a following error's band is
  measured, not a bound; a simulated motion's masses' band is first order,
  not a bound.

What is not modelled:

- Static friction: no sticking. A simulation stops where friction would
  hold the mechanism still, saying so.
- Contact and hard stops: a joint's limits are claims the motion may
  fail, not stops; nothing collides.
- Damping: none is stated, so a ringing does not decay (the integrator
  damps it slightly, within the band).
- Joints' friction; a coupling's loss where it gives no efficiency (its
  efforts, and its train's reactions and tensions, refused); the force a
  mesh's sliding friction puts along its line of centres, and what reacts
  a routed belt's friction through its other parts' supports.
- Steppers are linear about their holding stiffness: a torque past that
  line, and so a lost step, is not seen.
- Springs are linear and lossless: no tooth or mesh stiffness, no
  backlash, no sprung coupling's loss, no belts' own mass in modes or
  springs, no gravity's or tensions' stiffness.
- Rigid and sprung parts mixed: a rigid line between joints the springs
  release is refused, and an actuator whose coordinates move what the
  springs release has its effort refused.
- Bodies whose shape moves with no `carried` describing them
  (`deforming`) are left out of every effort.
- Actuators on one coordinate share by a fixed policy (The effort); shares
  statics leaves open are refused (Supports).
- A loop placed by closed-form formulas gives no loads at its pins: a
  reaction on its links is unknown. (No form closes a loop by solving it.)
- A loop placed by closed-form formulas takes no notice of a pose where
  it cannot close (`sqrt (max a 0)` keeps its formulas defined there, and
  efforts are given for a stretched link): claim its closure (the
  catalog's linkages' `coupler-closure`).
- A load averaged over a motion for life (a gearbox's average torque, a
  power mean weighted by speed and time; a bearing's mean equivalent load)
  has no value: a motor's RMS does (`effort-rms`), and a bearing's life is
  read at each moment's load.

Known limits of what is found:

- A single coupling's loss is on the power through its ends' own joints:
  where its ends' joints turn in frames moving apart and the power enters
  through the one with no joint (a planet's sun held still, its carrier
  driven), none is charged.
- A simulation runs to the moment asked, past the call's time limit;
  under `ink mcp` such a call is stopped at the limit with no answer.

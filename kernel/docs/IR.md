# ink language reference

A document is a sequence of top-level forms in S-expression syntax. `;` starts
a comment.

A code block here that begins `(import "std")` is a whole model: it loads
as written. Any other is an excerpt, read with std and the names the text
around it gives.

## Tour

A name in parentheses, like (Claims), is a section with the detail: `ink
spec --section Claims` reads it.

A design is `.ink` text files of forms in parentheses; a comment runs
from `;` to the line's end. `ink COMMAND FILE` runs a command on one and
prints JSON.

**Shapes.** A shape is a field: one function f(x, y, z) of position in mm,
negative inside, a regular solid (Fields): two blocks placed flush are
one solid, and a slit needs a width.
`(min a b)` is
the union of a and b, `(max a b)` their intersection and `(neg a)` a's
outside. Every expression is a solid, so no operation can fail.

std, the standard library (`(import "std")`, written in ink), gives:

- `(subtract a b)`: a with b cut away;
- primitives centred on the origin: `(box sx sy sz)`, `(cylinder-z r h)`,
  `(sphere r)`;
- 2D profiles such as `(rect-2d sx sy)` or `(polygon ...)`, made solid by
  `(extrude-z g z0 z1)` or, built in, `(revolve-z g)`;
- `(move f dx dy dz)`, `(rotate-z f deg)`;
- fillets, threads, gears, sheet metal and materials.

`ink std` lists its contents; `ink std --section NAME` reads a section.

**Names and numbers.** `(param name default)` is a dimension any command
can set (`--param hole=6.2`); `(param name "Ada")` a text, such as a name a
customer types (`--param name=Grace`; Text params). `(def name expr)` names
a value, `(fn name (args) expr)` a function. A param's name and an import's namespace are
plain names: letters, digits, `-` and `_`, starting with a letter.
Namespaces join with `.`.

A bare number is in mm, g, N, MPa, N mm, °C or s, or names its unit,
metric or imperial: `0.25in`, `50lbf`, `68degF` (Units). Angles are
degrees, but `sin`, `cos` and `tan` take radians.

**Bodies.** A body is one part: a shape, attributes and named frames.
`(attr material "steel")` gives `ink mass` its density and packages
(add-ons such as dynamics) its other properties. A frame is an origin and
axes, `(at (rot-x 180) (move 0 0 5))`: turn, then move. A frame where
parts touch sits on the touching face, z out of the part.

**Assemblies.** An assembly places parts frame to frame, not at
coordinates. `(instance b bolt (frame-of p top) seat)` puts the bolt's
frame `seat` on frame `top` of `p`, an instance placed before it, face to
face: z axes opposed, x along x. `(shape-of root/b)` is a placed part's
shape; `(shape-of root/b shank)` is the piece of it built as `(tag shank
...)`. A joint gives an instance one motion, `(turn z DEG)` or `(slide x
MM)`, in params and the time `t` in s. `(t-span 0 2s)` is the time that
claims about the motion cover.

**Claims.** `(require NAME CLAIM)` states a requirement: `disjoint`,
`clearance`, `touches`, `subset`, `value`, `volume`, `integral` and more
(Claims). `ink require FILE` answers each yes, no, not settled or
not available; under that, the prover's `status` is `holds` (met, to its
accuracy), `fails` (not met) or `unknown` (undecided, often with a hint
saying what would decide it). The accuracy is how far off an answer may be: 0.05
mm by default, a part's `(attr tolerance x)` or `(accuracy x)` after the
claim to change it; `null` marks
an `unknown` or an estimate. `res` is how finely ink looks, a cell size
(Resolution).

```lisp
(import "std")
(param hole 6.6)                                 ; M6 clearance, mm
(body plate (shape (subtract (box 40 40 5) (cylinder-z (/ hole 2) 10)))
  (attr material "aluminum-6061")
  (frame top (at (move 0 0 2.5))))               ; on its top face, z up out of it
(body bolt
  (shape (min (move (cylinder-z 5 4) 0 0 2)      ; the head, above its seat
              (tag shank (move (cylinder-z 3 12) 0 0 -6))))
  (attr material "steel")
  (frame seat (at (rot-x 180))))                 ; under the head, z down out of it
(assembly root
  (instance p plate)
  (instance b bolt (frame-of p top) seat))
(require seated (touches (shape-of root/p) (shape-of root/b)))
(require free-fit (clearance (shape-of root/p) (shape-of root/b shank) 0.2))
```

Both claims hold. With `--param hole=6.2`, free-fit fails: its witness is
two points 0.15 mm apart.

**The loop.**

1. Write the model.
2. `ink check`: what loaded, and errors and warnings that say how to fix
   it.
3. `ink render --out FILE.png`: four labelled views. Fix until it looks
   right.
4. Claim what the design depends on, and run `ink require`.

A wrong placement passes check, so claim where each part must be
(Placement by example).

`ink help` lists the commands, among them `interfere` (parts that
overlap, at rest or moving) and `mesh --out FILE.stl` (or .3mf, .glb). A parts catalog,
installed apart from ink, adds parts: `(import "parts/FAMILY/NAME.ink" as
b)` names its bodies `b.NAME` and its record `b.spec` (Top-level forms: Libraries).

## Units

A bare number is a length in mm (so an area in mm² and a volume in mm³),
a pure number, or, where a pose is taken (`rotate`, a turn's drive), an
angle in degrees. Every other quantity says its unit where its number is
written: `(attr mass 12g)`, `(force 5N)`, `(t-span 0 2s)`,
`(speed-limits -90deg/s 90deg/s)`. A bare number where one is read is an
error that names the unit to write; a bare 0 and an infinity stand, being
no quantity's. A value worked out says its unit by a product, `(* 1g (get
spec "mass-g"))`, or, as a temperature must (a point on a scale is not
multiplied into), by a zero of its unit added to it: `(+ 0degC (* 2 x))`,
`(+ 0deg/s rate)`. Beside a quantity that is not a length or an angle a
bare term is an error too (`(- t 300)`, `(+ p 0.2)` beside a pressure):
nothing says its unit, so say it (`(- t 300ms)`). So is a value whose unit
a function lost (`(exp 1g)`, the sine of a length): divide what goes in by
its unit first. A placement's numbers are lengths and angles: `(move 2s 0
0)` and `(rot-z (* 360 t))` are errors. A param whose default is such a
quantity is given with its unit (`--param m=50g`, an import's value, a
tolerance's ends), whatever unit the default was written in. The product says
the unit of a pure number only: a formula that has already mixed a bare
number with a united one (`(/ rate 60rpm)`, rate bare) has a dimension
of its own, and is refused where another quantity is read, so say each
number's unit where it is written (`(/ (+ 0deg/s rate) 60rpm)`).

`t` is a time, in seconds, and a motion's formula carries it as any
quantity: a rate times it is a place, `(* 360deg/s t)`, `(* 20mm/s t)`;
a time is taken from it with its unit, `(- t 0.3s)`, and a share of one is
a pure number, `(/ t 0.5s)`; a sine takes an angle, `(sin (* 180deg/s
t))`, or a pure number (radians). A turn's drive is an angle and a
slide's a length: a bare number is one (a pose, in degrees or mm), a
number times t is a time and is refused, saying what to write. The rate
of a drive over t, `(diff angle t)`, is a turning rate or a speed (of a
pose held, a bare 0: `(+ 0deg/s (diff angle t))` reads either as a
turning rate). A formula that works in bare numbers takes t out as its
number of seconds, `(/ t 1s)`.

The unit a form's answers are in, and `(* 1UNIT x)` reads x in: masses
g, densities g/mm³, forces N, pressures and stresses MPa, torques N mm,
temperatures °C, times s, frequencies Hz, powers W, speeds mm/s, turning
rates deg/s, gravity m/s², conductivity W/(m K), specific heat J/(kg K),
convection W/(m² K).

A number that names another unit is converted: `(box 4in 2in
0.25in)`, `(param tq 12lbf*in)`, `(attr density 7.85g/cm3)`, `(* 3000rpm
t)`. The units:

```
length       mm cm m um in ft thou mil
area         mm2 cm2 m2 in2 ft2           volume   mm3 cm3 m3 mL L in3 ft3
mass         g mg kg lb oz                density  g/mm3 g/cm3 kg/m3 lb/in3 lb/ft3
force        N kN lbf                     torque   N*mm N*m lbf*in lbf*ft in*lbf ft*lbf
pressure     MPa Pa kPa GPa bar psi ksi Msi
temperature  degC degF K                  time     s ms min h
frequency    Hz kHz                       power    W mW kW hp BTU/h
speed        mm/s m/s in/s ft/s ft/min    turning  deg/s rpm rad/s
angle        deg rad                      turning acceleration  deg/s2 rad/s2
acceleration m/s2 mm/s2 ft/s2 in/s2       flow     mm3/s mL/s mL/min L/s L/min in3/s gpm cfm
viscosity    Pa*s mPa*s cP                conductivity  W/m/K BTU/h/ft/degF
specific heat J/kg/K BTU/lb/degF          convection    W/m2/K BTU/h/ft2/degF
heat flux    W/mm2 W/in2                  stiffness     N/mm lbf/in
inertia      g*mm2 kg*m2 lb*in2
```

- A power-of-ten conversion is exact (`2.7g/cm3` is exactly 0.0027);
  others are a rounded multiply. The inch, pound and pound-force are the
  international ones.
- `68degF` is the temperature 20 °C. A temperature difference has no
  literal: write it in °C. An accuracy or tolerance in degF or K, alone
  or in a sum, is an error.
- Only numbers carry units. Give a bare number its unit by multiplying:
  `(* n 1in)`.
- Where a form takes a known quantity (a package form's forces and
  temperatures, a claim's bound on a metric, `clearance`, `volume`, `res`,
  a joint's limits, a body's density), another quantity's unit is an
  error saying what it takes and where that number is written. A number
  keeps its unit through names, lets, calls, the branch an `if` takes and
  a `for`'s elements, and through sums, whose terms with a unit must be
  one quantity, every one of them checked (a bare term uses the
  language's unit).
- A product, a quotient, `square`, `sqrt` and a whole-number `pow` are
  the quantity their operands make: bare operands are pure numbers and
  only scale (`(* 2 p)`, `(* n 1in)`), and operands with units combine,
  so `(* 50lbf 15mm)` is a torque, `(/ 10N 2mm2)` a pressure and `(* 10kg
  9.81m/s2)` a force, 98.1 N (the number is put in the quantity's unit:
  ink's grams and m/s² alone would make it 98100). Units that cancel
  leave a pure number (`(/ 2in 4mm)` is 12.7, `(* 1mm (/ 1 1mm))` is 1).
  What they make that is no quantity above (a force times a time, four
  lengths) has no unit and is taken wherever a bare number is, its
  number in the base units g, mm, s, K and deg whatever its operands
  were said in and however they were grouped (`(* 2N 3s)` is 6000000
  g mm/s), and its powers are kept: it is never a pure number to be given a unit again,
  and is a quantity again when what is left is one (`(/ (* 2N 3s) 1s)` is
  2 N). The arithmetic is one, however it is reached: written, applied,
  folded, reduced, mapped.
- A sum, a difference, `min`, `max`, `neg`, `abs`, `floor`, `select`, a
  name (`tag`) and a promise (`lipschitz`) keep the dimension of what
  they are given. A sine, cosine or tangent of degrees, and `atan2` of
  two values of one dimension, are pure numbers. Any other function of a
  value with a unit (an exponential, a derivative, the sine of a length)
  is no quantity, marked so, with no powers: not a pure number, and a
  sum, a `min` or a `max` with such a term is no quantity either: it
  never takes another term's unit. (A bare term in a sum is taken in the
  other's unit, as it always was: every shape subtracts a length from a
  bare x.)
- An angle is a quantity, in degrees (`90deg`, `1.5708rad`); a form that
  takes an angle reads a bare number as it always did (degrees; radians
  for `sin`, `cos`, `tan`), and one said with its unit as the angle it
  is, so `(sin 30deg)` is 0.5. A turning rate times a time is an angle,
  degrees per mm times mm/s a turning rate, `(/ 360deg 60rpm)` a second.
  Arithmetic never drops or converts an angle: its power is carried
  like any other, so a product is one value however it is grouped,
  ordered or reached. `(* 1N*m 60rpm)` is "power with one power of
  angle", no quantity above, and stays that.
- A radian appears only where a value is read as a quantity: in a form
  that takes one (a speed limit, a mass), in a `fact`, in a printed
  unit. There a value that is that quantity times a power of angle is
  read with the angle in radians, as in every formula that has one:
  torque times a turning rate as power (`(* 1N*m 60rpm)` reads 6.28 W), a
  radius times one as a speed, a radius times an angle as an arc's
  length, power over a turning rate as a torque. To take the angle out
  inside an expression (a sum needs one dimension), divide by `1rad`:
  `(+ (/ (* r a) 1rad) 5mm)`. A bare number is a pure number, so
  `(/ 360 60rpm)` reads 57.3 s, 360 radians' worth: write `360deg`.
- A frequency and a turning rate are never read as each other: a cycle
  is not a radian, and a speed over a length is a shear rate as often as
  a spin. Each is refused where the other is taken, saying what to
  write, and the two ways a per-second arises differ by 2π: a speed over
  a radius is radians a second, `(* (/ v r) 1rad)` (a drum's or a rolling
  wheel's rate); a count of turns a second is `(* 360deg f)`; and a rate
  as turns a second is `(/ w 360deg)`. A
  root that halves an angle's power (`(sqrt (/ k I))` of a spring's
  stiffness per radian over an inertia) is no quantity, with no powers:
  it is refused against a frequency or a rate. Say the radians:
  `(* 1rad (sqrt (/ (* k 1rad) I)))` is the turning rate, and that over
  `360deg` its cycles a second.
- A joint's motion over t says its units, t being seconds: a rate where
  the number is written (`(* 60rpm t)`), an amplitude and the sine's
  turning (`(* 45deg (sin (* 180deg/s t)))`), a speed and a radius theirs
  (`(* (/ v r) 1rad)` is the turning rate). A number times t is a time,
  and a turn's or a slide's drive refuses it. A pose held is a bare
  number, in degrees or mm.
- A temperature in K, degC or degF is a point on its scale: multiplying,
  dividing or raising it is an error, in every scale alike. A temperature
  less a temperature is a difference, an amount in kelvins whatever the
  scales (`(- 100degF 50degF)` is 27.8): it scales and multiplies like
  any quantity, `(* 1W/m/K (- 300K 290K))`. A difference added to a
  temperature is a temperature. A form that takes a temperature cannot
  tell a difference from a point.
- When a parameter default contains a number with a unit or its value
  carries a unit, give units on command-line overrides, import arguments,
  selected cases, and tolerance deviations. For example,
  for either `(* 200 1GPa)` or `(+ 200GPa 0)`, write
  `--param modulus=210GPa` to override the default with 210 GPa.
  Bare values are accepted when the default's known unit and all its
  written units are native, as in `(param wall 2mm)`. Bare zero is also
  accepted on scales with a common zero. The zero exception applies to a
  scalar or a zero-width range at zero. If either command-line endpoint
  has a unit, write one on both. A bare default accepts an override in the quantity
  expected by the consuming form.
- An explicit pressure can exceed 1 MPa in magnitude. Bare pressure
  must be strictly between -1 and 1 MPa. A bare nonzero term in a sum,
  difference, minimum or maximum makes the result count as bare. A bare
  constant zero preserves the other term's unit in these operations.
- A bare value uses the form's documented native unit. Ink cannot detect a
  mistaken conversion, or a wrong unit of the same quantity, if the result
  is plausible. For example, bare `200` in an `E` clause means 200 MPa.
- Body, study and mass-objective densities are checked in g/mm³. A uniform
  density must be zero or in 5e-8..0.023. A graded or ranged density can
  fade to zero. It must stay in 0..0.023 and either be zero throughout
  or reach 5e-8 somewhere in the material and parameter range. A fluid
  density must stay in 5e-8..0.023. A published mass's implied density is
  checked when measured. A claim needs its measurement bounds inside the
  allowed range. `ink mass --density NAME` checks the named
  attribute on each body it weighs. A fluid's viscosity is checked in
  1e-7..1e9 Pa·s, thermal conductivity in 1e-5..1e5 W/(m K), and
  specific heat capacity in 0.001..20000 J/(kg K). During assertion,
  a material check on a placed field follows it to its body (a field a
  claim itself moves, rotates or mirrors is not followed). An invalid constant or a range on a setting fixed at load
  time can prevent the file from loading.
- Command options that are lengths or times take units (`--res 0.01in`,
  `--param wall=0.25in`, `--density 0.284lb/in3`). A `--param` value or
  range in a unit (a range's ends one quantity, as a sum's terms) is the
  param's, wherever the param goes: checked against its default's quantity
  and every form that uses it, at every value a run tries. `--points` and
  `--rays` are mm.

Answers are in the same units and don't name them: areas are mm²,
volumes mm³, inertia g mm², angles deg, and a joint's values deg (a turn)
or mm (a slide), its speeds deg/s or mm/s. Some keys name their unit
(`mesh_ms`). A claim's `value` and `accuracy` are in its quantity's unit
(a clearance's mm, a volume claim's mm³, a metric's), else the model's
own numbers.

## Angles

Angles are degrees wherever a form takes or reports one: placement steps
(`rot-x`, `rot-y`, `rot-z`, `rot`), std's `rotate-x/y/z` and `twist-z`
(degrees per mm), gear pressure angles, a turn joint's drive and limits.
Turns are right-handed: 90 about x takes +y to +z, about y +z to +x, about
z +x to +y. A turn by a constant multiple of 90 is exact.

Radians are the exception:

- `sin`, `cos` and `tan` take radians; `atan` and `atan2` return them.
  Of an angle said with its unit, `(rad->deg a)` is its number of
  degrees and `(deg->rad a)` its number of radians, pure numbers:
  `(rad->deg 1rad)` is 57.29578.
  `(deg->rad d)` converts, with pi exact where d reads x, y, z or t
  (constants and params are their floats, as everywhere but a value
  claim on constants): an angle
  wrapped at a whole turn (`(mod d 360)`) reads the same either side of
  the wrap, so motion across it is proven continuous (`(* d (/ pi 180))`
  rounds pi: not). A wrap counts as whole turns only where it is 360 k
  degrees exactly (1e20 is not).
- `path-arc`'s sweep, and the angles std's belt and chain tools return
  (`sync-belt-theta`, `roller-chain-theta`, `belt-route-entry-angle`,
  `belt-route-sweep`).

`sin`, `cos` and `tan` take a pure number, in radians, or an angle with
its unit, which they convert: `(sin 30deg)`, `(sin (* 360deg/s t))`. A
time is neither: `(sin (* 2 t))` is an error that says to write the rate.

## Fields

A field is a real-valued expression over the inputs `x y z t` and any
declared parameters. A **shape** is a regular solid: the interior of the
closure of where its field is `< 0`. So a 0 with the shape on every side
(faces placed flush, inside the material) is inside, a 0 with outside on
one side is the surface, and a region where the field is 0 that has volume
is outside: the closure does not reach across it. Two blocks written flush
are one solid, a cut that ends exactly on a face leaves no skin, and
`subtract` keeps a point where the cut is exactly 0 on a sheet with
material round it (a cut of no thickness removes nothing). `subtract` is `max(a, -cut)` on
the fields: the regular difference wherever the cut's field is not 0 over
a volume. Where it is, the result's field is 0 over that volume, which
by the definition is outside, and a claim reads it so: `ink check` warns
where a body has one, naming a point in it. A cut meant to be empty
somewhere is written positive there, not as `(min 0 g)`.
A slit is a gap with a width: `(subtract part (box 0.2 50 50))`, not two
blocks written flush. The shape need not be an exact distance, though distance-like fields
(|gradient| near 1) mesh best. The `clearance` claim is sound for any
field (Claims).

`(offset-shape f d)` grows f's shape by d: everything within d of it (for
d < 0 it shrinks it: what lies farther than -d inside). It is f - d,
exact where its construction shows f is that distance: x, y and z moved
by constants, slabs, spheres, cylinders and boxes along the axes (std's,
or a max of slabs, built as `emax`, the same shape), a union grown, an
intersection or subtract shrunk. Elsewhere `ink check` warns: with f's
slope at most L (the kernel's bound, promises included), a growth holds
everything within d / L and may reach farther (mitred corners, a 0 the
shape keeps out filled). Of samples of a distance (`import-stl`,
`redistance`) it is approximate, to within their error plus sqrt 3 / 2 of their
spacing, and `ink check` says so. `(- f d)` is the field's own growth.

The kernel's operations are the whole vocabulary; everything else is
library code:

```
(+ a b ...)  (* a b ...)  (min a b ...)  (max a b ...)   variadic
(- a b)  (- a)  (/ a b)  (atan2 y x)  (pow a b)
(neg a) (abs a) (sqrt a) (square a) (exp a) (ln a)
(sin a) (cos a) (tan a) (floor a)
(deg->rad a) (rad->deg a)  degrees in radians and back, pi exact
(sin-cos-deg a)         [sine, cosine] of a turn in degrees
(select c a b)          a where c < 0, else b
(== "a" "b")            1 for byte-identical strings, else 0 (case-sensitive, no normalization); broadcasts over vecs; a string and a number are a load error
(real a)                use the value of a as an operand of real arithmetic
(emax a b ...)          Euclidean max: sqrt(a^2 + b^2) where both are > 0, else max
(warp e fx fy fz [ft])  e with x:=fx, y:=fy, z:=fz (and t:=ft)
(lipschitz k f)         f, promised |grad f| <= k in x y z
pi inf x y z t          constants (inf: +infinity) and inputs
```

On shapes, `min` is the union and `max` the intersection. std's `(subtract a b)` is `(max a (neg
b))`. `emax` is the exact intersection of exact fields with orthogonal
gradients: `(emax qx qy qz)` is an exact box, `(emax g dz)` an exact
extrusion.

### Conditional estimates

`(estimate-when value condition "reason")` returns the numeric field
`value` with its units. A finite, nonzero condition (positive or negative)
admits an estimate. Zero or a nonfinite condition makes a claim reading
the result unknown with the given reason. Over a parameter range the
condition must be shown to hold throughout the range. The condition may
read parameters and time and must be independent of position. A time
`warp` substitutes for explicit `t` in the condition and value. A
derivative retains the condition and estimate status. Over a motion the
condition is shown over the span, or over each step between the moments
read (Motion, below): a step where it is not shown is not settled, naming
the step and the reason. What reads the value at other times than the
model's moment (a transient study on its own clock, a model that steps
the whole motion) has the condition shown at every one of those times,
or it has no answer, whatever the moment asked. Beside tolerances, the
condition is shown over their ranges too, at each moment and over each
span of t. One limit: a `warp` of a solved parameter can check the
condition at the wrong moment.

### The one promise

`(lipschitz k f)` promises f's slope in x, y, z is at most k, for every t
and param value. Proofs (clearance, distance) trust it instead of bounding
the slope from the formula's parts (interval arithmetic), and so does f's
bound over a box: within k times the box's reach of f at its centre. A false promise makes proofs unsound;
std makes it only for its exact primitives. A k that is not positive
promises nothing.

- A `warp` keeps a promise, scaled by how much it stretches space; rigid
  moves keep it as is. Any other rewrite under a promise drops it.
- A promise in a function holds for fixed inputs. Where the promised field
  reads an input that depends on x, y or z (`(sphere (+ 10 (* 9 (sin
  x))))`), it is dropped. Put a coordinate map in a warp, `(warp
  (lipschitz 1 (g x y)) s n z)`, not in the promised field.
- A promise holds for its constants as computed: a proof does not use it
  where one is infinite or comes of an undefined case (`inf / inf`).
- A k that reads x, y or z, under a warp whose coordinates jump (a
  `select`, a `floor`), is not used: over a box the frame's coordinates
  between the centre's and a point's may be ones the warp never reaches,
  where k is not read.

### Value rules

Every op is total and the same on every implementation: its IEEE-754
double result, except that there is one zero, `x / 0` is 0, and anything
else undefined is 0 (`inf - inf`, `ln` of a negative). `sqrt` of a negative
is 0 too, but defined: `sqrt(max(a, 0))`.
Infinities are values: overflow gives ±inf, `ln 0` is -inf. JSON has
none, so an answer gives an infinite number as `null`. The
transcendentals are Rust's `libm` everywhere. conformance/README.md
states the rules exactly.

`ink check` warns where its sampled survey of a body's field finds it
passing through an undefined case, naming the construction (a 0 in a
shape is not inside it). `ink eval` lists them per point under
`undefined`; one in a number used while the file loads is reported at its
line. Only a branch taken counts. Proofs cover every value an undefined
case could take.

An expression of constants and params only is its f64 value; one reading
a ranged param, x, y, z or t is proven as a real number. Arithmetic on
x, y, z or t that cancels as written, `(- x x)`, `(* 0 x)` or
`(* 0 (cos (* 90deg/s t)))`, is an exact 0: a constant to everything that
reads it, folded with other constants only where the float is the real
value, else into a band about the float that covers it, so
`(- (+ 1e16 (+ (- x x) 0.5)) 1e16)` is 0.5 within its rounding, never 0.
`ink eval` says inside by the field's outward-rounded bounds at the
point, as a probe claim does: where they reach across 0 it gives them as
`band` and `inside` null.

Use `(real top)` before computing a phase length such as
`(/ (real top) acceleration)` to keep subsequent products and quotients
in real arithmetic. Wrapping an already rounded result keeps that value.
At an affine time switch the rate checker can compare the two one-sided
values exactly, even when the switch lies between representable times.
Any nonzero gap still precludes a finite speed.

Std's `(trapezoid-position elapsed acceleration top-speed distance)`
gives a rest-to-rest move with positive acceleration and speed, and
distance at least `top-speed² / acceleration`. It holds zero before
the start and the requested distance after deceleration.

### Derivatives

`grad`, `diff` and `ink eval`'s gradient take each op's
derivative on the branch its value took, by the chain rule, and agree bit
for bit. At an exact tie (`min` or `max` of equal values, `abs` at 0) a
derivative is the average of both sides. `select`'s is its branch's,
`floor`'s 0. Proofs bound a derivative over a region by its
formula's range (Resolution: where a slope switches branch).

## Values

- **field**: a number or number expression. A literal outside a double's
  range (`1e999`) is an error.
- **string**: `"text"`, with escapes `\n`, `\t`, `\"` and `\\`; any other is
  an error (write a path with `/`).
- **vec**: `(vec a b ...)` of any values; `(nth v i)` picks one.
- **data**: `(data ...)`, opaque structured data kept verbatim as JSON
  (symbols become strings, lists arrays).
- **frame**: a rigid placement, an origin and unit x, y, z axes. A point,
  an axis (its z line) and a plane (its xy plane) are all frames.
  - `(at step ...)` builds one. Steps: `(move dx dy dz)`, `(rot-x deg)`,
    `(rot-y deg)`, `(rot-z deg)`, `(rot ax ay az deg)`, or any frame.
    `(at)` is the identity. Steps may read params, t and local names, never
    x, y or z. Each moves what is placed about the outer axes and origin,
    first step first. `(at (rot-x -90) (move 0 6 20))` turns a part, then
    puts it at (0, 6, 20). `(at (move 0 6 20) (rot-x -90))` puts it there,
    then turns it, offset and all, about the outer x axis, to (0, 20, -6).
    So turn first, then move, unless the turn should carry the offset (a
    bolt circle: `(at (move 20 0 0) (rot-z 60))`). Likewise `(at (move 0 0
    5) F)` is 5 mm along F's own z; `(at F (move 0 0 5))` is F moved along
    the outer z.
  - `(axes (origin ox oy oz) (x xx xy xz) (z zx zy zz))` writes one by its
    origin and directions, in any order, as `ink check` prints it. z keeps
    its direction (a face's normal). x is made square to it (x along z is
    an error). y is z cross x. Neither needs unit length. A clause takes 3
    numbers or a vec of them, in params, t and local names (A part on a
    tilted face).
  - Convention: a frame where parts touch sits on the touching face, z out
    of the part. x is the part's reference direction (a bolt pattern, a
    rail), as its family states.
  - `(frame-of PATH [NAME])`: frame NAME of the body or assembly PATH ends
    at, placed (Assemblies). Without NAME, PATH's placement.
  - `(attr-of PATH NAME)`: attribute NAME of the body PATH ends at, as
    written, read at load time (a catalog part's record: `(attr-of
    root/screw spec)`).
  - `(met-by PATH)`: the name of the child frame placing instance PATH
    (its CFRAME, Assemblies). An error for one placed by none.
  - `(get F "origin")` is a vec of 3 fields; so are `"x"`, `"y"`, `"z"`.
  - `(select c F G)`: F where c < 0, else G, c of params and t. A part
    that changes hands: `(instance cup cup (select (- t 1) (at (move 40 0
    0)) (frame-of gripper grasp)))` rests at (40, 0, 0), then rides the
    grasp frame from t = 1.
  - `(place SHAPE F)`: SHAPE, written in F's coordinates, put at F (a
    keep-out at a shaft tip, a hole at a frame).

## Load-time programming

Loading a file (load time) builds a finite math graph. The language that
builds it is small and functional:

```
(lambda (a b) body)            function value; captures local names; (lambda (a rest ...) body)
                               takes the arguments after a as the vec rest
(map f list ...)               element-wise (equal lengths)
(reduce f list)                balanced tree; f must be associative (min, max, +, * or any fn)
(fold f init list)             left fold, sequential
(apply f list)                 call f with the list as arguments
(range n) (range a b [step])   numbers known at load time (std's (linspace a b n): n of them from a to b)
(len l) (concat l ...) (slice l a [b]) (flatten l) (nth l i)    negative indices count from the end;
                               concat of strings is a string; counts and indices are whole numbers,
                               never rounded; a range has at most 1e6 elements
(char-codes "text")            a string's characters as their Unicode code points
(if c then else)               c must be known at load time; use select for spatial choice
(< a b) (<= a b) (> a b) (>= a b) (== a b)                    1 or 0 (as fields)
(read-json "file.json")        numbers, strings, arrays; objects become data (paths: relative to the file
                               the expression is written in)
(get data "key" [default])     default (evaluated only then) where data has no such key
(bounds f [default])           proven bounding box, as (vec (vec x0 y0 z0) (vec x1 y1 z1)); default
                               (evaluated only then) where f is empty, else that is an error; of a
                               field that reads t, over its file's (t-span) stated before it and t = 0,
                               and of a ranged param, over its band and value: runs beyond are refused
(grid-bake f [(res mm)])       sample a field onto a grid (trilinear; a cache for expensive fields)
(import-stl "part.stl" [(res mm)] [(region (x0 y0 z0) (x1 y1 z1))])
                               a triangle mesh as a signed distance field (sampled; with a
                               region: only the part of the mesh inside it, sampled finer)
(offset-shape f d)             f's shape grown by d, everything within d of it (Fields)
(redistance f [(res mm)])      signed distance to f's surface, sampled on a grid (approximate)
(diff f v)                     derivative with respect to x, y, z, t or a param, or a def that is one:
                               symbolic, a field like any other; std's (grad f) is (vec dfdx dfdy dfdz)
(unit v)                       v / |v|, the unit vector along a vec of fields
(check TEST value)             value, with a condition on what built it, which ink check reports where it fails
(error "text" [v ...])         a load error saying text, then each v, where it is evaluated
(profile-reach f)              how far a 2D profile reaches from its origin, proven (to 1/1000 of its size); over t and bands as bounds
(ranged lo hi [at])            a variable every claim reading it ranges over [lo, hi] (Tolerances)
(polygon pts [drafts z0 z1])   a closed profile of lines, arcs and rounded corners (Standard library, Curves)
(revolve-z g)                  a 2D profile turned about z (Standard library, Turned profiles: revolve-z)
```

- Function names are values: `(map square l)`, `(reduce min ds)`. So are
  the built-in forms that are functions (`vec`, `map`, `range`, ...). The
  syntax forms (`let`, `if`, `lambda`, `tag`, `read-json`, ...) are not:
  wrap one in a `lambda` to pass it.
- Kernel operations broadcast over vecs element-wise: `(- here a)`, `(* v
  2)`, `(dot a b)`. Vec lengths must agree.
- A `fn` or `lambda` cannot repeat an argument name. A sequential `let` may
  rebind one.
- There are no runtime loops and no recursion: functions expand at load
  time, so every model can be bounded and differentiated. For generative
  designs, write a program that emits IR.
- Load-time computations use the current `--param` values; the fields stay
  differentiable in their params. What a ranged param cannot be used for
  is under Tolerances.
- A direction must not be zero: where it reads no x, y or z, `(unit v)`,
  a `rot` axis, a `(normal ...)`, `(axis ...)` or `(slide ...)` refuse a
  zero length. One over space is a field, its zero reported at samples
  like a division by zero.

**Conditions on what a function builds.** A valid field can still be a
wrong thing (a hose through itself, a hole that leaves no wall). `(check
TEST value)` is `value`, with a condition on the inputs it was built
from. `ink check` warns where one fails, and a
model passes only with no warnings. TEST takes numbers (constants, params,
t):

- `(value E CMP V)`: E compared with V (`<`, `<=`, `>`, `>=`);
- `(path-thickness PTS >= R)`: the path PTS (a list of `(vec x y z)` or
  `(vec x y)`) is at least R thick, R being how far what is swept along it
  reaches from it.

A condition is evaluated once, at the params and moment `ink check` runs
at (`--param`, `--t`), not over ranged params or a motion: check a
moving part at the corners of its travel. It is reported in the file that
states it and in every importer.

```lisp
; A tab with a hole: the hole must leave 1.5 mm of wall each side.
(fn tab (w t hole)
  (check (value (- w hole) >= 3)
    (subtract (box w w t) (cylinder-z (/ hole 2) (* 2 t)))))
```

The warning names the form, the call and the params' values.

**Paths.** A body swept along a path overlaps itself where the path is
thinner than the section reaches. std's `tube`, `pipe`, `sweep` and
`sweep-varying`, and the catalog's path builders, state this rule with
their section's reach, so what they build is checked. Build your own
path-following body on them, or state the rule in it, with R the farthest the
body reaches from the path (`profile-reach` gives a 2D profile's). A
path's thickness is the least of half the distance between two of its
stretches that face each other and, at each turn, what keeps the mitres
at its segment's ends apart (a U-turn needs segments of at least 2R). The
check allows half the model tolerance on the section's reach (0.025 mm at
the default tolerance). A path whose first and last points coincide is a
loop. The warnings say where it fails and how to fix it. For a path in a plane and a body
extruded square to it, R is the reach in the plane:

```lisp
; A flat strap along a path in the xy plane: its section is 1.2 mm thick.
(fn strap (pts width)
  (check (path-thickness pts >= 0.6)
    (extrude-z (- (polyline-distance-2d pts) 0.6) (* -0.5 width) (* 0.5 width))))
```

**Make fields boundable.** Interval arithmetic must be able to prove where
a shape is not, out to infinity (Where analysis looks). Products of
signs, selects between far-reaching pieces, infinite extrusions,
and a mix `a + (b - a) k` (use std's `blend`) can defeat it. Intersect
such fields with a containing shape, and clip each piece a `select`
chooses to its own region. `ink check` warns about bodies not proven
bounded.

### Sampled fields

`grid-bake`, `redistance` and `import-stl` make sampled fields: values on
a regular grid, read by trilinear interpolation, with exact interval
bounds and derivatives of the interpolant. They behave like any field.
Outside their grid they continue as the
boundary value plus the distance to the grid. Sampled data enters a model
only this way.

- `redistance` is approximate, with no error bound. It needs a bounded
  field whose surface its grid resolves (default res: the largest extent
  / 128).
- `import-stl` imports the solid the mesh encloses, by even-odd crossings:
  overlapping shells cancel, so import overlapping parts separately and
  take their union. The mesh should be closed; an open one is sampled as
  its facets face, and `ink check` warns. Binary and ASCII STL are read,
  ASCII strictly (anything else is an error at its line). `(region ...)`
  samples only that part of the box, finer, and cuts the solid to the
  region: use it where a fit matters, beside a coarse import of the whole
  part.
- A `(res h)` here (and a package form's) is a cell size known at load
  time. A grid's spacing is at least 1e-4 of its box's largest extent (a
  finer res is coarsened, with a warning), so sample a small fit over a
  small `region`. A grid holds at most 64 million values, as f32.

A claim that reads a sampled field is an estimate (Results: `samples`);
one that reads a result whose form measures a band it did not measure for
it is unknown, saying so.
Sampled data cannot shape the model at load time (a `for` list, an `if`,
a count, a form's setting): that is an error. `(bounds f)` of it bounds the samples, its
numbers measured like them: size and place shapes with them. Found empty,
it is an error, default or not.

## Repetition: `for`

`(for var list form ...)` repeats forms, once per element of `list`
(evaluated at load time), with `var` replaced by it and `$var` in names by
its text. A nested loop's list can read the outer variable. In its body,
the inner binding shadows an outer one of the same name. `$i2` names
`i2`; punctuation can follow a placeholder (`bolt-$i-head`). Text inserted
for a placeholder is used literally. It nests, and works:

- at top level;
- in assembly, body, relation and package form clauses;
- among an instance's clauses after its placement (WHERE), its list
  reading the assembly's arguments;
- among an `except`'s paths (Paths).

A top-level `for` expands where it is read: its list reads the
names defined before it.

```lisp fragment
(def n 6)
(for i (range n) (param bolt-d-$i 5))                         ; params, defs, bodies, asserts
(assembly root (for i (range 6) (instance bolt-$i bolt (at (move 20 0 0) (rot-z (* 60 i))))))
```

## Top-level forms

```
(import "file.ink")                    load another file, by a path relative to this one or absolute;
                                       else from $INK_PATH's directories (Libraries, below)
(import "file.ink" as p)               ... under a namespace: its names read p.name here (p.bracket,
                                       p.root/arm, --param p.size=4), so files never clash. std stays
                                       global; files p imports without `as` join p's namespace
(import "file.ink" (name v) ...)       ... setting its params' defaults for the whole load, so one value
                                       per name, known at load (a param here passes its value); (name v
                                       (tol lo hi)) sets the band too. --param overrides
(import "std")                         load the standard library (see: ink std)
(members spec [its-name])              the file draws a family, a part for each member of its table in
                                       families.json; spec is the member's record (Families)
(param name default)                   design parameter; override with --param name=v
(param name default (tol lo hi))       ... with a tolerance, proven over (Tolerances)
(param name default (cases v ...))     ... and other values ink require asks every claim at (Cases)
(param name "text" [(chars S)] [(len lo hi)] [(cases "..." ...)])
                                       a text param (Text params)
(def name expr)                        named value; names are immutable. A def whose value is a
                                       lambda is called by name, like a fn
(fn name (a b ...) body)               function, expanded when the file loads
(fn name (a rest ...) body)            ... whose last argument collects the rest, as a vec
(let ((a e1) (b e2)) body)             local names, bound in order; usable anywhere
(t-span A B)                           the motion: require checks this file's claims that read t at
                                       every moment from A to B; interfere checks pairs that read t over
                                       the least span holding every file's (Motion)
(body name (shape f) (attr k v) ...)   a part: one shape and any attributes (Bodies)
(body (name (p default) ...) ...)      ... with arguments, which instances pass
(assembly name clauses...)             a tree of placed instances (Assemblies)
(assembly (name (p default) ...) ...)  ... with arguments
(require name claim [options])          a checkable design claim; options (res r), (accuracy x) (Claims)
(require name (within-ratings PATH [(excluding P ...)]))
                                       asks whether every rating of every part at or under PATH holds,
                                       less the parts P left out on purpose (Facts and requirements)
```

Packages add more forms (Data for other packages).

**Sizing a part.** A part is sized by its arguments: a file's size
params are its body's and root's arguments, defaulting to the params,
`(body (rail (length length)) ...)`. An instance passes another size,
`(instance r (p.rail 470))`, so a part at two sizes is two instances of
one import. An import's value sizes the whole file, its claims and
studies too: `(import "rail.ink" (length 470))`.

Write dimensions from the params, `(probe 0 0 (+ base-t 1))` not `(probe
0 0 6)`, so claims follow a new value. x y z t pi inf cannot name a
param, def or fn. A file's own may take an operation's or form's name
(`tag`, `len`): what the file writes after it reads its own, and std and
its imports the built-in. Reading the built-in where it is written
before it (a function's body too), or an import's own one unnamespaced,
is an error.

**Reading order.** Each form sees only the names defined before it: in the
file, in files imported before it, and in std. A name read before its
definition is an error saying where it is defined.

**Names are lexical.** A `fn` or `lambda` reads the names of its own file,
its imports and std. Local instance paths in assembly clauses belong to
that assembly. Functions and default arguments read the frames in their
defining scope. When a frame path is read outside the scope of the
assembly being built and its first name denotes both an instance of that
assembly and a body or assembly in the defining file, ink reports an
ambiguous frame. To use the instance's placement, read its frame in the
assembly and pass it as an argument. To use the file's frame, bind it in
a top-level def or rename the conflicting instance. So a part's `(box ...)` is std's
box whatever its importer calls `box`, and a path a function reads
(`read-json`, `import-stl`) is relative to its own file. A file may define
a name std defines, once: it shadows std's for that file and its
importers, while std's functions keep reading std's. Any other name is
defined once. An import cycle is an error.

Each namespace loads its own copy of an import, including its parameters,
cases, and values read from files beside the import path. Repeated imports
in one namespace refer to one copy. If two paths to the same source have
different directories for relative reads, give the second copy its own namespace.

**Libraries.** An import not found beside its file is read from the
first directory of `$INK_PATH` (separated as `$PATH`'s are) holding it,
and named by its path there in errors, claims' sources and `ink check`'s
imports. A parts catalog is `parts/` in such a directory.

### Families

One file can draw many parts: a family, each part a member of its
table. `(members spec)` in `F.ink` makes it one: the `families.json`
beside it holds, under `F`, its table, each member by its name and
record, and `spec` is the record of the member the file is drawn as
(data, as `read-json` gives it). A member is named by a path of its
own: `DIR/NAME.ink`, where `DIR/families.json` names member NAME, is
that member, the family file loaded with `spec` its record, and
`DIR/NAME.json` is its record to `read-json`. So a member is imported,
checked and found like any part (`(import
"parts/fastener/ISO4762-M3x8.ink" as s)`, `ink check
parts/fastener/ISO4762-M3x8.ink`), its path its identity and its name in
errors and imports, though no file holds it. A name is one part: a
file of a member's name (`NAME.ink` or `NAME.json`) beside it is an
error. The family file loaded as itself is its table's first member.
`(members spec its-name)` also names the member's own name, a string
(`"ISO4762-M3x8"`): what a file of its own would have written out.

```lisp fragment
; cube.ink
(import "std")
(members spec)
(def size (get spec "size"))
(body part (shape (box size size size)) (attr spec spec))
```

<!-- not run: the families.json beside cube.ink -->
```json
{"cube": {"record": {"family": "cube", "size": 1},
          "members": {"cube-S": {"name": "small cube"}},
          "groups": [{"record": {"size": 3}, "members": {"cube-L": {"name": "large cube"}}}]}}
```

`cube-S.ink` is a 1 mm cube, `cube-L.ink` a 3 mm one. families.json
holds an entry for each family file beside it, by its name less
`.ink`. A table is a group: `record` (what each of its members' records
has), `members` (each name's own) and `groups` (groups within it), each
optional. A member's record is its groups' records, outermost first,
then its own, merged: objects key by key, the outer's keys first and in
their order, an inner value over an outer one; any other value
replaced. Its members come in order: a group's own, then its groups'. A
member's name, and a family's, is a file name less `.ink` in that
directory; a member's is once there and no family's. What a member is read from (its hash, `ink parts --export`) is
its family file and its own record, not its table's others.

### Text params

`(param name "Ada" (chars "ABC...xyz -") (len 1 12))`
is a text read while the file loads: a name, a label, what a customer
types. `--param name=Grace` sets it (a value a number param takes is text
here too: `--param code=042`); a text holding a comma, or starting with
`"`, is written as a JSON string, `--param 'name="Lovelace, Ada"'`
(over MCP, `{"param": ["name=\"Lovelace, Ada\""]}`). An import sets one as
it does a number, `(import "sign.ink" (name "Grace"))`. `(chars S)` is
the characters it may hold, `(len lo hi)` how many; a value outside them,
the default included, is refused at load, saying which character or how
long. `ink check` lists each with its `value`, `chars` and `len`, which a
client asking a customer checks first. It is a string everywhere it is
read (`char-codes`, std's lettering), never a field, so no claim ranges
over it: each answer is at its value, and `ink require` says which
(`texts`). A claim about every text it allows reads a bound over all of
them (std's `lettering-max`: the widest any n of the allowed characters
can be set), or is asked at worst cases the file names (Cases).

### Cases

`(cases v ...)` after a param's default (a number param's numbers, a
text param's texts) names other values `ink require` asks every claim at:
a count or a variant read while the file loads, the widest name a sign
takes. A claim holds where it holds at the param's value and at every
case, each param at its value or a case (every combination, at most
256), fails where it fails at one, else is unknown. Its row is its
answer at the values, its `status` and `accuracy` the worst of all, and
`cases` each case's answer, `case` the values there; each case is its
expression's value there, every param set (one may read another); the
answer lists the combinations asked (`cases`), and those not, where a
text param does not allow the text it takes there, its default too
(`cases_not_allowed`), and those not loaded before the call's time
limit (`cases_not_loaded`: every claim unknown there). Cases go on
params every combination declares, as many at each (their values may
differ): on one a `for` over a param makes at some alone, they are
refused (it takes its value there), as are cases a combination declares
more or fewer of; a claim made at some alone is reported of those
(`at_values` false). A param
`--param` sets is asked at that value alone. Cases are values checked, not a proof over the values between
them: a range a proof covers is a tolerance or a `(ranged lo hi)`.

```lisp
(import "std")
(param n 1 (cases 2 3))            ; how many pegs
(def pegs (apply min (map (lambda (i) (move (cylinder-z 2 10) (* i 10) 0 0)) (range n))))
(require fits (within pegs (-3 -3 -6) (27 3 6)))   ; holds at 1, 2 and 3
```

### Bodies

```lisp
(import "std")
(body bracket
  (shape (box 40 20 5))
  (attr color (vec 0.1 0.1 0.12))  ; sRGB numbers or fields in 0..1, or "#RRGGBB"
  (attr metallic 0)              ; 0 (default) .. 1
  (attr roughness 0.4)           ; 0 polished .. 1 matte, default 0.4
  (attr transmission 0.9)        ; 0 opaque (default) .. 1 clear
  (attr density 1.24g/cm3)        ; any name; numbers or fields (density: g/mm^3)
  (attr material "pla")           ; a std/materials.ink name, or data
  (attr material-note "matte black, 40% infill")
  (attr tolerance 0.2mm)          ; as made: its claims' accuracy (Claims)
  (attr spec (data (tensile 50 "MPa")))
  (frame mount (at (move 0 0 5)))  ; a named frame, in the body's coordinates
  (mark scribe (vec (- x 20) (box 40 20 5))))  ; a curve drawn on it (Marks)
```

An attribute has any name and a number, string, data or field (a graded
material) as its value. One read as a number (`(attr NAME)` in a claim,
`by-body`, `ink mass`'s density) is the body's own `(attr NAME ...)`, else
its material's NAME: a body given only `(attr material "pla")` has pla's
density. An envelope (a body modelling only a part's outside) may give its
published `(attr mass 199g)` in place of a density, and an assembly the
published mass of the bodies it places (Assemblies). A body with none of
these has no mass: `ink mass` says so, and a claim reading its `(attr
density)` is an error naming it.

Material names resolve to the grades and conditions in `std/materials.ink`;
accepted legacy aliases are listed in `std/README.md`. `(material NAME)`
returns numeric properties plus `name` (as requested), `canonical-name`,
`condition`, `evidence` (a `std/material-sources.json` reference), and
`strength-kind` where a strength is supplied. Each property's adjacent
comment and manifest entry give its source, kind and test conditions.
Specification minimum yield applies only to the named temper, product form
and size range; physical values and maker coupon strengths are typical.
Brittle strengths retain their test name (`bending`, `tensile`,
`compressive`); none is interchangeable with ductile `yield`.

A strength (`yield`, `tensile`, `compressive`, `bending`) read by `get`
whose `strength-kind` is not a specification minimum is a maker's typical
coupon result: a claim that reads it is an `estimate`, naming it in
`typical_strengths`, so it never holds as a guarantee. Typical physical
values (E, density, conductivity) do not mark a claim so: inspect their
evidence before resting a limit on them.
Printed defaults apply only to their
stated process and direction; stock-shape data does not qualify prints.
Explicitly `UNSOURCED` physical reference
values remain unqualified inputs.

**Looks.** `color`, `metallic`, `roughness` and `transmission` only change
how a body is drawn. Each may be a field, clamped to 0..1, that follows
the instance's placement: `(attr metallic (select (- z 20) 0 1))` is a
plastic housing with a steel shaft in one body; std's `paint` gives each
piece of a body its own look. A number outside 0..1 (`(vec 200 30 30)`) is
an error. Without a colour, a body gets a palette colour. A colour
reading a package form's result shows how the result varies over the
body.

**Marks.** `(mark NAME CURVE ...)` draws curves on a part that are not
its shape: a sheet's bend lines, a scribe line. A CURVE is `(vec g r)`,
where field g is 0 inside field r < 0 (r: where the curve runs), or a vec
of them. A mark is where g changes sign within r < 0, as a part's
outline is where f does: a g that only touches 0, or a mark lying in the
slice's plane, draws nothing there. `ink slice` gives a mark's slice within
the part's bounds as curves, open or closed, written to a DXF dashed on
layer NAME, which no target's layer may share.

**Arguments.** A body takes arguments as an assembly does: `(body (spring
(length 30)) (shape (compression g length)) ...)`. Its shape, attributes
and frames read them like local names. An instance passes values of params
and t: `(instance s (spring (- 30 travel)) ...)`. So a part that deforms
with a mechanism (a spring compressed by its joint, a belt running with
its pulley) is one shape driven by its input. An argument that reads t
makes the shape read t (Motion). Instances are shared and named as an
assembly's are (Assemblies: Arguments). A body's clauses read only its
arguments, and their `frame-of` paths start from the top level. See
`examples/spring-press.ink`.

### Assemblies

```lisp fragment
(assembly name
  (instance NAME CHILD)                       ; CHILD: body or assembly, at the identity
  (instance NAME CHILD WHERE [CFRAME] [MOTION] [(limits LO HI)] [(speed-limits LO HI)] [CLAUSE ...])  ; CLAUSE: a package's
  (instance NAME (CHILD ARG ...) ...)         ; a body or an assembly given arguments
  (frame FNAME F)                             ; a frame, in this assembly's coordinates
  (contact A B ZONE)                          ; overlap of A and B inside ZONE is intended
  (relate A B KIND)                           ; what A and B are to each other: a zone, claims, data
  (attr mass 280g))                           ; the published mass of the bodies it places
```

For a published assembly mass, `ink check` reports `met: null` when
the known bodies' bounds still straddle the total after refinement, or an
input has no proven bounds. An enclosing assembly uses the proven density
bounds of each inner published mass. If an inner mass's feasibility is
undecided, the enclosing assembly's check also reports `met: null` and
names the inner mass it depends on.

**Placement.** WHERE is any frame in this assembly's coordinates: `(at
...)`, `(axes ...)`, `(frame-of ...)`, a def. Without CFRAME the child's
own axes line up with WHERE: `(instance n c (at (rot-z 90) (move 10 0
0)))` turns the part, then moves it 10 mm along this assembly's x. CFRAME
names a frame the child declares, met face to face with WHERE, exactly:
origin on origin, z opposite z, x along x. An assembly exports a frame by
declaring one, usually `(frame tip (frame-of l2 tip))`.

### Facts and requirements

An instance may state numbers about the part
it places, among its clauses: `(fact WHAT VALUE)`, or `(fact WHAT WHERE
VALUE)` on a frame the child declares, on a part under it, or on that
part's frame (`shaft`, `arm/motor`, `arm/motor/shaft`), e.g. `(instance
motor motor.drive (fact radial-load shaft door-weight) (fact speed shaft
22rpm))`. VALUE is any number of params and units, not of x, y, z or t.
(A file that defines a `fact` of its own keeps it: the clause is not read
there.) A part carries its ratings as data (`(attr ratings ...)`, written
by its generator). `(require NAME (within-ratings PATH))` at the top
level asks whether each rating of each part at or under PATH holds
against the facts stated where the part is placed. Its answer says which
parts under PATH carry ratings and, first, which of those the model
loads carry none (a motor, a screw and its nut, a named support): while
one does, the file is not all yes. A part whose maker states no limit
says so, `(attr no-rating "WHY")`, and is listed apart: nothing could
check it. The file is all yes only where every rated or loaded part of
the design is under some requirement or left out by name
(`parts_under_no_requirement`). `(within-ratings PATH (excluding P
...))` leaves the parts P under PATH out of the question on purpose,
and the answer says so. What a part works out, a formula's fact or a
duty's (`mean-speed`, `duty-load`, `rms-torque`), is not the design's to state: an
error, where it is written. A frame an assembly
writes as exactly `(frame NAME (frame-of INSTANCE ITS))` is that
instance's frame ITS (a fact on a motor's `shaft` is on its rotor's); a
frame and an instance of one name cannot both be meant, and a fact that
could be on either is an error.

A fact has the unit of what its value is, by the language's rules
(Units): a number written with its unit, a name of one, a sum of such,
and a product or quotient, which is the quantity it makes (`(* door-weight
pitch-radius)`, a force times a length, is a torque). A bare number, a
sum with a bare term and what is no quantity of the language's have
none, and are refused against a rating, with the line and why; so is
another quantity's. Never read in the rating's unit. A rating its data marks `sizes` reads its facts without
their sign. A fact that reads a ranged param (a tolerance, `(ranged lo
hi)`, `--param name=lo:hi`) is a band over it. A fact whose value is no
number is refused: an operation in it has no result there or may have
none over its range (a division by 0, ln of a negative, a negative to a
fractional power, the square root of a negative), a number in the model
is outside ink's range, or it reads a value ink works out at a run (a
solve, an estimate). Facts and requirements of ratings are data:
loading and `ink check` do not work them out, and a file that
requires no ratings pays nothing for them. `ink check` lists them:
`facts` (each stated fact's name, place, value and unit as read, its
line, and `read`: whether a rating or formula at its place takes that
name) and `takes` (by place, the names taken there). A fact nothing
reads is a warning there, and under each requirement over its place in
`ink require`: a misspelled name changes no answer, so it is said, with
the names that place does take. `ink require FILE` answers
every requirement, of ratings and of claims (Claims), the claims its
limit clauses and relations make among them: yes, no, not settled or
not available, each with why and what it rests on (how: CLAIMS.md,
beside this file), a claim with the prover's witnesses and hints (Claims).

A part may work out facts from other facts about it: its body's `(attr
derives (vec (vec NAME FN UNIT ...) ...))`. The fact NAME, on the part
where it is placed, is FN of the facts its parameters name there
(`(lambda (radial-load speed) ...)` reads those two), each read in the
quantity of the UNIT written for it after FN (`1N 1rpm`; a bare `1`, a
pure number; none for a function of nothing, a number from the part's
record). It is worked out when a requirement asks for it, from each
fact wherever that was stated (the instance that places the part, an
assembly above it, an earlier entry): open, with the reason, where one
is not stated or is stated two ways, and refused where one is in
another quantity than its parameter reads (a speed in Hz, a bare
number). FN is the language's own arithmetic, so the fact has its
units, a band from the bands of what it reads, and its refusals. FN
says where it holds by giving `(vec VALUE OK WHY)`: where OK is not 1
for every value read the fact is refused, and says WHY. FN is called
once, on variables: one that cannot take its units (a sum of a force
and a length) or decides on a fact at load (`if`, where a `select` is
meant) is an error when a requirement is asked. A design does not state
a fact a part works out, and a rating has neither a stated fact's name
nor a worked-out one's at its place: each is an error, since one would
answer for the other.

**`frame-of` paths.** In an assembly's clauses a path starts at one of its
instances placed so far (`b`, `l1/servo`), in its coordinates. Elsewhere
it is a body (its coordinates) or a path from a top-level assembly
(`root/arm/l1`, world coordinates).

**Joints.** MOTION is the connection's one degree of freedom, about or
along an axis of WHERE, applied after the faces meet:

- `(turn z a)` spins about the face normal (a shaft, a hinge);
- `(slide x s)` slides in the face plane (a carriage on a rail);
- `(slide z d)` moves along the normal.

The value is an expression of params and t: a constant is a fixed
setting, a param a slider, `(* 360deg/s t)` motion (an angle for a turn,
a length for a slide: Units). Other fixed offsets go in
WHERE: `(at (move 0 0 -8) (frame-of m shaft))`. One motion per instance.

**Arguments.** An assembly takes arguments as a function does, each with
a default: `(assembly (motor (angle 0)) ...)`. Its clauses read them like
local names. An instance passes values in order, `(instance m (motor (*
3000rpm t)) ...)`, reading params and t but not x, y or z; omitted ones
take their defaults. A part's interface is its argument list: catalog
motors, steppers and servos export `(drive (angle 0))`, rod actuators
`(drive (extension 0))`.

Instances passing equal values to the same body share one build.
Instances of a parameterized assembly also share a build when their
arguments are passed in the same file copy. Across file copies they get
separate builds so each checks its generated claims over the motion at
its argument site. A new build is named after the first instance path for
that combination (`root/m`, so its joint's claim is `limits:root/m/rotor`). `ink
check` lists each one's `params` and what it is `of`.

```lisp fragment
(assembly (motor (angle 0))                   ; a rotor in a housing, driven by angle
  (instance housing housing)
  (instance rotor rotor (frame-of housing axis) base (turn z angle) (speed-limits 0 6000rpm))
  (frame shaft (frame-of rotor shaft)))
(def motor-deg (* 3000rpm t))                 ; t in seconds, 3000rpm in deg/s
(assembly root
  (instance m (motor motor-deg))
  (instance p pulley (frame-of m shaft) hub)  ; rides the shaft
  (instance g gear (at (move 40 0 0)) (turn z (* -0.5 motor-deg))))  ; a 2:1 pair: a multiply
```

**Limits.** `(limits LO HI)`, numbers or params in degrees or mm, adds
the claim `limits:ASSEMBLY/NAME`: LO <= value <= HI, a value claim
checked over the tolerances and the motion (Claims: Bounds and
comparisons), exact unless it states `(limits LO HI (accuracy x))` in the
joint's units. A failure gives the t and params where it fails, and its
`terms` the drive, LO and HI there. A joint outside its limits fails the
claim; it is never a load error and never clamped.

`(speed-limits LO HI)` bounds the joint's rate, `(diff drive t)`, in deg/s
or mm/s: the claim `speed-limits:ASSEMBLY/NAME`, checked the same way. A
motor states its speed range on its rotor's joint, so every assembly
driving it faster fails that claim. An assembly built with arguments has
its limits claimed by the file and line that passed them.

**Published mass.** `(attr mass M)` on an assembly, such as a part's
record mass on its root, is the total mass of the bodies it places, each
counted as often as it is placed. Bodies with a density of their own weigh
what it says; the rest share the density that makes up M, found from
their volume at the params' defaults, or at those `(attr mass-params
(read-json "mass-params.json"))` gives; params the mass expression reads
keep their current values there. So mass scales with volume. Where a body
sharing it deforms (its shape reads t), the mass is published at t = 0:
each body keeps the mass it has there, a deforming one's density its mass
over its volume at each moment, so a spring that shrinks does not pass its
mass to the rigid base beside it. An
envelope body's `(attr mass M)` works the same way. A mass less than the
bodies with their own density weigh is an error; one that gives no body a
density is unused, and `ink check` warns. The derived density is an
estimate (Geometric queries: `ink mass`). M as written is exact, as every
number is; where M is nearly all known bodies a small error in it is a
large one in what is left (a 4.19 g assembly really weighing 4.1909 g
moves a 1.2 mg remainder to 2.1 mg), so write what it is known to:
`(param m 4.19g (tol -0.5mg 0.5mg))` and `(attr mass m)`. A claim reading
the shared density is then asked over the tolerance, the density found
again at each mass (a tolerance on a param a shape reads is not carried:
the claim is refused, as for any solve); `ink mass` and `ink check`
report at M's value.

**Contact.** `(contact A B ZONE)`: A and B are instance names or paths
into sub-assemblies (`left/shelf-l`); ZONE is a field in A's own
coordinates. `interfere` leaves overlap of A with B inside ZONE out of
their pair's proof, as for a screw's thread in a nut: `(contact s1 n1
(shape-of screw thread))`. A zone that is A's own feature (a tag, or an
offset of one) leaves out exactly that feature, at no cost. Any other
zone must hold all of the overlap with margin: use a simple shape clearly
larger than the overlap, clipped to it, `(contact s1 plate (max
(cylinder-z 2 30) (- z 1)))`.

Keep zones tight: a zone hides any clash inside it. It must be proven
bounded at every moment and every value of the ranged params.
`disjoint` and `clearance` claims ignore zones. A contact holds for any
`--target` holding both its ends.

Catalog tapped holes and nuts have thread clearance, so a catalog screw
needs no zone when its thread is in phase with theirs, as it is when
placed on their thread frame (a part's `tap-0` ..., a screw's and a nut's
`seat`): `(instance s screw (at (along-thread h p) (frame-of plate tap-0))
seat)` seats a screw h out from the tap, turned with the thread; in
phase, its engagement is proven clear outright on one cross-section a
pitch tall (two threads of one pitch on one axis are the same at every
turn). Threads placed by hand can be out of phase and overlap (an M3
screw in its nut, turned 180 degrees, 0.08 mm deep). Faces seated by
CFRAME need no zone: touching faces are proven to the accuracy.

**Relations.** `(relate A B KIND)` says what A and B (named as in
`contact`) are to each other: claims and data, checked where it is used,
and a zone where the parts are meant to overlap. Gears in mesh are
related by std's `gear-mesh`, `bevel-mesh` or `rack-mesh` (Standard
library, Paths and power transmission). A mesh has no zone, so the teeth
are proven with the rest of the pair, over the pinion's angle (Motion).

```lisp fragment
(instance w wheel (at (move (mesh-distance 15 45 1) 0 5)) (turn z (mesh-turn 15 45 0 motor-deg)))
(relate p w (gear-mesh 15 45 1 20 0.2))          ; std: teeth, module [, pressure angle [, backlash]]
(relate p (g (at (rot-y 90))) (gear-mesh 16 50 1.27 14.5))  ; g's axis is its own x
```

KIND is a function written in `.ink` (std's or a part library's); ink
knows no kinds. It is called with A's and B's placements in this
assembly's coordinates and their shapes, each in its own coordinates;
one of six arguments also gets their bodies' attributes, a record each
(an end that is an assembly is refused: relate its body by its path),
as `(get attrs "gear" (vec))` reads them: numbers at load (a param's
value: a run ranging it refuses the relation's claims), strings, lists
and data; any other attribute is an error. A number that is not known at
load (a colour painted by x y z, a density solved from a published mass)
is text there saying why, so a kind reading it as a number is refused and
one reading only other keys is not. An
end written `(PATH FRAME)` is read at FRAME, a frame its part declares or
one written there. A relation never places or moves anything. The kind
returns `(relation NAME [(zone WHAT ZA [ZB])] (claim CNAME VALUE LO HI
[(accuracy x)]) ... [(data KEY VALUE ...)] [CLAUSE ...])`:

- zone: overlap inside both ZA (in A's coordinates) and ZB (in B's;
  omitted, everywhere) is intended, as a contact's. WHAT says what it
  leaves unchecked.
- claim: the value claim `relate:ASSEMBLY/A~B/CNAME`, LO <= VALUE <= HI,
  checked like a joint's limits, exactly or to `(accuracy x)`. VALUE
  reads numbers, params, t and the frames.
- data: shown with the pair.
- CLAUSE: a package's (Data for other packages).

`ink check` lists an assembly's relations; `interfere` names them in the
pair's row with what a zone left unproven (`"skipped"`).

**Held parts.** Nothing checks that a part is held: a part touching
nothing passes check, require and interfere. Where holding matters, place
it on the frame of what holds it and claim it, `(touches a b)`.

**Interference.** `ink interfere` claims that no two instances overlap
outside declared contact zones, each pair a `disjoint` claim to its
parts' accuracy as a claim's (Claims), or to `--accuracy x`, the
answer's `accuracy`. A pair passing only by that allowance says how deep
it came (`held`, `near`). Without root it checks every body, unless an
assembly places one: then it refuses; name one with `--target`.

`pairs` lists the pairs whose bounding boxes overlap, those failing
first, then those unknown, each row the answer of its `disjoint` claim
(Results) with the contacts it used. The
summary counts `overlapping_pairs` (`overlapping_pairs_estimated` of them
estimates), `near_pairs` (holding, with a `near` depth above 0) and
`undecided_pairs`. Pairs placed alike (every copy of a
sub-assembly) are asked once.

`ink check` lists each body's and assembly's frames (origin and axes at
the current params and t) and each instance's joint (kind, axis, drive as
written and its value, limits).

#### Paths

Paths name placed parts: `root/left/bracket`.

- `(shape-of path)`: the world-space field of a body, an assembly (the
  union of its parts) or an instance path, for claims and other shapes.
- `(shape-of path name)`: a named construction inside it (a def or tag
  of the body), placed where the instance is: `(shape-of root/bearing
  inner-ring)`.
- `(shape-of path (except P ...))`: the union of its parts but those at
  or under each P, a path from the same top, such as a tool's
  surroundings but the screw it turns. A `for` may write the paths:
  `(except (for q (vec "a" "b") root/$q))`.

#### A part on a tilted face

1. Write a frame on the face, z out of it, as `(axes ...)` with the tilt
   in its directions (as ink prints frames).
2. Build the body's slanted piece flat in that frame and put it there
   with `place`, so face and frame cannot disagree.
3. Meet the part's own frame with it (CFRAME).

From examples/phone-stand.ink, a phone leaning back at 65 degrees:

```lisp fragment
(def phone-height 149.61)                          ; the phone's body: frames back, front
(def a (deg->rad 65))
(def face (axes (origin 0 0 8) (x 1 0 0) (z 0 (neg (sin a)) (cos a))))  ; its y runs up the plate
(def plate (min (move (box 80 100 5) 0 50 -2.5)    ; in the face's frame, the face at z = 0
                (move (box 80 6 12) 0 -3 1)))       ; a lip below the phone
(body stand
  (shape (min (move (box 90 110 6) 0 15 3) (place plate face)))
  (frame seat (at (move 0 (+ (/ phone-height 2) 3) 0) face)))
(assembly root
  (instance stand stand)
  (instance phone phone (frame-of stand seat) back))   ; back on the face, x with x
(require apart (disjoint (shape-of root/stand) (shape-of root/phone)))
```

`(at (move ...) face)` moves along the face's own axes (up the plate).

#### Placement by example

Claim where each part must be, as these examples/ do:

- motors-on-radius.ink: motors on pads round a circle (Values: turn
  first, then move). Each motor `touches` the frame, and a probe along a tap's frame, `(place (probe 0 0 2) (frame-of
  root/motor-0 tap-0))`, is `disjoint` from the frame where the holes
  line up.
- two-link-arm.ink: servos met face to face on seats; value claims hold
  the tip frame's origin where the angles put it.
  ```lisp
  (instance e (servo elbow) (frame-of upper seat) mount)
  (require elbow-on-link (touches (shape-of root/e) (shape-of root/upper)))
  ```
- frames.ink: a bearing seated by its face on a shoulder `touches` the
  housing below the seat (failing afloat or sunk) and is `within` the
  bore.
- sliding-drawer.ink: a box rides a slide's moving frame; the other slide
  is the same part met to a seat facing the other way.

#### Tool access

Tool reach and swing are claims about the tool placed on its fastener: a
catalog fastener's `drive` frame (z toward the tool) met by a hand tool's
`engage` frame:

```lisp fragment
(def out (ranged 0 50))                     ; withdrawn 0..50 mm along the axis
(def swing (ranged 0 60))                   ; turned through 60 degrees
(assembly root
  (instance axis axis.root)
  (instance key key.part (at (move 0 0 out) (frame-of axis/rail-screw-0 drive)) engage-short (turn z 60)))
(assembly turning
  (instance axis axis.root)
  (instance key key.part (frame-of axis/rail-screw-0 drive) engage-short (turn z (+ 60 swing))))
(require reach (disjoint (shape-of root/key) (shape-of root/axis)))
(require swing (clearance (shape-of turning/key)
  (shape-of turning/axis (except turning/axis/rail-screw-0)) 0.05))
```

Reach is the tool withdrawn along the socket's axis; swing the tool
turned, against every part but its fastener. The hand is not checked.

#### Assembly order

A build order is checked step by step: each part comes in from outside
the parts placed before it, along a straight path clear of them. The
order is a list of std's

- `(build-step PART FROM [INTO ...])`: the part (a path under the
  assembly, as text), the direction it comes in from (a unit vec), and
  the parts its thread turns into;
- `(build-laid PART)`: a part laid in place by hand (a belt), checked
  only where it lies.

The claims are one loop over the order: copy it from
examples/build-order.ink. Per step it claims the part, moved out along its
path by every `out` in 0..1 (`build-along`, `build-length`), is disjoint
from the parts placed before it but those it turns into (`build-later`,
`build-left-out`), and disjoint from each of those but std's
`(thread-band d p)` about the screw's axis (`build-into`), so a screw
longer than its blind tap fails.

A failing step gives `out` and the point where its part meets one
placed before it. The order is for an assembly at rest. ink checks the
order and directions written; it searches for neither. A part that turns
into place, or one joined off the assembly, is not checked as such.

## Claims

`(require NAME CLAIM [(accuracy x)] [(res r)])` states a requirement; no
two share a name. (Before 0.2 the word was `assert`: a file that still
writes it is refused, with the form to write.) `ink require` decides it by the cheapest method that does: a witness
point, an interval proof refined only where undecided, or a measurement.

```
(empty f)                     f has no interior
(disjoint a b)                a and b do not overlap
(subset a b)                  a lies inside b
(within f (x0 y0 z0) (x1 y1 z1))
                              f lies inside the box (corners: constants and params, each
                              written (x y z) or a vec of three, as (nth (bounds g) 0))
(clearance a b d)             no point of b within d of a; (clearance a b 0) is (disjoint a b)
(touches a b)                 a and b sit against each other
(overlaps a b)                a and b overlap
(range f lo hi region)        lo <= f <= hi everywhere inside region
(value expr < v)              a number from params and t (also <=, >, >=)
(volume f < v)                the integral of 1 over f (any comparison)
(integral target g < v)       the integral of g over target's bodies, summed
(mean target h [(weight w)] < v)   the (w-weighted) mean of h over target
(METRIC FORM < v)             a package's metric: (min-thickness part > 1), (cavities part < 1)
```

A claim is a form written in a `require`, not a value, so a function
cannot return one; it may compute what a claim reads (its shapes, points
and bounds). `for` writes many claims: `(for p (vec (vec 1 2 3) (vec 4 5
6)) (require solid-$p (subset (probe (nth p 0) (nth p 1) (nth p 2)) B)))`.

**Accuracy.** A claim's accuracy is how far its answer may be from the
truth. For a claim about shapes it is a length: the claim's `(accuracy
x)`, else the finest of the parts its shapes are made of, each its
`(attr tolerance x)` (how far the made part may differ from it, a
constant), or 0.05 mm where it declares none. A shape is made of every
body it reaches at any depth, through unions, assemblies, placements,
warps and other bodies' shapes; a named piece (`(shape-of p tip)`) is of
its body; a shape of no body (a box written in the claim) counts as
declaring none. A part printed to 0.2 against a shaft machined to 0.01
is checked to 0.01, two printed parts to 0.2, a part declaring none
against either to 0.05 or finer. `accuracy_given_by` says which
(`require`, `--accuracy`, or `tolerance` with each part's `tolerances`).
A param's tolerance is its band (Tolerances). `holds` and `fails` are
shown under the claim's accuracy rule; `unknown` means ink could show
neither, and is never a wrong `holds` or `fails`. The rule, by claim:

- `(clearance a b d)`, d above 0, is exact on d. It holds where no point
  of b is within d of a, and fails on a pair of points closer than d. The
  accuracy is how fine the search goes, so a gap within about it of d may
  be `unknown`: at 0.05, `(clearance a b 0.09)` fails with a gap of 0.06,
  holds with 0.2, and may hold or be unknown at 0.095.
- `disjoint`, `subset`, `within`, `empty` and `touches`' overlap allow a
  violation no deeper than the accuracy. Parts touching, or overlapping by
  less, are disjoint; across two flat faces that is an interference up to
  twice the accuracy thick. A pin may stand proud of its bore's face by
  up to twice it and be a subset. So a `within` shows no gap (exact on
  its corners, a part flush with them would be unknown); claim a
  clearance for one.
- `(overlaps a b)` (shapes that must meet, as in a weld) holds by a
  witness point whose cube of half-side the accuracy lies inside both.
  It fails only where no point is inside both: a press fit 0.04 deep at
  0.05 is `unknown`; faces that only touch fail (claim `touches`).

A claim is asked once, at its accuracy. So a shape no thicker than twice
the accuracy can pass as absent: at 0.05, `(empty ring)` holds for a ring
0.04 thick, and `(overlaps f (box 0.02 0.02 0.02))` cannot hold, since no
point of the box is that deep (the `hint` says so). State a smaller
`(accuracy x)`, under half the thinnest part's thickness, or check points
(`probe`). The accuracy is a length whatever the fields' scale.

**Points.** `(probe x y z)` is a point (constants and params). Check a
dimension with a probe pair, one point that must be material and one
that must be clear: `(subset (probe x y z) part)` and `(disjoint part
(probe x2 y z))`.

A point is decided exactly (`"accuracy": 0`): inside where the field is
below 0 there, outside where above. Where it is 0, inside where the shape
is on every side of it (faces placed flush), outside where air is beside
it (a point on a surface), read within a few billionths of its size.
Its `value` is the field's bounds there; it is `unknown` only within
rounding of the surface, or at a 0 neither reading shows.

`(place (probe x y z) F)` is the point at (x, y, z) in frame F. A `warp`
of a point (std's `move` and `rotate-z` are warps), a min of points or
arithmetic on one is an error: use `place`.

**The other claims.**

- `(touches a b)`: parts that sit against each other, a nut on its
  washer. It holds where they overlap no deeper than the accuracy and
  come within it of each other (`contact`), and fails where they overlap
  deeper or float off. Over ranged params it and `overlaps`
  only fail or stay `unknown` (Tolerances).
- `(range f lo hi region)` bounds a field. Slopes and normals (a moulded
  part's draft) are a `range` of `(diff f z)` or a `grad` component,
  right where the field is an exact distance. Overhangs are read by signs
  instead, for any field (Manufacturing checks).
- `(value expr < v)`: any number from params and t (ratios, design
  rules, solved values), exact unless it states an accuracy.
- `(integral target g < v)` and `(mean target h < v)`: mass properties
  and any other integral over material. The target is a body, an
  assembly, an instance path (its bodies, placed) or a field. In g, h and
  w, `(attr NAME)` is the attribute of the body being integrated. h may be
  a `vec` of fields, each one's mean claimed in one pass: a centre of mass
  over a convex footprint is one claim, one field per side
  (examples/desk-lamp.ink).

  ```lisp
  (require light (integral root (attr density) < 2500))           ; mass
  (require low (mean root z (weight (attr density)) < 80))        ; centre of mass height
  (require over-foot (mean root (vec (- x 40) (- -40 x) (- y 30) (- -30 y))
                       (weight (attr density)) < 0))            ; over a 80 x 60 footprint
  ```

  Integrals, means and volumes are estimates (`"accuracy": null`),
  measured as `ink mass` measures a mass (Geometric queries), or at cells
  `(res r)` sets. Over ranged params they are estimated at evenly spaced values
  (`sampled_values`) and fail where one fails (`params`). A body not
  proven bounded is not integrated: `unknown`, `unbounded`.
- `(METRIC FORM < v)`: a package's metric of one of its forms, or of a
  shape: `(min-thickness part > t)`, `(cavities part < n)` (Manufacturing
  checks). The form must come before the claim. The package measures it
  by running the form, so the claim is an estimate, and `estimate` says
  what it rests on: its cells and its grid error
  (`grid_error`, the number) over three grids, each twice as coarse as the last, the coarsest the
  coarsest of 4, 2 and 1 times its cells that reads the part (so a part
  too thin for coarser cells is checked, and answered, on finer ones:
  `res` says which). One that holds or fails by less than its grid error
(none where the three move apart the same way by more than one grid's
reading is known to: a value moving away as the cells shrink may go on
moving) and what its
read may miss where its package states it (a peak stress's, SIM.md),
or that
  no such grids read, is `unknown`; so is one whose run failed, with the
  reason in `unavailable`. A metric
  measured with a band (where the value may lie, an estimate) is answered on the band widened by
  that grid error on each side, whether it holds or not: it holds only
  where both ends hold and fails only where both fail; else `unknown`,
  naming the part that dominates.

**Bounds and comparisons.** A claim's bound (the `v`) is any expression
of constants and params, at the current `--param` values (solved params
at their solved values): `(min-thickness part > (* 2 nozzle))`. It may
not read x, y, z or t, and only a value or mean claim's may read a ranged
param.

A comparison (`value`, `volume`, `integral`, `mean`, a metric) may state
an accuracy in its own units: how far past its bound it may be (a void
count's is a length on its shape: Manufacturing checks). `(require
reach (value angle < 90) (accuracy 0.01))` holds where `angle - 0.01 <
90` everywhere. Generated value claims (`limits`, `speed-limits`, a
relation's `claim`) take one too. An exact value claim that meets its bound (a joint at
its limit) over a motion or a ranged param may stay `unknown`; with an
accuracy it decides. A value claim on constants is decided in real
numbers, not by their floats: a decimal and a unit written as one are
exact fractions, and their sums, differences, products and quotients
stay exact (`(value (/ 1ft 12in) > 1)` is a proven no, `(value (+ 0.1
0.2) <= 0.3)` a proven yes); anything else (pi, a degree, a sine, a
square root) is the interval its rounding cannot leave, and a
comparison within it is `unknown` (`(value (sin 30deg) >= 0.5)`, and
`(value (* (sqrt 2) (sqrt 2)) > 2)`: an interval bounds, it does not
know that sin 30 deg is 1/2). Everywhere else a constant is its float: in
a shape, and where the model needs a number when it loads (a loop's
count, an `if`, a tolerance's or a range's ends, a case). A number
computed by the program, not written in the file, is its float too. A
value claim reading a param takes it at its one value. A param written
in digits (its default or a case) is the fraction they write, as a
literal is, and so are the numbers written beside it: `(- (+ big 1)
big)` with big 1e16 is exactly 1, not a float's 0, and 334 mm is
exactly twenty links of 16.7. One computed through something inexact (a square root) is
its float, its arithmetic rounded outward as a range's is: not settled
there.
Points and ranges are exact and refuse `(accuracy x)`.

### Results

`ink require` gives each requirement a row: its `verdict` (`yes`, `no`,
`not settled`, `not available`; `why`, for the last two), `strength`,
`covers` (over how much of each variable it reads it was shown), `from`
(`file:line`), `notes` and `left_out` (how: CLAIMS.md). A claim's row
carries the prover's evidence beside them: its `status` (`holds`,
`fails`, `unknown`) and `accuracy`:

- a number: guaranteed to it, under the claim's meaning (0: exactly).
  Points, ranges and value claims are exact by default.
- `null`: `unknown`, or an estimate with no guaranteed bound: package
  metrics, integrals over sampled values of a ranged param (an integral,
  mean or volume claim is otherwise decided by its proven bounds, and
  exact to them), a claim measured at moments of a motion, and any
  claim that reads a sampled field or a solved value (not a certified
  root), directly or through defs, bodies and placements.

Its other fields say what the answer rests on:

| field | when | meaning |
|---|---|---|
| `witness`, `witness_pair` | a shape claim fails, except by being apart (a failing `overlaps`, or a `touches` too far off, has none) | the violating point (or pair). For an emptiness claim (the `depth` row's), every point within the accuracy of `witness` violates the claim. A failing clearance guarantees only its pair: two points, one in each shape, `gap` apart. A probe's is `point` |
| `at_witness` | a shape claim or a probe fails at a point | for each shape it compares (its first, then its second), the placed bodies there, each by the path of the instance placing it: those whose shape holds the witness (a clearance: the nearest to each point of its pair). Which bodies meet, where the point itself moves whenever the proof does: compare an expected failure by it |
| `value` | fails; a number claim | what the claim measured: a clearance's `gap`; a shape claim's field at the witness, negative (`-depth`); a value, volume, metric, integral or mean claim's number (a mass is an integral). Of several fields, each is under `integrals` or `means`, with its status |
| `gap` | a clearance fails | the distance between `witness_pair`: the first pair found closer than the clearance, so the smallest gap may be smaller (`ink distance` measures it) |
| `depth` | an emptiness claim fails (`empty`, `disjoint`, `subset`, `range`, a clearance of 0) | how deep `witness` lies in what violates the claim (inside both shapes for disjoint; inside a and outside b for subset), as the fields measure distance, walked down the field as deep as it leads: the overlap is at least this deep (two blocks overlapping 2.5 mm: 1.25) |
| `held`, `near` | an emptiness claim holds, guaranteed (over a motion, with t ranged; an `interfere` pair as placed, with its zones) | `held`: `outright` (nothing violates it) or `to the accuracy` (some cells held only by its allowance). Then `near`, where measured: `depth`, the deepest violating point found in those cells, by bounds of the field as placed there at one moment and value of the ranged params (to a µm: seated faces read 0), and where above 0 `at` and `params`, those values (t's over a motion) |
| `params`, `t` | fails over ranged params or a motion | the params' values and the moment there |
| `hint` | unknown, where it can say; measured over a motion | what would decide it: a smaller `(accuracy x)`, a finer `(res r)`, a probe. Over a motion, why ranging t did not, and that search's hint (`motion.t_ranged`) |
| `undecided_cells`, `undecided_volume`, `undecided_cell` | unknown, by a search that kept them | how many cells it left undecided, their total volume, and one of them (the cell that ended the proof, where one did): `centre` and `size` |
| `hunt` | unknown, after a hunt | the cells (a value's pieces) searched for a violation once the proof ended, none found (`cells`), and what stopped it (`stopped_by`): `cap` (as many as the proof tested), `budget`, `time limit`, or `refinement limit` (every cell left at it) |
| `accuracy_asked` | unknown | the accuracy it was undecided at |
| `accuracy_attained` | unknown shape claim, where found | the accuracy it does hold to: `(clearance a b 0.29)` unknown with 0.003 means `(clearance a b 0.287)` holds. `accuracy_attained_in_window` where only the window was examined |
| `unbounded`, `in_window` | a shape not proven bounded | only the window was examined, and what was settled there |
| `res`; `res_given`, `res_given_by` | a claim searched on cells (not a value or a probe); where a res is given | the finest cells the search used, coarse where their bounds decided it before refining; the cell size asked, and by `require` (the claim's own) or `--res` |
| `accuracy_given_by`, `tolerances` | the accuracy is not the default | `require` (its own), `--accuracy` (interfere's), or `tolerance`: its parts' (Claims) and each one's |
| `band` | a metric measured with one | the interval the claim was decided on: `low`, `high` (widened by `grid_error` where one was measured), each part's half-width (`half_widths`, `grid` among them), the largest `dominant`, and `basis` |
| `samples` | reads a sampled field | each field's source, where, spacing, its error bound where one can be derived (`grid-bake`'s; none for `import-stl`, `redistance` or a package form's result) and the error measured at probe points (`measured_error`) |
| `solved_inputs`, `unsettled` | reads a solved value | the solves it rests on; those short of their accuracy, and why |
| `param_spans` | reads a ranged param | the params and their ranges; not a `(ranged)` no def names (`hidden`'s views) |
| `texts` | the model has text params | the answer's, at the top: each one's value, which every claim was asked at |
| `cases` | a param declares cases (Top-level forms: Cases) | in a row, each case's answer, `case` its values; at the top, the combinations asked |
| `motion` | over a span | how the motion was checked (Motion) |
| `unavailable` | unknown, not checked | why: a solve with no value, or a model outside ink's range (docs/CORE.md Answers; interfere, distance and bounds refuse there) |
| `declared_at` | placed motor effort claim | `file:line` of the underlying rating |
| `ms`, `cached` | always | wall time; whether the proof was reused from disk or from an identical generated check in this run |

The summary gives `all_yes` (every requirement yes, and none of them
shown at one value only of what it reads: `yes_at_one_value` names
those; and no requirement of ratings leaving a part its model loads
unrated, a stated load unapplied or the rest of a motion unread:
`parts_not_rated`, `parts_with_no_limit_stated`, `loads_not_applied`,
`motion_goes_on`, as
docs/CLAIMS.md says), `all_proven` (and each proven, none an estimate or only as
stated) and `counts` by verdict. A declared limit that a package could not check gets an
`unknown` row and appears in `not_checked`, with the reason in
`not_checked_note`. When declared cases are run, each default row keeps
its original status and accuracy in `at_values_result` alongside the
overall verdict and the individual case results.
The list covers the requested claims across all allowed cases. Any such
limit, or a case the time limit prevented loading, makes both summary
fields false and the exit status 3. An empty request gives null for both
fields (`nothing_required`) and exits 0. All yes (including estimates)
exits 0, one that is not exits 3, and an error exits 1 (`ok` false).
`ink require FILE` asks what FILE itself requires, including every generated
`effort-limits:root/...` for an actuator placed under its root. This includes
placements without arguments and roots inherited from imports. Each uses the
run's gravity, loads, motion and parameter cases. Unchecked limits under that
root are reported as unknown. The claims its imports require are counted
in `imported_not_run`, and asked too with `--imports`. `--name NAME` asks
one requirement.

`ink check` lists what the file defines, each def with its value
(shortened where long), and each import in brief (`--imports`: all they
define). `has_cases: true` means a parameter in the loaded design declares
cases, including in an import. `--name NAME` gives one def's value whole: to see a value, def
it. Warnings cover every body, imported or
not. `ink check --load` skips surveying the bodies' surfaces: a quick
look at a large model.

Each claim `ink check` lists has `cost_estimate`, what `ink require` will
spend on it, said before it runs: `basis`, the work counted (cells, and
the instructions an evaluation runs there; each case it is asked at), and
`class` (`under a second`, `seconds`, `tens of seconds`, `minutes`,
`over ten minutes`, on this machine's threads) where what that work costs is
measured: a body's thickness (a ray from a surface sample per cell of
its res, on three grids; over a motion, at 17 moments), a value at one
point, a refusal. That measure (benchmarks/cost.md) puts 46 of 58 real
thicknesses in their class and none two off; it was fitted on those same
claims, on each one's first case; none took minutes, so the classes from
`minutes` up are extrapolated. `class` is null where the work is counted
but its time is not known: a proof (its cells are where its bounds
cannot yet tell, which no cheap look predicts: seated faces settle in a
few, a curve on a curve at about the accuracy; its basis gives its
budget and its field's size), a value over ranged params or the motion
(a search), an integral, and a thickness of a shape no body is (check
samples only the bodies' surfaces). The key is left out where nothing is
estimated (a package's solve, a shape not proven bounded) and with
`--load`. A thickness estimated at minutes is cheaper at a coarser
`(res r)`.

**The cache.** The kernel's answers (proofs, searches, masses, walls'
thickness, meshes, solves, `ink check`'s survey) are cached in
`~/.cache/ink`, keyed by the question as the kernel sees it and a hash of
the ink binary: an unchanged question is answered from the cache whatever
else changed (`cached`); a rebuilt binary recomputes everything.
It is kept to 8 GB, the least recently used out first; `INK_CACHE_MB=n`
keeps it to n megabytes. `INK_NO_CACHE=1` disables it.

**Time limits.** `--timeout S` (seconds, or a time such as `2min`) stops
searches after `S` seconds so a long command (`ink interfere`, `ink
require`, `ink mesh`, ...) answers with what finished, marked `stopped`
(on `ink interfere`, `holding_pairs_estimated` counts moving pairs whose
yeses stayed measured rather than upgraded to proven over ranged `t`).
The same call with a larger `--timeout` goes on from the cache.

**`ink mcp`.** The MCP server offers every command but `mcp` as two tools:
`ink_write` runs a call that writes a file (mesh, share, or a file option
such as `out`), `ink` any other, and each refuses the other's.
A call stopped at its `timeout_s` answers with what finished, marked
`stopped`; the same call with a larger `timeout_s` goes on from the
cache. A load it stops (sampling a grid: `redistance`, `grid-bake`) is
an error there, marked `stopped`. A `min-thickness` under a time limit
measures its grid-error cells (4h, 2h) before h, so one stopped answers
on the finest that finished (`res`; `stopped_before_res` is h) with their
`grid_error` from the next coarser, holding only where that error cannot
cross the bound (a `>` failing on its value alone: Walls). Of a `require` or `interfere`,
`finished.left_estimated` lists up to five claims or pairs left (running,
at any case, or not begun), costliest first, each with its
`cost_estimate`, and the message gives the one running's. A call's
estimates never cost it time: a thickness whose body `ink check` has not
surveyed (cached) has class null, "no survey yet; ink check computes it".
An answer too long for clients comes back shortened, still whole JSON;
`shortened` names the file holding all of it.

### Resolution

Without `--res`, a search halves its cells down to half the claim's
accuracy: 0.025 mm at the default 0.05, whatever the model's size. Faces
that touch by design (a lid on its base) are proven in a few cells,
whatever their area; curved contacts (teeth in mesh) at cells of about
the accuracy. A cell undecided at the refinement limit (space at the
finest cells; a ranged param halved 8 times in a shape claim's cells, 32
in a value claim's pieces) ends the proof: the claim is not proven at that
refinement; so does a search that would test over four million cells.
The search then hunts a violation in the rest, the largest cells first,
until it has tested as many cells again as the proof did: one found fails
the claim; else it is `unknown`, with `undecided_cell`,
`accuracy_attained`, and `hunt`. A hunt that finds no violation does not
show there is none: the cells it did not reach may hold one.

- `--res r` (to `require`, `interfere`, `distance`) lets the search halve
  undecided cells down to r, if r is finer than half the accuracy. Cells
  settled earlier need no refinement. A coarse `--res` costs no answer.
- `(res r)` after a claim does the same for that claim, or the run's
  `--res` where finer. It is for claims checked on cells; a package
  form's cells are its own `(res ...)`.
- Small things (MEMS, optics, a fit in hundredths) ask a finer
  `(accuracy 0.0005)`.
- A finer res decides a cell whose bounds narrow with its size. It does
  not decide a claim at its boundary (its value reaching the bound), nor a
  slope over a cell where a min, max, abs or clamp switches, nor across a
  `floor` or atan2's cut, unless the bounds see the field continuous
  there: an angle read by its sine or cosine, turns (degrees / 360) less
  their floor, a wrap read by its distance from 1/2 (std's thread-iso,
  a belt's teeth), however many of its steps the cell crosses. A
  `select` whose condition reads only t or params takes one branch over
  the cell at each moment.
- Elsewhere `--res` is the sample size, used as given: `ink mesh`,
  `render`, `slice` and `components` sample at it, `ink mass` measures
  at it.

### Where analysis looks

Claims and queries cover all of space: each shape's bounds come from a
search that cannot miss material, so a part placed 20 km away is found.
A shape the search cannot prove bounded (a half-space, an unclipped
lattice, a field interval arithmetic cannot bound at infinity) is
examined in the window, the 1 m cube centred on the origin. A result
that needs all of the shape is then `unknown`, `unbounded`, with what was
found there (`in_window`):

- emptiness, subset, `within`, `range` and clearance proofs, and
  `overlaps` disproofs, need all of it; a witness does not (a point
  inside disproves emptiness, proves an overlap);
- volumes, integrals, means and shape metrics are `unknown`.

A 2D profile (a field of x and y alone) is a prism, the same at every
height: a claim of profiles alone that is the same at every height
(emptiness, `subset`, `disjoint`, `overlaps`, `touches`, `range`,
clearance, `min-thickness`, `cavities`: none) is asked in the plane, of
the profiles' extent in x and y, so `(subset letters face)` and
`(min-thickness letters > 0.8)` of 2D letters cost a search of the plane.
One shape read with z (an extrusion) makes it a search of space; an
emptiness proof (emptiness, `subset`, `disjoint`, `overlaps`) or a
clearance of extrusions costs about the plane's inside their spans, a
`min-thickness` a few times the plane's.

The claim's own field must be bounded: `(empty (max x ball))`, a
half-space cut by a ball, proves outright. A union is bounded by its
members' bounds, an intersection by its cheapest member's. To limit a
claim to a region, intersect with it.

A z twist preserves height and distance from its axis. Ink can use a
cylinder enclosing the source to prove the twisted part fits the window.
If the cylinder is too wide, a tighter containing shape can help.

On a shape not proven bounded, `ink bounds` and `ink mass` say
`unbounded`; `ink distance` proves a lower bound only where a target is
bounded; `ink interfere` checks it against every other item
(`unbounded_items`); `ink ray` searches only to the window's edge;
load-time `(bounds f)`, `grid-bake` and `redistance` are errors.
`ink mesh`, `render`, `slice` and `components` refuse the command if
any selected target lacks a proven bound. A closed surface or outline
at the edge of the window can hide missing material, and a sealed cavity
can lie entirely outside it. These commands exit with an error before
writing a file or reporting whole-part counts. To specify a finite part,
intersect its field with a bounded shape that encloses the intended design.

"Proven" holds for real numbers: bounds round outward, so rounding cannot
fake a proof or a witness (CORE.md).

## Names

- Every `def` of a spatial field (one reading x, y or z), every body and
  every `(tag name f)` names a construction. The name is a symbol, a
  string, or an expression giving a string (`(tag (concat "rim-" n) f)`).
  `(tag name f opaque)` hides the names inside it (a blend is one piece).
- std's blends are named: `fillet-union` and `fillet-intersect` name
  their rounded edge `fillet`, `chamfer-union` its bevel `chamfer`,
  `smooth-union` its blend `blend`.
- At a surface point, the kernel follows the min/max/select branches that
  decide the value and reports the names passed (`ink ray`).
- `ink check` warns about a named construction proven to form none of a
  body's surface, and about one it never saw deciding a patch of it where
  it searched (undecided: it may form only an edge, where it ties with
  what covers it, or a blend's weighted-out side, or a face the sampled
  search missed; the warning says how many cells it searched). A body
  that moves (reads t) is judged over the t-span: hidden is proven at
  every moment, and a name is seen if seen at any of up to 16 moments
  (the span halved again and again). Another body used in this one
  (`shape-of`) counts as one construction.
- `(shape-of path name)` gives a named construction inside a placed part
  (Assemblies).
- `ink check` also warns about a body with no interior: proven empty, or
  a sheet of no thickness, as from intersecting two shapes that only
  touch (`(max a b)` written for a union). `ink render` and `ink mesh`
  warn about a target that meshed to no triangles.
- Mathematically identical constructions share names.

## Geometric queries

Each but `ink distance` (which takes `--a` and `--b`) takes `--target`:
a body, assembly, def or instance path. The default is root, else every
body (`ink render`: else the first assembly none places). `ink help`
gives each command's options.

- `ink bounds`: each target's proven box (`min`, `max`, `size`) and the
  size of the cells its search ended on (`cell`).
- `ink eval --points '[[x,y,z], ...]'`: at each point, the field's value,
  inside, gradient and rates (along t and each param), the undefined
  cases met (`undefined`), and attributes `--keys` names (`attrs.color`,
  `attrs.spec.mechanical`).
- `ink ray --from '[x,y,z]' --dir '[dx,dy,dz]'` (or `--rays '[[[o],[d]],
  ...]'`, many in parallel): where each ray first enters material, none
  skipped: the point, distance, normal (null where the field has no
  gradient) and the name there (`name`, the innermost; `path`, outermost
  first). A ray that only touches or runs along a surface does not hit
  it, but where only rounding decides which side of a surface it runs
  (a face turned at an angle, the ray in it) its hit is where that
  begins, with `rounding`: how far from 0 the field may be along it. One
  starting inside material answers its origin (`starts_inside`),
  so start outside the material and inside the window (Where analysis
  looks). `--points`, `--from`,
  `--dir` and `--rays` also take `@file.json`, or `-` for stdin.
- `ink distance --a T --b T`: the minimum distance as a bracket, to the
  model's accuracy, by the clearance search. `upper_witness` is the
  distance between a witness pair of points, each confirmed inside its
  target, with the names there (`name_a`, `name_b`); null where no pair
  is confirmed (a target with no inside). `lower_proven` is what the
  search proves: 0 where neither target is proven bounded, or they read
  inexact inputs. Overlapping targets report distance 0 with an overlap
  witness.
- `ink slice [--axis x|y|z] [--at VALUE] [--accuracy MM]`: the loops
  where the plane (default: through the middle of all targets) cuts each
  target, each with its kind (`outer`, counter-clockwise in `plane_axes`,
  or `hole`), area, perimeter and bounds. `--keys kind,points` gives their
  points, simplified to `--accuracy` (default half the cell). `--out
  FILE.dxf` writes the loops as an outline to cut (DXF R12, mm, a layer
  per target), and the targets' marks (kind `mark`, Bodies) dashed on
  their own layers. A sharp corner is cut across by up to a cell.
  `unresolved` (where there are any): samples reading exactly 0 that
  nothing near them decides (no surface crossing them, no air beside
  them by bounds, and not material all round them on a 0 of no volume,
  where pieces meet flush, drawn as material); drawn as air.
- `ink mass [--density NAME|NUMBER]`: each target's mass, volume and
  centre of mass; with `--keys` naming them, its `inertia` tensor about
  the centre of mass, `principal_moments` and `principal_axes`. `total`
  sums the targets.
  - The density is the bodies' `density` attribute (Bodies) unless
    `--density` names another attribute or gives a number.
  - A target with no density has null `mass` and `center_of_mass`, its
    volume measured; `total.without_mass` names such bodies, `total.no_mass`
    what would give them one.
  - Mass, volume and a claim's integrals are estimates within proven
    bounds (`mass_bounds`, the volume's `bounds`; none for a target not
    proven bounded, whose mass is the window's). For a published density,
    the mass bounds include the uncertainty in the known bodies' mass and
    in the shared volume. When subtraction leaves little mass, the known
    bodies are refined toward 2% of the remainder, within the integration
    budget. A quantity depending on a measurement with no
    enclosure has no proven bounds. Integration uses boxes of the body's
    box, the widest bounds first, each halved along the axis the field
    changes most along, a box the surface or a density's jump crosses
    bounded by that field's slope over it (a flat face exactly, however
    thin the wall), until the mass (a claim's: each part's first
    integral) is proven within 2% of its estimate; `res` the smallest box
    side. Where the
    rates the bounds narrow at project a proof over 2^18 boxes twice
    running (a field that jumps, edges everywhere), three times where the
    rate is still rising (features smaller than its boxes), it is given up
    after about 32 or 64 thousand, and the mass estimated on cells, a 100th of the body's box
    and halved until a halving moves it under 1% (`sampled`; the bounds
    the proof reached). The mass and volume given are estimates, proven to
    lie within their bounds (`mass_bounds`, the volume's `bounds`), which
    are the only error stated. One that did neither says so (`unsettled`),
    as does every answer reading a published mass's density that did not.
    A claim measures to 32%, 8%, then 2%, and stops where its bounds
    decide it; one they do not decide is `unknown`, whatever its estimate
    (shown as its value) says.
    `--res` halves every box whose bounds differ to sides of at most it
    instead.
- `ink components`: each target's cavities, as the `cavities` metric
  counts them (Manufacturing checks; `proven`: its bounds, null where they
  differ), and a mesh's shells (an estimate): `solids`, each shell's
  volume, area and box, `watertight`, and what the mesh missed
  (`unresolved`).

**Answers as tables.** Each list of samples (eval's points, ray's rays, a
slice's loops) is a table, `{"columns": [...], "rows": [[...], ...]}`.
Objects open into columns (`d_dparams.r`); arrays stay whole (`grad`). A
ray that misses has null `nearest` and empty `hits`. `--keys
value,grad,d_dparams.wall-t` keeps only those keys (a dot steps into an
object); some appear only when named (a slice's `points`, a mass's
`inertia`). `--compact` prints a row a line.

**Threads.** `RAYON_NUM_THREADS` caps the threads (default one per core);
results do not depend on it.

**Options.** An option a command does not read is an error naming the
commands that read it, else suggesting a spelling; so is an extra
argument. Every JSON command takes `--compact` and `--timeout`.
`render --out` writes PNG. `ink spec` and `ink std` take no
model file.

## Meshes

`ink mesh` (STL, 3MF or GLB), `ink render` and `ink components` work from
a triangle mesh of each target. A mesh approximates the field, and its
report says how far to trust it.

- **Sampling.** The field is sampled on a grid of cells `--res` wide and
  meshed by manifold dual contouring: watertight, sharp edges and corners
  kept, a thin wall's two sides apart. A feature thinner than res can go
  missing. Every cell is shown by bounds to hold nothing the mesh does not
  (one sign throughout, a 0 where pieces meet flush read as the material
  round it, or the surface one sheet across it; bounds split where a value
  both pieces read crosses a level one compares it with), or, where the
  surface does not cross it, searched finer as far as a budget allows;
  each report counts every cell holding material or a gap it leaves out
  (`material`; `thinnest` across), every one undecided
  (`undecided_at_least`: nothing ruled out down to samples `searched`
  apart; of them, at a select switching between two fields, `seams_at_least`;
  at least, as a crossed cell is left `untried` where its brick of 8^3
  cells already reports `thinner_than` 2 res and no search reads it),
  and of them those
  beside a cell the surface crosses (`beside`) and the undecided ones not
  (`away`), with one `at`, and `regions`: boxes covering them (every cell
  in one; touching bricks of 8^3 cells one region, more than 16 merged
  into larger covers), each with its counts, the coarsest spacing
  searched, the thinnest found and a point in it; and warns. Where
  the surface crosses a cell, one sheet over the cell grown by half a cell
  (where its vertex may be placed) counts as drawn: the surface there lies
  within that box's diagonal, 2 sqrt 3 res (`drawn_within`), of the
  cell's vertex, and so does what else lies in that box; elsewhere only a
  sign kept throughout does. A select is continuous where its condition
  is x - y and its sides one expression of x and of y; one switching
  between two other fields (a sweep's mitres) bounds cannot show
  continuous. Each region
  and the report say how thin what may hide is (`thinner_than`: it holds
  no ball that wide): in a searched cell, sqrt 3 times the spacing
  searched; across a crossed cell's box, a seam's reach (its sides' gap
  over their least rate along the direction checked), else the box, 2
  res; round a sample on 0 drawn as air, 2 res. A mesh whose grid
  would take more bricks of cells' samples than its budget pays for is
  made at twice the res until it does not (`refused`: the res asked).
  `--max-res R`: no coarser than R.
- **Drawn whole.** A mesh made to an accuracy (not `--accuracy 0`) is
  not sampled where bounds show the surface flat: in slabs of blocks of
  cells where the field's rate along an axis is never 0 and its faces
  across that axis have opposite signs (one sheet, a height field,
  nothing else there), its slopes spread little enough over the block,
  and paint keeps one value. Slabs beside each other are joined into
  regions, meshed only in the cells round their rims and triangulated
  inside. A face of a wall under two cells thick is drawn so where
  bounds show the rest of its cell one sign; a brick holding an edge is
  sampled. Their triangles lie within `drawn_whole.within` of the
  surface (from the slopes' bounds and the rim's vertices), no more than
  half the accuracy; `drawn_within` is the larger of it and 2 sqrt 3
  res.
- **A mesh for a file** (STL, 3MF): res is the target's size / 256, no
  finer than all the targets' / 2048 (mesh a small part beside a large one
  alone, or give `--res`). It is simplified to `--accuracy` (default
  res / 30; `0`: uniform), topology kept and every vertex of the uniform
  mesh within it. `components`: size / 160, made again finer where it
  leaves material out.
- **A mesh for seeing** (GLB, render): each body's res gives its surface
  as many cells as a cube 128 across has (a render's, its box: its size /
  128), coarsened so the whole model fits a triangle budget, and made
  again finer under the
  thinnest material it leaves out, most material drawn per triangle first,
  while the budget lasts (a render's no finer than half a pixel); `--res`
  replaces it. It
  is simplified to `--accuracy` (default res / 4); shading turns at the
  shape's creases only.
- **The report** (each object of `ink mesh`): `res`;
  `uniform_triangles` (before simplifying; a region drawn whole as
  triangulated) and `triangles` (in the file);
  `watertight`; `unresolved`; and `deviation`:
  - `max`: the largest distance from the mesh to the surface, sampled at
    points at most res apart, and where (`at`). A sample, not a bound. It
    includes the sampling's own error, largest where edges, steps or walls
    lie closer than res.
  - `simplified_within`: what simplifying keeps to, bounding only what it
    adds; `simplification`, the farthest it moved a vertex.
  - `beyond`: samples with no sign change within 2 res (`beyond_at`). The
    mesh strays there, or two solids meet face to face.

  The STL header (`ink NAME: res R dev D (sampled)`) and the 3MF's
  metadata say the same.
- **GLB** (`--out FILE.glb`): the assembly, moving, for a person to view
  in any glTF viewer (three.js, Blender).
  - One mesh per body in its own frame; one node per instance, named by
    its path (`extras.body` its body), placed by its placement's numbers.
    In mm with Z up, under a root node `ink` that scales to metres and
    turns to glTF's Y up.
  - An instance whose placement reads t is animated over the file's
    `(t-span)` (else 0 to 1). A body whose shape reads t is meshed at
    `--t`, and the answer says so.
  - `extras.ink` holds the file, the span and each body's attributes.
- **`ink render`** draws orthographic views, `--size` pixels across
  (default 1024, clamped to 64..4096).
  - Default: iso, top, front and right in a 2x2 grid. `--view NAME` draws
    one (iso, top, front, right, back, left, bottom; front looks from -y);
    `--view '[x,y,z]'` looks from that direction.
  - Views are labelled (`TOP: FROM +Z`). Lines mark silhouettes, creases
    over 30 degrees and steps in depth.
  - It draws meshes for seeing, each body's res by its box (a GLB's by
    its surface). `legend` gives each target's res. Where a mesh may leave
    material out or close a gap (thinner than its res, or undecided), the
    picture draws a magenta box round each such region, with a caption
    under each view, and `left_out` lists each target and what it may
    leave out (`ink mesh`'s GLB carries the same list in its scene's
    `extras.ink.left_out`); `warnings` say it too.
    Naming no target, it names the assemblies it left out
    (`other_assemblies`).
  - `transmission` is drawn as surface transparency, without refraction.
  - `--cut AXIS[=AT]` cuts every target by the plane AXIS = AT (default
    the middle) and removes the half toward the camera (`cut` says which).
    Cut faces show in the part's colour, not hatched.

## Manufacturing checks

A part's makeability is a few claims on its field, the same for any
process, each taking the process's limit as its number (published ones
below). std gives the regions they read (`ink std --section
manufacturing`, `printing`). e is the claim's accuracy: how far off its
answer may be.

| Check | Claim | Answer |
|---|---|---|
| walls, webs, fins | `(min-thickness part > t)` | estimate: rays from surface samples |
| gaps, slots, holes | `(min-thickness (air-around part m) > g)` | estimate, the same rays across the air |
| undercuts from a set of directions | `(empty (hidden part stock (clear-air part stock lo) views))` | guaranteed for undercuts deeper than lo plus twice the accuracy |
| overhangs, printed along +z | `(empty (overhangs part z0 a h))` | a pass proves continuous support in a tested direction (below) |
| bridges at most l long | `(empty (unbridged part z0 a h l w))` | supports must touch the roof. Along x or y, or not settled; `unbridged-along` takes a direction and can fail. A failure may be up to 2w under l (below) |
| enclosed voids | `(cavities part < 1)` | proven to the accuracy, or unknown (below) |
| a flat outline to cut | `ink slice --out FILE.dxf` | its loops (Geometric queries) |

examples/machining-access.ink, mould-undercut.ink, print-checks.ink and
sls-hollow.ink show each holding and failing.

- **Walls.** `(min-thickness f > t)`, a shape metric of the design
  package: rays along the inward normal from surface samples; `value` is
  the thinnest, `at` where. Corners take their faces' thickness, so sharp
  edges are not thin. A thin spot between samples is missed: `(res r)`
  sets how fine a feature it sees, and walls no sample falls on make it
  `unknown`. A field very near zero inside material can start a ray there
  and report too small a thickness, including for a comparison labelled
  as decided without a grid check.
- **Gaps.** `(air-around f margin)` is the air out to margin past f's
  box. Its thickness is f's gaps, slots and holes narrower than margin,
  so give a margin above the gap claimed. A through hole's length is no
  gap. Between separate parts, `(clearance a b d)` is exact.
- **Line of sight.** `(hidden f stock band views)` is the points of band
  with material of f ahead of them along every view (a direction `(vec
  dx dy dz)`), looking past the stock's far side.
  - The band is `(clear-air f stock lo)`. The point and six axial probes
    at distance lo must all be air. A thin wall between probes can put
    points closer than lo in the band.
  - `(empty (hidden ...))` holds where every point can be reached from a
    setup, or pulled along the mould's pull (views +z and -z; a side
    action is one more view). A failure's witness is a hidden point.
  - Only signs are read, so the answer does not depend on how the field
    is written.
  - A line of sight has no width: cutter fit is not checked.
- **Overhangs**, built up along +z from the bed at z0. For each point
  inside the part, `(overhangs f z0 a h)` checks straight down and eight
  directions tilted a degrees from vertical, spaced 45 degrees apart.
  Each segment descends h, stopping on the bed if it reaches it first.
  A point is unsupported if every segment meets air above the bed.
  A probe landing on a floor beyond a gap gives no support. A volume
  of zeros is air. A shared face with material on both sides can carry a path.
  Use h > 0 and 0 < a < 90.
  - `(empty (overhangs part 0 45 1))` requires a proof that every
    material point has at least one solid segment. The same applies
    when clipping with `max` or comparing the mask with another shape
    using `disjoint` or `subset`. The requested
    accuracy controls refinement. Even a roof thinner than twice the
    accuracy needs support. Both print checks follow curves by bounding
    how their fields change along a probe. A path can pass through the
    overlap between two union members.
  - To draw a flush joint, use the same face on both members. For example,
    end the lower member with `(slab-z 0 h)` and begin the upper with
    `(slab-z h top)`, where their footprints overlap. For a tilted face,
    define `(def face (- (+ z (* 0.125 x)) h))` and clip the two members
    with `(max a face)` and `(max b (neg face))`. The checker can drop
    the two opposing face constraints where all the other constraints of
    both members hold. It recognizes an axis coordinate, or its sum or
    difference with an expression independent of that axis, on a continuous
    field it can analyze. Each point on the shared face there has material
    on both sides. Touching along an edge alone supplies no such patch.
    Face positions must agree exactly. For instance, a box 0.8 mm tall
    centred at 5.4 mm starts just above 5 mm in the represented numbers.
    Using a common height directly avoids that gap.
  - For a rectangular member the vertical probe can also certify the
    whole open prism from the bed to its bottom inside another member.
    A hole under the footprint prevents that certificate.
  - A failure gives a point in material and proves air above the bed
    across the entire cone of allowed directions for a positive thickness.
    At a chosen depth d up to h the cone's radius is d tan a. The checker
    encloses that disk from depth 31d/32 to d. A direct check can report
    the heights of this layer as `air_layer_z`. Field operations use a
    square enclosing the same disk. A single plane of zeros cannot
    establish a failure.
  - Slants between directions, small projecting lips, unmatched faces and
    unresolved curved joins can be `unknown`. At 45 degrees the nearest
    tested direction can lose about 2.3 degrees on a flat face. An
    undecided cell or an unproved path marks where to look. A failure
    requires a confirmed witness.
- **Bridges.** `(unbridged f z0 a h l w)` is where an overhang's ceiling
  is not a bridge at most l long, read in the ceiling. A ceiling whose
  row along x or along y reaches support both ways within l holds. One
  whose rows along both reach air, or run longer than l, first is not
  settled, never failed: a bridge may run another way.
  `(unbridged-along f z0 a h l w dx dy)` reads rows along (dx, dy) alone,
  the way the part's bridges run, and fails where one does not bridge (a
  failure's witness is there). A support must be solid for the
  full 0.3 mm below the row. Air at any depth in that column breaks
  contact. Where a row meets air is a param the claim searches.
  - A failure needs evidence of air under the ceiling and a broken or
    overlong row along the direction asked. Its angle probes read h and 0.3 mm
    below the row. Each height uses eight horizontal directions, 45°
    apart. Halfway between directions, a planar face only has the reach
    for atan(tan(a) × cos(22.5°)), or 42.7° at a = 45°. A rotated
    45° chamfer can therefore be undecided. For a pass, an angled probe
    counts as support only if its whole segment stays in material. Near the bed the segment
    ends halfway down to the bed. A floor across a slot or a wall across
    a side gap cannot supply that proof.
  - The pass must cover every possible unsupported row, including layers
    thinner than the requested accuracy. For a rectangular roof the
    checker can also prove that two bearing walls overlap its underside
    across its full width, including bearings that start at the bed under
    a shallow groove. A roof over walls at any angle to x and y (a
    slanted channel, a part turned with rotate-z) can pass too: where
    the part is upright prisms written in its own numbers, no param, the
    checker reads each roof's footprint by the rows' rule along the
    directions asked. A sloped or curved face can pass when its angled
    segment stays inside, including across overlapping solids.
    The same shared-face rule as for overhangs can certify a flush
    bearing or sloping shelf. Complicated joins can remain `unknown`.
  - A hole's diameter alone does not determine its bridge span. A
    teardrop can have support along both flanks. In a round hole, a
    horizontal row just below the crown can enter the opening before
    its far bearing. Even a 3 mm hole can fail a 10 mm bridge check.
    A row above the opening can bridge if it has bearings on both sides.
  - Support is local. A hanging rib can supply the 0.3 mm column under a
    row. Its own underside is checked in turn. A pass does not certify a
    continuous supporting path all the way to the build plate. Use the
    bridge region as the argument to `empty`.
  - Rows are read for supports at steps of
    w, so a failure may be a span over (floor(l / w) - 1) w (at l 10, w
    0.8: 8.8), or one resting on a wall thinner than w (FDM: w 0.8, the
    thinnest wall).
  - A ceiling over walls at an angle to x and y is not settled, along
    any direction: its support columns stay undecided beside them.
- **Enclosed voids.** `(cavities f < n)` counts f's enclosed voids, air
  no path leads out of (trapped powder or resin), with air deeper than
  the accuracy (a length, as for shapes: Claims). A shape metric of the
  shape package, counted by interval bounds on cells halved where they
  do not decide, down to its box / 20000, or `(res r)`, within a budget of
  1 M cells and the call's time: a void sealed by a wall at least a cell
  thick is found, however steep its field. A proven count is exact.
  - Air no deeper than the accuracy (a slit thinner than twice it, faces
    placed 100 nm apart or meeting flush inside the material) is set
    aside, as `empty` sets aside a shape so thin: the claim holds `to the
    accuracy`, `set_aside` saying where. A smaller `(accuracy x)` counts
    thinner voids.
  - Where they leave it undecided, the count is bounds: the claim is
    `unknown` where they straddle its bound (`undecided-cells`: where, the
    undecided cells' sizes, and what stopped it). At the finest cells: a
    wall or gap thinner than them. A flush mouth (0 across it, air either
    side) is air: cells across it split where both pieces read one value
    (a coordinate, a turned radius) at a level one compares it with. At
    the budget (or forecast past it)
    or the time limit, a finer res does not help: count a smaller region.
  - An escape hole's width is a gap.
- **Draft** is a `range` on a slope (Claims: The other claims).

Published minimums, in mm (read 2026-10-31):

| Process | Walls | Holes | Parts that move | Overhangs, bridges | Escape holes | Smallest feature |
|---|---|---|---|---|---|---|
| FDM | 0.8 | 2 | 0.5 | 45 degrees; 10 | | 2 |
| SLA | 0.5 supported, 1 unsupported | 0.5 | 0.5 | supports always | 4 | 0.2 |
| SLS, MJF | 0.7 | 1.5 | 0.3 (0.1 for a connection) | none needed | 5 | 0.8 |
| Material jetting | 1 | 0.5 | 0.2 | supports always | | 0.5 |
| Binder jetting | 2 supported, 3 unsupported | 1.5 | | | 5 | 2 |
| DMLS | 0.4 supported, 0.5 unsupported | 1.5 | | supports always; 2 | 5 | 0.6 |
| CNC milling | 0.8 metal, 1.5 plastic (0.5, 1.0 feasible) | 2.5 for standard tools (drilled: over 1) | | | | |
| Laser cutting | | half the sheet's thickness, at least 0.38 | | | | |

Sources:

- Protolabs, Design Rules for 3D Printing
  (https://www.protolabs.com/media/viknghfd/pl_3dp_design_rules_en.pdf).
- Protolabs Network (Hubs), How to design parts for CNC machining
  (https://www.hubs.com/knowledge-base/how-design-parts-cnc-machining/):
  also pockets at most 4 times as deep as wide, inside corners rounded to
  a third of the depth.
- SendCutSend, minimum geometry
  (https://sendcutsend.com/faq/what-are-best-practices-for-minimum-geometry/).

## Motion

Time `t` is an input like x, y and z, in seconds. Placements and fields
may read it: `(instance arm arm (at) (turn z (* 360deg/s t)))` turns once a
second. A joint's rate is in deg/s (turn) or mm/s (slide); `3000rpm` is
18000 deg/s, and a speed in rpm is a multiply, `(* r 1rpm)`. `--t T`
evaluates any command at one moment.

`(t-span 0 1s)` in a file makes `ink require` and `ink interfere` check
every claim and pair that reads t at every moment of the span;
`--t-span 0:1` or `--t T` replaces every file's span. A claim takes its
own import's span. In `ink check`, `import_t_spans` lists them by file;
a single copy with a span gives a pair of numbers, e.g. `[0, 1]`. Two or
more copies with spans give an object keyed by namespace, e.g.
`{"short":[0,1],"long":[0,10]}` (the empty key names a plain import). An
`interfere` pair, the hull of the model's spans, whatever the target. A
claim reading t under no span is asked at one moment (its `one_moment`,
named in the answer's `note`; `ink check` warns, with `--load` of claims
only). One that does not read t, or whose fields do not change over the
span, is checked once: at `--t` where in the span, else at its start
(`ink check`: `t_span`).

A claim that reads t is asked with t ranged over the span, as a
tolerance band is (Tolerances; `param_spans`: `t`): a hold the search
decides covers every moment, a failure gives a moment that violates it.
What it leaves undecided in the call's time, or cannot range, is
estimated (`"accuracy": null`) at sampled moments (`steps`, `step`,
`basis`: 17 moments evenly spaced, the span's ends included, as a band's
17 values; or a package's own steps; sampling stops at a failure or the
time limit), and the `hint`
says why (`motion.not_ranged`, `motion.t_ranged`). A failure at a sampled
moment is guaranteed only when its `accuracy` is a number; a moment where
what the claim reads has no value (a solve that does not solve) makes it
`unknown`.

A hold at every sampled moment is never a yes by itself. Each step
between two of the 17 moments is asked with t ranged over that step
alone, where a part travels little and the question is whether the room
the claim has at the moments is more than what it reads can move within
the step; a step left undecided is halved and asked again. The work is
bounded, the same on every machine: 64 asks in all, each search given
8192 boxes, and no more asked once the steps left cannot all be halved
within the 64 or none of the first 4 is shown. No step is asked where
the search over the whole span itself ended at its limit of boxes (a
step's, given fewer, would end sooner), and none in which a body travels
over 16 times its own size. Every step held: the
estimate stands (`motion.not_ranged` says so).
A witness in a step: a failure, as any t ranged gives. A step left: the
claim is `unknown` (not settled), `motion.loose` the steps and the `hint`
where and why, and `ink require` exits 3. `ink interfere` asks its moving
pairs the same once every pair has an answer. A claim read at a package's
own steps is the package's to cover between them (it marks the steps it
does not follow, and they are not settled); what such a claim reads of t
besides is shown over each step by interval bounds, or that step is not
settled.
A motion that repeats is asked over one period of it. Each shape that
turns about its own axis repeats by its turn symmetry (a gear every
tooth, whether its teeth are a fold of its angle or the sine and cosine
of its tooth count times it); a wheel that repeats but for a part of it
(a bore's flat, a set screw's hole) is asked as the part that repeats and
as what is added to it, each alone, the wheel lying within the two;
shapes carried round one axis together (planets on their carrier)
repeat so in the frame that goes round with them, where each only turns
in place; a repeat along a path (a belt's teeth, a floor of its shift
over its pitch) repeats every pitch, and where its period and a wheel's
are a rounding apart (a pitch radius is a rounded number) the difference
over the motion is charged to the accuracy. A period whose search ends at
its limit of boxes is asked again in 16 pieces, each a search of its own
share. `ink interfere` asks a pair whose own field moves (a belt) alone
so, once every pair has an answer.

A claim about shapes read at a package's own steps (a part a simulation
moves) is asked with t ranged over each step, the part moving about that
moment as the package says it does: a witness there is a failure, and a
step not shown is not settled.

An overlap is shown over a range (of t, of a toleranced param) by a
point: one found in both shapes that is in both at every value of the
range, by interval bounds of the fields there, in the frame the first
shape rests in. A point in both at one value alone says nothing of the
rest.
Not shown over any step, so not settled over a motion where its fields
change: a contact (`touches`: claim what can hold at every moment, a
clearance or an overlap); a metric, integral or mean of a shape that
changes with t (it deforms, or its parts move apart), or of a body but
as below; a value a package sets on t that chooses no steps.

A body that only moves rigidly is asked in its own frame: its field is a
shape that reads no t, put where it is by placements and warps (`move`,
`rotate-z`, an instance's `turn`) each a rotation exactly with a finite
translation at every moment of the span (a turn through a fixed angle of
floats is no rotation exactly, nor a mirror or a stretch). Its `motion`
says `how`: `rigid`.
- A thickness (`min-thickness`) or a count of voids (`cavities`), the
  same however the body is turned or moved, is asked once, at a moment of
  the span; its answer, an estimate or a proof, is every moment's.
- An integral or a mean whose weight reads no t (and, of a body that
  moves, no x, y or z) and whose integrand is affine in x, y and z (`1`,
  a density, `x`, a direction's dot with the point; its coefficients may
  read t; read as real numbers, each finite over the span and each
  divisor apart from 0) is its own frame's integrals measured once, carried by its
  placement: the mean of x over `(move (box 10 10 10) (* 10 (/ t 1s)) 0 0)`
  is 10 t. Bounded with t ranged over pieces of the span, each halved
  until it decides (4096 in all): proven, or failing at the middle of a
  piece every moment of which fails (its `t`). Undecided, it is asked as
  above. `(mean b (* x x) ...)` is not carried: not settled over a motion.

- **Name the moving part first.** Two shapes (`disjoint`, `clearance`,
  `subset`, `overlaps`, `touches`) are searched in the frame the first
  rests in, so gears or a carriage on its rail are one search and a pin
  turning in its hole does not move: a pin's `subset` of its hole, a
  link's `disjoint` of the sprocket. Where in the second's frame, or the
  world's, nothing of the claim moves, there instead (a shaft turned
  about its own axis, any axis, however placed, in a still bore: asked
  once). Two shapes riding one moving
  instance rest together in its frame and are asked there: neither moves,
  and what moves the instance is not read. Witnesses are given in the world. A part
  turning about a point the motion holds still, though placed through
  another angle (a coupler on its rocker's pin, met by the loop's
  closure), is searched turning about that pivot, its bounded drift
  charged to the claim (a clearance d asked at d plus it; an emptiness at
  its accuracy less it). A part that does not turn, standing a fixed way
  from the frame's origin (a belt seen from the wheel at its far end), is
  the frame's turn seen backwards about that one point, exactly, with no
  drift. What rides with the frame cancels: a sum in which the motion
  drops out in real numbers (a belt's teeth read at a jaw's point plus
  its travel, less the same travel) is read as the sum left, so the run a
  jaw clamps does not move in the jaw's frame. A turn rounded off in
  floats (a CAD converter's) may leave a pair measured.
- **Speeds** (`speed-limits`, `diff` along t) require a bounded drive
  rate. A claim reading a speed is `unknown` where the drive may jump,
  including at a single requested moment or either end of a span. A
  further derivative or arithmetic cancellation still reads that speed.
  A warp checks it at the original drive's local coordinates and clock.
  Any speed used as a warp coordinate must itself be defined. A probe
  checks the speed wherever the probe can be. A simulated coordinate and
  its speed are continuous over completed integration steps. Their claims
  are estimates and still check the commands used to reach that state. A turn's speed limit
  reads the turn, so its angle wrapped by `mod` at a whole turn (`(mod (* 360deg/s t) 360)`), or switched by a
  `select` between branches a whole turn apart, is decided. Write
  a trapezoidal position profile with `(trapezoid-position time acceleration top-speed distance)`
  from std. Its phases meet exactly in real arithmetic.
- **Repeating motions** are searched over one period: a claim reading t
  only through a drive (t, or one value all its reads of t pass through:
  a spin-up's angle) and, through it, only through the sine and cosine
  of turns at constant rates of it, over a span where the drive travels
  past their common period, taken from their exact values (`(* 3600deg/s t)`
  and `(/ (* 2 drive) 3)` share a turn; floats such as 0.3 and 0.7 share
  only a period of about 1e16 s). A shape turned about its own z that is
  the same turned by 1/n of a turn (a gear's teeth) repeats whenever the
  drive turns it so: a gear pair at mesh-turn by one pinion tooth, its
  rounding (about 1e-13 mm) charged to the accuracy. Through a drive
  other than t only a hold is an answer; else t is ranged.

The answer's `motion` block gives `how` (`t ranged`, `rigid` or `measured`), the `undecided` t ranges, and for `measured` the moments
checked (`steps`, and `step` or a package's `basis`). A failure gives its
moment, `t`.

```lisp
(import "std")
(t-span 0 1s)
(body crank (shape (move (box 40 6 4) 20 0 0)))
(body post (shape (move (cylinder-z 3 20) 35 35 0)))
(assembly root
  (instance post post)
  (instance c crank (at) (turn z (* 90deg/s t)) (limits -10 100) (speed-limits 0 100deg/s)))
(require clear (clearance (shape-of root/c) (shape-of root/post) 1))
```

`clear` and the joint's `limits:root/c` and `speed-limits:root/c` are
checked over the whole turn.

## Data for other packages

Packages add top-level forms, values, metrics, commands and clauses. The
`shape` package owns the metrics of Manufacturing checks. Others (simulation, dynamics) come in distributions
built with them: `ink spec` lists theirs after this reference, `ink help`
their commands. They read a body's attrs (`(attr spin-inertia J)`) and
clauses of their own: in an assembly, on an instance with a motion after
its limits (an actuator's `(effort-limits LO HI)`), or in a relation (a
belt's `(couple ...)`).

A form, clause or command no package of the running `ink` defines is an
error naming it. In a relation or an instance, `(optional CLAUSE)` is dropped
where no package owns CLAUSE; std's gear meshes give their coupling so,
and their claims stand without it, and the catalog's actuators their
`effort-limits` and `stiffness`, so the kernel alone places their joints. Only a
package's clause may be optional: a kernel clause in it, a `for`'s
included, is an error. Where a loaded package reads that kind of clause
(dynamics: joints' and relations'), an optional clause none reads warns,
naming it (a misspelling would drop its check). A dropped optional joint
clause ending in `-limits` appears in `not_checked` for each placement
under the run file's root, including in a kernel-only build.

## Tolerances

A tolerance, how far a real part may differ from its drawing, is declared
on the param as deviations from its value, as a drawing writes it:

```lisp
(param bore-d 8 (tol 0 0.022))            ; 8 H8: 8.000..8.022
(param shaft-d 7.98 (tol -0.008 0.007))   ; 8 f7: 7.972..7.987
```

The sizes it allows are the param's band. Give a value inside it, which
commands that work at one value measure at: a value outside its band
(`(param shaft-d 8 (tol -0.028 -0.013))`) is warned of by `ink check` and
flagged `outside_band` by the others.

- `ink require` and `ink interfere` prove every claim for every value of
  each ranged param (one with a band) it reads; other commands work at
  the values. `--param bore-d=8.1` moves the band with the value;
  `--param k=lo:hi` gives an absolute range; `--param k=v:v` pins k;
  `--nominal` drops the declared bands. `ink check` lists each param's
  `tol` and `band`, and of one an import sets, its `written` default and
  `set_at`.
- `(ranged lo hi [at])` is a range of the question, not the part (how far
  a tool is withdrawn, every angle from 0 to 90). It is ranged as a band
  is, but `--nominal` keeps it, and `require` and `interfere` refuse
  `--param` on it, so a claim over it is never decided at one value.
  `eval`, `render`, `mesh` and `distance` take it at `at` (default lo),
  or at `--param NAME=V` within the range. It is named by its def, else
  by its form and a number (`--param "a 0=0.25"`).
- Each claim lists `param_spans`, the ranged params it reads and their
  bands. A failure names `params`, each one's value at the witness.
- By claim kind:
  - shape claims, ranges and values hold for every value
    (`(value (square p) >= 0.25)` over p in [-1, 1] fails at p = 0);
  - `overlaps` and `touches` fail where apart at every value (`touches`
    also where too deep at one), else are `unknown`: a contact found at
    one value says nothing of the others, so claim a fit by its sizes,
    `(value (- pin-d hole-d) > 0)`;
  - volumes, integrals and means are estimated at evenly spaced values;
  - shape and package metrics whose form reads one are `unknown`: check
    them at one size (`--param name=v:v`, or `--nominal`).
- Refused, with an error naming why, never a proof at one value labelled
  as covering the band:
  - a ranged param used at load time, while the file is read (an `if`, a
    count, a sampled field, a package form's settings), or beyond the band
    a `(bounds f)` was found over;
  - a direction not proven nonzero over the bands and the motion;
  - a bound that reads one, except in value and mean claims.
- Ranges are boxes: a round positional zone written as dx, dy in [-r, r]
  is a square. Write round zones in polar params.
- Statistical tolerancing (root-sum-square, Monte Carlo) is not covered.

## Standard library

std is written in ink (std/std.ink, std/materials.ink).

- `ink std` lists its sections and definitions' signatures, then the
  forms built into the language.
- `ink std --section NAME` prints, with comments, the section whose
  banner (heading comment) is NAME, or the definitions it names; a
  built-in form's name (`polygon`) gives this reference's entries for it
  and the sections writing it, the one writing it most whole; else what
  holds the name (`gear`: every gear builder). `--all` prints all of std.
- `ink spec --section NAME` finds this reference's sections the same way,
  by title or path (`Claims/Results`), else looks
  in std and the text of both.
- A call of an unknown function suggests std's name for the CAD word
  (`hole`: subtract a `cylinder-z`; `union`: `min`;
  `extrude`: `extrude-z`) and near spellings.

Sections: scalar helpers (`clamp`, `mix`, `blend`, `atan`,
`acos`, `asin`); combination (`subtract`, smooth blends, fillets,
chamfers, shells); appearance (`paint`, `color-map`); warps (`move`,
`rotate-z`, `scale`, `twist-z`, `extrude-z`, mirrors,
repeats and patterns, thread frames); primitives (2D and 3D, exact
distances); lattices; manufacturing; printing; sheet metal; paths; belts
and chains; power transmission (gears, threads, pulleys).

### Lettering

Solid letters from outline fonts, for names, signs and labels to print,
cut, emboss or engrave: `ink std --section lettering`.

- `(font NAME [(AXIS v) ...])` (built in): `"sans"` (Lexend), `"rounded"`
  (Fredoka), `"serif"` (Bitter), `"script"` (Pacifico, joined), vendored
  under open licences (std/fonts), or a .ttf or .otf file, at a point of
  its variation axes: `(font "sans" (wght 700))` is bold. A variable
  font's overlapping and self-crossing contours are joined into what it
  fills.
- `(lettering text F h)`: a layout, cap height h, centred;
  `(lettering-styled text F h tracking leading align)` sets the spacing,
  lines (`"\n"`) and alignment. The text is shaped as HarfBuzz sets it
  in that font: ligatures, a script font's joining forms, accents composed
  or placed on their letters, kerning. `(lettering-2d L)` is its letters,
  a 2D profile as `polygon` is (each glyph's lines and arcs within 4e-4 em
  of the font's curves); `(lettering-box L)` their box. For a layout std
  does not make, `(glyph-run F text)` (built in) is the shaped text: each
  line `(vec glyphs advance clusters)`, each glyph `(vec outline data)`:
  its outline a 2D profile at its pen, its data its `"x"` and `"y"` from
  the line's start, its `"cluster"`, `"advance"` and `"box"`, in ems.
- Styles compose as fields: raised or engraved (`extrude-z`, `subtract`),
  `extrude-z-chamfer` and `extrude-z-round` tops, outlined (`shell`:
  ink check warns that its corners mitre),
  grown into a backing plate (`(- g d)`), along a path (`lettering-along-2d`:
  the layout's x = 0 at a station, its y = 0 on the path), round a circle
  (`(lettering-around-2d L r a out)`: centred at a degrees, tops outward
  reading clockwise, or inward), round a cylinder (`wrap-cylinder-z`, any
  flat design).
- Claims of the 2D letters are asked in the plane (Where analysis
  looks): that they fit (`(subset g face)`), their strokes
  (`(min-thickness g > 0.8)`, an estimate), their gaps (`clearance`). For
  every text a text param allows, `(lettering-max chars n F h tracking)`
  bounds the box of any line of up to n of its characters, soundly, over
  every glyph the font has (whatever it substitutes), their advances,
  outlines and all the font's kerning: loose, never short. An error, saying
  why, where it is not counted simply (a combining mark, a script with its
  own shaping rules, mark, cursive or contextual positioning, a table or
  format it does not read): ask such texts as cases. A specimen of them all,
  `(lettering-styled chars F h h 1 0.5)`, is every letter at once.
  examples/name-ornament.ink.
- Lines run left to right: no right-to-left or vertical text.

### Curves

- `polyline-distance` (3D and `-2d`): exact.
- `(polygon pts)` (built in): a closed profile of lines and arcs, any
  simple one, either winding. Each point is:
  - `(vec u v)`, a corner;
  - `(vec u v b)`: the edge to the next point is an arc of bulge
    b = tan(sweep / 4), positive counter-clockwise (DXF's convention);
  - `(vec u v b r)`: the corner rounded to radius r.

  Its numbers are read at load time (a param in them is baked, not
  ranged). Its value is 1-Lipschitz, never over the distance. The profile
  must be simple: edges that meet other than at a shared corner are
  refused, saying where. It is exactly its corners, arcs and rounds (a
  circle's centre and radius known within bounds of their rounding, near
  1e-14 of the circle, where a claim is unknown): a corner is dropped only
  exactly on the line between its neighbours, and edges that touch are
  refused. Except: an arc of half a turn or more, one whose triangle of
  its chord and end tangents meets another edge, one whose circle would
  need a number outside ink's range, or a round whose centre or tangent
  points ink cannot bound (inside its edges), makes the whole
  profile read by its allowance: a corner or an arc's bulge within 1e-9
  of its points' bounding-box diagonal of the line that would replace it
  is dropped, and that line used (drops in a row add up); beside a shared
  corner, a meeting whose sliver between the two edges is within 100
  times 1e-9 of the diagonal (a tangent joint's rounded data) is allowed.
  The shape so read carries that allowance as a band: a claim on it, or
  a value read at a point, holds or fails proven where it does so over
  all of the band, and is not settled where the band could change it,
  saying it reads a shape built within an allowance. An integral, a
  mean, a metric, or what reads the shape through a solve (a motion's
  efforts) is an estimate and does not carry it.
- `(polygon pts drafts z0 z1)`: a moulded or drafted wall. pts is the
  profile at z = 0; each edge's wall leans outward by its draft (degrees;
  one for all, or one per edge) as z rises; rounds keep their radius;
  straight beyond z0 and z1.
  - Straight walls are exact: each wall the plane through its drawn
    edge, each corner where two meet, a round between two moving with its
    corner. Arcs too: each arc's wall a cone (its circle growing or
    shrinking with z), each corner beside one where its walls meet, its
    segment read as undrafted (at its chord, or its cone within a kite
    moving with it). What it cannot build so (a round beside an arc, a
    flat arc, a kite it cannot place, and as undrafted an arc of half a
    turn or more or two corners) is read by its allowance, as above.
    Drafts all 0 are the profile undrafted.
  - Load time proves the construction holds from z0 to z1, else refuses
    the profile, naming the corner: a round no longer fits or an edge
    would vanish (where CAD would trim the walls), or "ink cannot
    certify" (coordinates far larger than the profile, walls meeting
    nearly tangentially).
- `(polygon-distance pts)`: corners only (a bulge or rounding is
  refused), no `lipschitz` promise, exact outside, so an outward offset
  rounds convex corners.
- `bezier3-points`, `spline-points` (Catmull-Rom), `helix-points`,
  `regular-polygon-points`.

### Turned profiles: `revolve-z`

`(revolve-z g)` (built in) turns a 2D profile g(x, y) about the z axis:
g read with x the radius, sqrt(x² + y²), and y the height, z. Draw the
profile from x = 0; the axis is inside where the profile touches it.

- An edge on the axis is no surface. Read at the radius (never below 0)
  it is a term at most 0, and 0 only on the axis: -k x + c with k > 0 and
  c <= 0 (its coefficients proven, not rounded), or k x + c, c >= 0,
  negated. Such terms are left out of the profile's maxes (and their
  negations out of its mins), read through names, `lipschitz`, `neg`,
  `min`, `max`, multiplying or dividing by a positive constant, and
  `select`. So `polygon`, a moved or scaled one, `polygon-distance` and
  `subtract`s of them read inside along the axis. Every sign off the axis
  is the profile's; only values inside, near an edge left out, change
  (they read the remaining walls). Values are the profile's kind
  (`polygon`: never over the distance), not distances to the turned
  surface.
- A profile reaching the axis some other way is not inside along the
  axis: an edge on x = 0 written with `abs` (a `rect-2d` moved to touch
  it) reads 0 there, so a probe is not inside and a ray along it meets a
  surface; one offset inward, `(+ g r)`, reads r, a hole of radius r.
  Write it as a `polygon` from x = 0, or carry the profile past the axis,
  to x = -1.

### Sheet metal

A sheet-metal part is designed flat and folded from it, so its folded
body and flat pattern agree.

- `(sheet outline holes th k bends)`: the flat (outline, holes,
  thickness, K-factor, bends).
- `(sheet-bend a b angle r)`: the material left of a -> b folds by angle
  toward +z at inside radius r. At an angle to an edge, a and b lie past
  it, so the bend spans the material.
- `(sheet-folded S)`: the folded part. `(sheet-flat S)`: the flat pattern
  (`ink slice --out FILE.dxf` is its cut file). `(sheet-bend-lines S)`:
  each bend's centre line, a to b within the outline, a body's marks:
  `(mark bend (sheet-bend-lines S))` puts them in the cut file.
- `(sheet-frame S u v side)`: a frame at a point of the flat, carried to
  the folded part.
- Bend rules, regions a design leaves empty: `sheet-short-flanges`,
  `sheet-holes-near-bends`, `sheet-unrelieved`, `sheet-collisions`; and
  `(value (sheet-min-radius S) >= th)`. `ink std --section 'sheet metal'`
  gives each rule's source.
- Not built: crossing bends, closed corners, flattened hems (a hem is a
  180 degree bend of r > 0). A bend runs straight across its flange.

```lisp
(import "std")
(def S (sheet (polygon
 (vec (vec 0 0) (vec 60 0) (vec 60 40) (vec 0 40)))
              (vec (move (circle-2d 2.5) 48 20 0)) 2 0.4
              (vec (sheet-bend (vec 30 40) (vec 30 0) 90 2))))   ; x > 30 folds up
(body bracket (shape (sheet-folded S)) (frame upright (sheet-frame S 48 20 -1)))
(body flat (shape (sheet-flat S)) (mark bend (sheet-bend-lines S)))   ; the cut file
(require flanges (empty (sheet-short-flanges S 4)))
(require reliefs (empty (sheet-unrelieved S)))
```

See examples/sheet-metal.ink.

### Paths and power transmission

- **Paths of lines and arcs**, what a belt, chain or cable runs along:
  - pieces `(path-line a b)` and `(path-arc c u r f sweep)` (sweep in
    radians), each tangent to the next, which may read t;
  - `(path pieces up)` open, `(path-loop pieces up)` closed;
  - `(path-length P)`; `(path-at P s)`, the frame s mm along it (x along
    the path, z `up`); `(path-point P s)`;
  - `(path-link-at P s len)`, a rigid link from s to s + len: instances
    at `(path-link-at P (+ (* i pitch) shift) pitch)` with a shift that
    reads t circulate (examples/chain-drive.ink);
  - along a planar path, `(path-sweep SECTION P)` bends a 2D section;
    `(path-repeat F P PITCH SHIFT)` places F whole at every station
    k PITCH + SHIFT (a belt's teeth).
- `(attr carried (carried PROFILE P LENGTH SHIFT))` marks a body that is
  a profile carried along a path (a belt, a hose), so a viewer may bend
  one mesh along the moving path.
- **Belts and chains**: `(belt-along b P width side shift)`, a
  synchronous belt of any tooth profile along any planar pitch path. Its
  paths: two pulleys (`sync-belt-path`); wheels and back-side idlers
  (`belt-loop-path`, `sync-belt-around`); an open route (`belt-route`); a
  roller chain (`roller-chain-path`). `ink std --section belt`.
- **Gears**: `(involute-gear-2d z module pressure-angle backlash)`,
  `involute-gear-shifted-2d`, `involute-rack-2d`,
  `straight-bevel-gear`. A negative z is an internal gear's (a ring's)
  teeth, filling all beyond its roots: clip it with a circle. A gear
  body declares the numbers it is built from, `(attr gear (vec z module
  pressure-angle))`, z as a mesh writes it; a rack `(attr rack (vec
  module pressure-angle))`; a cluster (gears on one body) none.
  `(mesh-distance za zb module)` and `(mesh-turn
  za zb bearing a)` place a pair. In these, and in `gear-mesh` and
  `rack-mesh`, a gear (za, zb, z) may be `(vec teeth centre-teeth)`, one
  cut to run at another count's centre distance; `bevel-mesh` refuses
  it.
- **Gears in mesh** are relations, A the input:
  - `(gear-mesh za zb module [pa [backlash [efficiency]]])`: spur, each
    axis its gear's z, a ring's pair with the ring's negative count
    (`(gear-mesh 12 -36 1)`, mesh-distance and mesh-turn alike; two
    rings, or a ring of no more teeth than its pinion, refused); the
    efficiency is its coupling's, read by a dynamics package;
  - `(bevel-mesh za zb module face ...)`, `(rack-mesh z module ...)`.
  - pa is the pressure angle (degrees, default 20). Under 10 is an error
    (likely a backlash written in its place).
  - Its claims say the gears sit as a mesh, pitch surfaces parted by no
    more than the backlash allows: `centre-distance`, `parallel-axes`,
    `faces-meet`; a bevel's `pitch-cones`, `square-axes`, `apexes`; a
    rack's `pitch-line`. And the gear an end declares is the mesh's,
    every number exactly: `a-teeth`, `a-module`, `a-pressure-angle`,
    `a-centre-teeth` where either is shifted, B's alike (the terms give
    both). Numbers declared any other way (the other kind's, of another
    shape) are refused.
  - It has no contact zone: interfere checks the teeth over the pinion's
    angle, so a wrong ratio or phase is a clash. It reads no tooth's
    shape: a body that declares no numbers meets the mesh's only there.
- **Belt and screw relations**: `belt-drive`, `chain-drive`,
  `belt-span`, `cord`, `lead-screw`, `held`, `tensile-stiffness`: couplings read
  by a dynamics package; an ink without one refuses them (Data for other
  packages).
- **Threads**: `(thread-iso d p)`, an ISO metric external thread along
  +z. Its field is the distance to the profile, so `(- (thread-iso d p)
  clearance)` is a tapped hole's void. Pulley grooves
  (`gt-pulley-teeth-2d`, `htd-pulley-teeth-2d`) are nominal, with no
  running clearance.

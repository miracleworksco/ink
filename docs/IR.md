# ink language reference

A document is a sequence of top-level forms in S-expression syntax. `;` starts
a comment.

A code block here that begins `(import "std")` is a whole model: it loads
as written. Any other is an excerpt, read with std and the names the text
around it gives.

## Tour

A name in parentheses, like (Claims), is a section with the detail: `ink
spec --section Claims` reads it.

A design is `.ink` text files. `ink COMMAND FILE` runs a command on one
and prints JSON.

**Shapes.** A shape is a field: one function f(x, y, z) of position in mm,
negative inside; its surface is where it changes sign (0 where faces meet
flush, the same sign either side, is none; a region it reads as 0 is air).
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
can set (`--param hole=6.2`). `(def name expr)` names a value, `(fn name
(args) expr)` a function. A param's name and an import's namespace are
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
MM)`, in params and the time `t` in s. `(t-span 0 2)` is the time that
claims about the motion cover.

**Claims.** `(assert NAME CLAIM)` states a requirement: `disjoint`,
`clearance`, `touches`, `subset`, `value`, `volume`, `integral` and more
(Claims). `ink assert FILE` answers each `holds` (met, to its accuracy),
`fails` (not met) or `unknown` (undecided, often with a hint saying
what would decide it). The accuracy is how far off an answer may be: 0.05
mm by default, `(accuracy x)` after the claim to change it; `null` marks
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
(assert seated (touches (shape-of root/p) (shape-of root/b)))
(assert free-fit (clearance (shape-of root/p) (shape-of root/b shank) 0.2))
```

Both claims hold. With `--param hole=6.2`, free-fit fails: its witness is
two points 0.15 mm apart.

**The loop.**

1. Write the model.
2. `ink check`: what loaded, and errors and warnings that say how to fix
   it.
3. `ink render --out FILE.png`: four labelled views. Fix until it looks
   right.
4. Claim what the design depends on, and run `ink assert`.

A wrong placement passes check, so claim where each part must be
(Placement by example).

`ink help` lists the commands, among them `interfere` (parts that
overlap, at rest or moving) and `mesh --out FILE.stl` (or .3mf, .glb). A parts catalog,
installed apart from ink, adds parts: `(import "parts/FAMILY/NAME.ink" as
b)` names its bodies `b.NAME` and its record `b.spec` (Top-level forms: Libraries).

## Units

A bare number is in the language's unit for its quantity, in every form:
lengths mm, masses g, densities g/mm³, forces N, pressures and stresses
MPa, torques N mm, temperatures °C, times s, frequencies Hz, powers W,
speeds mm/s, turning rates deg/s, gravity m/s², conductivity W/(m K),
specific heat J/(kg K), convection W/(m² K).

A number may name its unit instead, and is converted: `(box 4in 2in
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
- Only numbers carry units. Scale an expression by multiplying: `(* n
  1in)`.
- Where a form takes a known quantity (a package form's forces and
  temperatures, a claim's bound on a metric, `clearance`, `volume`, `res`,
  a joint's limits, a body's density), another quantity's unit is an
  error saying what it takes and where that number is written. A number
  keeps its unit through names, lets, calls, the branch an `if` takes and
  a `for`'s elements, and through sums, whose terms with a unit must be
  one quantity (a bare term is in the language's unit, and is allowed); a
  product's units combine, and it is not checked. A param's tolerance
  deviations with a unit are its default's quantity, where that has one.
- A density outside what real materials have (5e-8 to 0.023 g/mm³) is an
  error that says how to write it. 0 is massless.
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
  `(deg->rad d)` converts, with pi exact where d reads x, y, z or t
  (constants and params are their floats, as everywhere): an angle
  wrapped at a whole turn (`(mod d 360)`) reads the same either side of
  the wrap, so motion across it is proven continuous (`(* d (/ pi 180))`
  rounds pi: not). A wrap counts as whole turns only where it is 360 k
  degrees exactly (1e20 is not).
- `path-arc`'s sweep, and the angles std's belt and chain tools return
  (`sync-belt-theta`, `roller-chain-theta`, `belt-route-entry-angle`,
  `belt-route-sweep`).

`ink check` warns where `sin`, `cos` or `tan` gets what reads as degrees
(`(sin 30)`, `(sin (* 360 t))`). Angles take no unit: `30deg` is an error.

## Fields

A field is a real-valued expression over the inputs `x y z t` and any
declared parameters. A **shape** is where a field is `< 0` (Shapes). It
need not be an exact distance, though distance-like fields
(|gradient| near 1) mesh best, and offsets (`(- f d)`: f grown by d,
shrunk for d < 0) are exact only for exact distance fields. The
`clearance` claim is sound for any field (Claims).

The kernel's operations are the whole vocabulary; everything else is
library code:

```
(+ a b ...)  (* a b ...)  (min a b ...)  (max a b ...)   variadic
(- a b)  (- a)  (/ a b)  (atan2 y x)  (pow a b)
(neg a) (abs a) (sqrt a) (square a) (exp a) (ln a)
(sin a) (cos a) (tan a) (floor a)
(deg->rad a) (rad->deg a)  degrees in radians and back, pi exact
(select c a b)          a where c < 0, else b
(emax a b ...)          Euclidean max: sqrt(a^2 + b^2) where both are > 0, else max
(warp e fx fy fz [ft])  e with x:=fx, y:=fy, z:=fz (and t:=ft)
(lipschitz k f)         f, promised |grad f| <= k in x y z
pi inf x y z t          constants (inf: +infinity) and inputs
```

On shapes, `min` is the union and `max` the intersection. std's `(subtract a b)` is `(max a (neg
b))`. `emax` is the exact intersection of exact fields with orthogonal
gradients: `(emax qx qy qz)` is an exact box, `(emax g dz)` an exact
extrusion.

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
a ranged param, x, y, z or t is proven as a real number. Spatial
arithmetic that cancels as written, `(- x x)`, is the constant 0, and
what is built on it folds in floats.

### Derivatives

`grad`, `diff`, `ink eval`'s gradient and the optimizer take each op's
derivative on the branch its value took, by the chain rule, and agree bit
for bit. At an exact tie (`min` or `max` of equal values, `abs` at 0) a
derivative is the average of both sides. `select`'s is its branch's,
`floor`'s 0. A solved value's derivative is the implicit function's
(Solving constraints). Proofs bound a derivative over a region by its
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
                               (evaluated only then) where f is empty, else that is an error; refused
                               for a field that reads t (found at load, it would hold at t = 0 only)
(grid-bake f [(res mm)])       sample a field onto a grid (trilinear; a cache for expensive fields)
(import-stl "part.stl" [(res mm)] [(region (x0 y0 z0) (x1 y1 z1))])
                               a triangle mesh as a signed distance field (sampled; with a
                               region: only the part of the mesh inside it, sampled finer)
(redistance f [(res mm)])      signed distance to f's surface, sampled on a grid: makes offsets,
                               shells and fillets work on fields that are not distances (approximate)
(diff f v)                     derivative with respect to x, y, z, t or a param, or a def that is one:
                               symbolic, a field like any other; std's (grad f) is (vec dfdx dfdy dfdz)
(unit v)                       v / |v|, the unit vector along a vec of fields
(check TEST value)             value, with a condition on what built it, which ink check reports where it fails
(error "text" [v ...])         a load error saying text, then each v, where it is evaluated
(profile-reach f)              how far a 2D profile reaches from its origin, proven (to 1/1000 of its size); one reading t is refused
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
at its segment's ends apart (a U-turn needs segments of at least 2R). A
path whose first and last points coincide is a loop. The warnings say
where it fails and how to fix it. For a path in a plane and a body
extruded square to it, R is the reach in the plane:

```lisp
; A flat strap along a path in the xy plane: its section is 1.2 mm thick.
(fn strap (pts width)
  (check (path-thickness pts >= 0.6)
    (extrude-z (- (polyline-distance-2d pts) 0.6) (* -0.5 width) (* 0.5 width))))
```

**Make fields boundable.** Interval arithmetic must be able to prove where
a shape is not, out to infinity (Where analysis looks). Products of
signs, selects between far-reaching pieces, infinite extrusions, a twist,
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

A claim that reads a sampled field is an estimate (Results: `samples`).
Sampled data cannot shape the model at load time (a `for` list, an `if`,
a count, a `(bounds f)` box): that is an error.

## Repetition: `for`

`(for var list form ...)` repeats forms, once per element of `list`
(evaluated at load time), with `var` replaced by it and `$var` in names by
its text. It nests, and works:

- at top level;
- in assembly, body and package form clauses;
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
                                       per name; (name v (tol lo hi)) sets the band too. --param overrides
(import "std")                         load the standard library (see: ink std)
(param name default)                   design parameter; override with --param name=v
(param name default (tol lo hi))       ... with a tolerance, proven over (Tolerances)
(def name expr)                        named value; names are immutable. A def whose value is a
                                       lambda is called by name, like a fn
(fn name (a b ...) body)               function, expanded when the file loads
(fn name (a rest ...) body)            ... whose last argument collects the rest, as a vec
(let ((a e1) (b e2)) body)             local names, bound in order; usable anywhere
(t-span A B)                           the motion: assert checks this file's claims that read t at
                                       every moment from A to B; interfere checks pairs that read t over
                                       the least span holding every file's (Motion)
(body name (shape f) (attr k v) ...)   a part: one shape and any attributes (Bodies)
(body (name (p default) ...) ...)      ... with arguments, which instances pass
(assembly name clauses...)             a tree of placed instances (Assemblies)
(assembly (name (p default) ...) ...)  ... with arguments
(assert name claim [options])          a checkable design claim; options (res r), (accuracy x) (Claims)
(solve name (unknowns ...) eqs...)     constraint solving (Solving constraints)
```

Packages add more forms (Data for other packages).

**Sizing a part.** A part is sized by its arguments: a file's size
params are its body's and root's arguments, defaulting to the params,
`(body (rail (length length)) ...)`. An instance passes another size,
`(instance r (p.rail 470))`, so a part at two sizes is two instances of
one import. An import's value sizes the whole file, its claims and
studies too: `(import "rail.ink" (length 470))`.

Write dimensions from the params, `(probe 0 0 (+ base-t 1))` not `(probe
0 0 6)`, so claims follow a new value. Built-in names (x y z t pi, and
those of operations and forms: `len`, `ln`, `at`, ...) cannot name a
param, def or fn.

**Reading order.** Each form sees only the names defined before it: in the
file, in files imported before it, and in std. A name read before its
definition is an error saying where it is defined. Solves are the
exception, once the names they share are declared first (Solving
constraints).

**Names are lexical.** A `fn` or `lambda` reads the names of its own file,
its imports and std, never its caller's. So a part's `(box ...)` is std's
box whatever its importer calls `box`, and a path a function reads
(`read-json`, `import-stl`) is relative to its own file. A file may define
a name std defines, once: it shadows std's for that file and its
importers, while std's functions keep reading std's. Any other name is
defined once. An import cycle is an error.

**Libraries.** An import not found beside its file is read from the
first directory of `$INK_PATH` (separated as `$PATH`'s are) holding it,
and named by its path there in errors, claims' sources and `ink check`'s
imports. A parts catalog is `parts/` in such a directory.

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

**Placement.** WHERE is any frame in this assembly's coordinates: `(at
...)`, `(axes ...)`, `(frame-of ...)`, a def. Without CFRAME the child's
own axes line up with WHERE: `(instance n c (at (rot-z 90) (move 10 0
0)))` turns the part, then moves it 10 mm along this assembly's x. CFRAME
names a frame the child declares, met face to face with WHERE, exactly:
origin on origin, z opposite z, x along x. An assembly exports a frame by
declaring one, usually `(frame tip (frame-of l2 tip))`.

**`frame-of` paths.** In an assembly's clauses a path starts at one of its
instances placed so far (`b`, `l1/servo`), in its coordinates. Elsewhere
it is a body (its coordinates) or a path from a top-level assembly
(`root/arm/l1`, world coordinates). Close a loop with a solve (Solving
constraints).

**Joints.** MOTION is the connection's one degree of freedom, about or
along an axis of WHERE, applied after the faces meet:

- `(turn z a)` spins about the face normal (a shaft, a hinge);
- `(slide x s)` slides in the face plane (a carriage on a rail);
- `(slide z d)` moves along the normal.

The value is an expression of params and t: a constant is a fixed
setting, a param a slider, `(* 360 t)` motion. Other fixed offsets go in
WHERE: `(at (move 0 0 -8) (frame-of m shaft))`. One motion per instance.

**Arguments.** An assembly takes arguments as a function does, each with
a default: `(assembly (motor (angle 0)) ...)`. Its clauses read them like
local names. An instance passes values in order, `(instance m (motor (*
3000rpm t)) ...)`, reading params and t but not x, y or z; omitted ones
take their defaults. A part's interface is its argument list: catalog
motors, steppers and servos export `(drive (angle 0))`, rod actuators
`(drive (extension 0))`.

Instances of a body or assembly passing equal values share one; other
values build it again, named after the instance path that first passed
them (`root/m`, so its joint's claim is `limits:root/m/rotor`). `ink
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
keep their current values there. So mass scales with volume. An
envelope body's `(attr mass M)` works the same way. A mass less than the
bodies with their own density weigh is an error; one that gives no body a
density is unused, and `ink check` warns. The derived density is an
estimate (Geometric queries: `ink mass`).

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
seat)` seats a screw h out from the tap, turned with the thread. Threads
placed by hand can be out of phase and overlap (an M3 screw in its nut,
turned 180 degrees, 0.08 mm deep). Faces seated by CFRAME need no zone:
touching faces are proven to the accuracy.

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
assembly's coordinates and their shapes, each in its own coordinates. An
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
nothing passes check, assert and interfere. Where holding matters, place
it on the frame of what holds it and claim it, `(touches a b)`.

**Interference.** `ink interfere` claims that no two instances overlap
outside declared contact zones, each pair a `disjoint` claim at a
claim's default accuracy, or `--accuracy x`; the answer's `accuracy`
says which. A pair passing only by that allowance says how deep it came
(`held`, `near`). Without root it checks every body, unless an
assembly places one: then it refuses; name one with `--target`.

`pairs` lists the pairs whose bounding boxes overlap, each row the
answer of its `disjoint` claim (Results) with the contacts it used. The
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
(assert apart (disjoint (shape-of root/stand) (shape-of root/phone)))
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
  (assert elbow-on-link (touches (shape-of root/e) (shape-of root/upper)))
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
(assert reach (disjoint (shape-of root/key) (shape-of root/axis)))
(assert swing (clearance (shape-of turning/key)
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

`(assert NAME CLAIM [(accuracy x)] [(res r)])` states a requirement.
`ink assert` decides it by the cheapest method that does: a witness
point, an interval proof refined only where undecided, or a measurement.

```
(empty f)                     f has no interior
(disjoint a b)                a and b do not overlap
(subset a b)                  a lies inside b
(within f (x0 y0 z0) (x1 y1 z1))
                              f lies inside the box (corners: constants and params)
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

A claim is a form written in an `assert`, not a value, so a function
cannot return one; it may compute what a claim reads (its shapes, points
and bounds). `for` writes many claims: `(for p (vec (vec 1 2 3) (vec 4 5
6)) (assert solid-$p (subset (probe (nth p 0) (nth p 1) (nth p 2)) B)))`.

**Accuracy.** A claim's accuracy is how far its answer may be from the
truth, not a tolerance (how far a real part may differ from its drawing:
Tolerances). For a claim about shapes it is a length: 0.05 mm unless the
claim states `(accuracy x)`. `holds` and `fails` are shown under the
claim's accuracy rule; `unknown` means ink could show neither, and is
never a wrong `holds` or `fails`. The rule, by claim:

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
  It fails only where the shapes do not meet: a press fit 0.04 deep at
  0.05, or faces that touch, is `unknown`.

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
below 0 there, else outside, so a point on a surface is outside. Its
`value` is the field's bounds there; it is `unknown` only within rounding
of the surface.

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
  (assert light (integral root (attr density) < 2500))           ; mass
  (assert low (mean root z (weight (attr density)) < 80))        ; centre of mass height
  (assert over-foot (mean root (vec (- x 40) (- -40 x) (- y 30) (- -30 y))
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
  by running the form, so the claim is an estimate; a run that failed
  makes it `unknown`, with the reason in `unavailable`.

**Bounds and comparisons.** A claim's bound (the `v`) is any expression
of constants and params, at the current `--param` values (solved params
at their solved values): `(min-thickness part > (* 2 nozzle))`. It may
not read x, y, z or t, and only a value or mean claim's may read a ranged
param.

A comparison (`value`, `volume`, `integral`, `mean`, a metric) may state
an accuracy in its own units: how far past its bound it may be (a void
count's is a length on its shape: Manufacturing checks). `(assert
reach (value angle < 90) (accuracy 0.01))` holds where `angle - 0.01 <
90` everywhere. Generated value claims (`limits`, `speed-limits`, a
relation's `claim`) take one too. An exact value claim that meets its bound (a joint at
its limit) over a motion or a ranged param may stay `unknown`; with an
accuracy it decides. Points and ranges are exact and refuse `(accuracy
x)`.

### Results

`ink assert` reports each claim's `status` (`holds`, `fails`, `unknown`)
and `accuracy`:

- a number: guaranteed to it, under the claim's meaning (0: exactly).
  Points, ranges and value claims are exact by default.
- `null`: `unknown`, or an estimate with no guaranteed bound: package
  metrics, integrals, a claim measured at moments of a motion, and any
  claim that reads a sampled field or a solved value (not a certified
  root), directly or through defs, bodies and placements.

Its other fields say what the answer rests on:

| field | when | meaning |
|---|---|---|
| `witness`, `witness_pair` | a shape claim fails, except by being apart (a failing `overlaps`, or a `touches` too far off, has none) | the violating point (or pair). For an emptiness claim (the `depth` row's), every point within the accuracy of `witness` violates the claim. A failing clearance guarantees only its pair: two points, one in each shape, `gap` apart. A probe's is `point` |
| `value` | fails; a number claim | what the claim measured: a clearance's `gap`; a shape claim's field at the witness, negative (`-depth`); a value, volume, metric, integral or mean claim's number (a mass is an integral). Of several fields, each is under `integrals` or `means`, with its status |
| `gap` | a clearance fails | the distance between `witness_pair`: the first pair found closer than the clearance, so the smallest gap may be smaller (`ink distance` measures it) |
| `depth` | an emptiness claim fails (`empty`, `disjoint`, `subset`, `range`, a clearance of 0) | how deep `witness` lies in what violates the claim (inside both shapes for disjoint; inside a and outside b for subset), as the fields measure distance: the overlap is at least this deep |
| `held`, `near` | an emptiness claim holds, guaranteed (over a motion, with t ranged; an `interfere` pair as placed, with its zones) | `held`: `outright` (nothing violates it) or `to the accuracy` (some cells held only by its allowance). Then `near`, where measured: `depth`, the deepest violating point found in those cells, by bounds of the field as placed there at one moment and value of the ranged params (to a µm: seated faces read 0), and where above 0 `at` and `params`, those values (t's over a motion) |
| `params`, `t` | fails over ranged params or a motion | the params' values and the moment there |
| `hint` | unknown, where it can say; measured over a motion | what would decide it: a smaller `(accuracy x)`, a finer `(res r)`, a probe. Over a motion, why ranging t did not, and that search's hint (`motion.t_ranged`) |
| `undecided_cells`, `undecided_volume`, `undecided_cell` | unknown, by a search that kept them | how many cells it left undecided, their total volume, and one of them (the cell that ended the proof, where one did): `centre` and `size` |
| `hunt` | unknown, after a hunt | the cells (a value's pieces) searched for a violation once the proof ended, none found (`cells`), and what stopped it (`stopped_by`): `cap` (as many as the proof tested), `budget`, `time limit`, or `refinement limit` (every cell left at it) |
| `accuracy_asked` | unknown | the accuracy it was undecided at |
| `accuracy_attained` | unknown shape claim, where found | the accuracy it does hold to: `(clearance a b 0.29)` unknown with 0.003 means `(clearance a b 0.287)` holds. `accuracy_attained_in_window` where only the window was examined |
| `unbounded`, `in_window` | a shape not proven bounded | only the window was examined, and what was settled there |
| `res`; `res_given`, `res_given_by` | a claim searched on cells (not a value or a probe); where a res is given | the finest cells the search used, coarse where their bounds decided it before refining; the cell size asked, and by `assert` or `--res` |
| `samples` | reads a sampled field | each field's source, where, spacing, its error bound where one can be derived (`grid-bake`'s; none for `import-stl`, `redistance` or a package form's result) and the error measured at probe points (`measured_error`) |
| `solved_inputs`, `unsettled` | reads a solved value | the solves it rests on; those short of their accuracy, and why |
| `param_spans` | reads a ranged param | the params and their ranges; not a `(ranged)` no def names (`hidden`'s views) |
| `motion` | over a span | how the motion was checked (Motion) |
| `unavailable` | unknown, not checked | why: a solve with no value, or a model outside ink's range (docs/CORE.md Answers; interfere, distance and bounds refuse there) |
| `source` | always | `file:line` of the claim |
| `ms`, `cached` | always | wall time; whether it was read from the cache, nothing computed (false where no cache is needed: a param's value) |

The summary gives `all_hold` (every claim holds, guaranteed or estimated),
`all_hold_guaranteed` (every claim holds to a number) and `counts` by
status. With no claim checked both are null (`no_claims`), not a vacuous
true. The exit status: 0 when every claim holds (an estimate too), 3 when
one fails or is unknown, 1 on an error (no answer; `ok` false). `ink
assert FILE --own` checks only the file's own claims, not its imports'
(`"scope": "own"`, `imported_not_run`); `--name NAME` checks one.

`ink check` lists what the file defines and each import in brief
(`--imports`: all they define). Warnings cover every body, imported or
not. `ink check --load` skips surveying the bodies' surfaces: a quick
look at a large model.

**The cache.** The kernel's answers (proofs, searches, masses, meshes,
solves, `ink check`'s survey) are cached in `~/.cache/ink`, keyed by the
question as the kernel sees it and a hash of the ink binary: an unchanged
question is answered from the cache whatever else changed (`cached`); a
rebuilt binary recomputes everything. `INK_NO_CACHE=1` disables it.

**`ink mcp`.** The MCP server offers every command but `mcp` as one tool.
A call stopped at its `timeout_s` answers with what finished, marked
`stopped`; the same call with a larger `timeout_s` goes on from the
cache. An answer too long for clients comes back shortened, still whole
JSON; `shortened` names the file holding all of it.

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

- `--res r` (to `assert`, `interfere`, `distance`) lets the search halve
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
  their floor, a wrap read by its distance from 1/2 (std's thread-iso).
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

The claim's own field must be bounded: `(empty (max x ball))`, a
half-space cut by a ball, proves outright. A union is bounded by its
members' bounds, an intersection by its cheapest member's. To limit a
claim to a region, intersect with it.

On a shape not proven bounded, `ink bounds` and `ink mass` say
`unbounded`; `ink distance` proves a lower bound only where a target is
bounded; `ink interfere` checks it against every other item
(`unbounded_items`); `ink ray` searches only to the window's edge;
load-time `(bounds f)`, `grid-bake` and `redistance` are errors.

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
  search missed; the warning says how many cells it searched). Another
  body used in this one (`shape-of`) counts as one construction.
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
  nothing near them shows beside air (no surface crossing them, no air
  read as 0 by bounds); drawn as air, a wall may be drawn there where
  pieces meet flush.
- `ink mass [--density NAME|NUMBER]`: each target's mass, volume and
  centre of mass; with `--keys` naming them, its `inertia` tensor about
  the centre of mass, `principal_moments` and `principal_axes`. `total`
  sums the targets.
  - The density is the bodies' `density` attribute (Bodies) unless
    `--density` names another attribute or gives a number.
  - A target with no density has null `mass` and `center_of_mass`, its
    volume measured; `total.without_mass` names such bodies, `total.no_mass`
    what would give them one.
  - Mass, volume and a claim's integrals are estimates: cells start at a
    100th of the body's size (a claim's: each part's) and halve until the
    `spread`, a sampled estimate of the change left (twice what the last
    halving moved it), is within 2% of the
    mass (a mean's weight), or leaves a claim's status as it is.
    Each target gives `res`, `spread` (relative), and the cells before
    (`coarse_res`) and the mass's move from them (`relative_change`);
    one that did not settle says so (`unsettled`), as does every answer
    reading a published mass's density that did not. `--res` gives the
    cells instead.
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
argument. `render --out` writes PNG. `ink spec` and `ink std` take no
model file.

## Meshes

`ink mesh` (STL, 3MF or GLB), `ink render` and `ink components` work from
a triangle mesh of each target. A mesh approximates the field, and its
report says how far to trust it.

- **Sampling.** The field is sampled on a grid of cells `--res` wide and
  meshed by manifold dual contouring: watertight, sharp edges and corners
  kept, a thin wall's two sides apart. A feature thinner than res can go
  missing. Every cell is shown by bounds to hold nothing the mesh does not
  (one sign throughout, or the surface one sheet across it), or searched
  finer as far as a budget allows; each report counts every cell holding
  material or a gap it leaves out (`material`; `thinnest` across, where it
  is not shown drawn), every one undecided (`undecided`: nothing ruled out
  down to samples `searched` apart), and of them those beside a cell the
  surface crosses (`beside`), with one `at`, and `regions`: boxes covering
  them (every cell in one; touching bricks of 8^3 cells one region, more
  than 16 merged into larger covers), each with its counts, the coarsest
  spacing searched, the thinnest found and a point in it; and warns. Where
  the surface crosses a cell, one sheet over the cell grown by half a cell
  (where its vertex may be placed) counts as drawn: the surface there lies
  within that box's diagonal, 2 sqrt 3 res (`drawn_within`), of the
  cell's vertex, and so does what else lies in that box; elsewhere only a
  sign kept throughout does. Undecided cells where a select switches
  between two fields, which bounds cannot show continuous (a sweep's
  mitres), are counted apart (`seams`). A mesh whose grid
  would take more bricks of cells than its budget pays for is made at
  twice the res until it does not (`refused`: the res asked).
  `--max-res R`: no coarser than R.
- **A mesh for a file** (STL, 3MF): res is the target's size / 256, no
  finer than all the targets' / 2048 (mesh a small part beside a large one
  alone, or give `--res`). It is simplified to `--accuracy` (default
  res / 30; `0`: uniform), topology kept and every vertex of the uniform
  mesh within it. `components`: size / 160, made again finer where it
  leaves material out.
- **A mesh for seeing** (GLB, render): each body's res gives its surface
  as many cells as a cube 128 across has, coarsened so the whole model
  fits a triangle budget, and made again finer under the thinnest material
  it leaves out, most material drawn per triangle first, while the budget
  lasts (a render's, no finer than half a pixel); `--res` replaces it. It
  is simplified to `--accuracy` (default res / 4); shading turns at the
  shape's creases only.
- **The report** (each object of `ink mesh`): `res`;
  `uniform_triangles` (before simplifying) and `triangles` (in the file);
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
  - It draws the GLB's meshes. `legend` gives each target's res; it warns
    where thinner features went missing.
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
| overhangs, printed along +z | `(empty (overhangs part z0 a h))` | guaranteed, to a few degrees (below) |
| enclosed voids | `(cavities part < 1)` | proven to the accuracy, or unknown (below) |
| a flat outline to cut | `ink slice --out FILE.dxf` | its loops (Geometric queries) |

examples/machining-access.ink, mould-undercut.ink, print-checks.ink and
sls-hollow.ink show each holding and failing. Bridges are not checked: a
bridge a printer can span fails as an overhang where wider than 2 h tan a.

- **Walls.** `(min-thickness f > t)`, a shape metric of the design
  package: rays along the inward normal from surface samples; `value` is
  the thinnest, `at` where. Corners take their faces' thickness, so sharp
  edges are not thin. A thin spot between samples is missed: `(res r)`
  sets how fine a feature it sees, and walls no sample falls on make it
  `unknown`.
- **Gaps.** `(air-around f margin)` is the air out to margin past f's
  box. Its thickness is f's gaps, slots and holes narrower than margin,
  so give a margin above the gap claimed. A through hole's length is no
  gap. Between separate parts, `(clearance a b d)` is exact.
- **Line of sight.** `(hidden f stock band views)` is the points of band
  with material of f ahead of them along every view (a direction `(vec
  dx dy dz)`), looking past the stock's far side.
  - The band is `(clear-air f stock lo)`, the air in the stock lo clear
    of f: what a 3-axis mill must remove, or what fills a mould's box.
  - `(empty (hidden ...))` holds where every point can be reached from a
    setup, or pulled along the mould's pull (views +z and -z; a side
    action is one more view). A failure's witness is a hidden point.
  - Only signs are read, so the answer does not depend on how the field
    is written.
  - Each view's ranged param (Tolerances) multiplies the search: one
    claim over all views can stop undecided where a claim per setup, its
    band cut to what that setup reaches, decides.
  - A line of sight has no width: cutter fit is not checked.
- **Overhangs**, built up along +z on a bed at z0. `(overhangs f z0 a h)`
  is the air with material h above it and none h tan a away in any of 8
  directions: what lies under a face leaning more than a degrees past
  vertical. Only signs are read, so it holds for any field.
  - A pass means no face leans more than asin(2 e cos a / h) past the
    limit, e the accuracy: at 45, h 1, e 0.05, 49 degrees may pass and 52
    fails.
  - At 45, a face up to 2.3 degrees under the limit may fail, between the
    8 directions.
  - A ceiling narrower than 2 h tan a counts as held.
- **Enclosed voids.** `(cavities f < n)` counts f's enclosed voids, air
  no path leads out of (trapped powder or resin), with air deeper than
  the accuracy (a length, as for shapes: Claims). A shape metric of the
  design package, counted by interval bounds on cells halved where they
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
    side) is air: cells across it split at a face square to an axis, or
    where both pieces read one face. At the budget (or forecast past it)
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
may read it: `(instance arm arm (at) (turn z (* 360 t)))` turns once a
second. A joint's rate is in deg/s (turn) or mm/s (slide); `3000rpm` is
18000 deg/s, and a speed in rpm is a multiply, `(* r 1rpm)`. `--t T`
evaluates any command at one moment.

`(t-span 0 1)` in a file makes `ink assert` and `ink interfere` check
every claim and pair that reads t at every moment of the span;
`--t-span 0:1` or `--t T` replaces every file's span. A claim takes its
own file's span, imported or not (`ink check`: `import_t_spans`); an
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
`unknown`. What a solve on t sets is ranged with t, enclosed over each
piece of the span by the steps it is followed by (Solving constraints).
Not ranged: a claim reading a solve on t that is not followed there, or
over ranged params; a metric, integral or mean; a body travelling more
than 256 times its own size over the span (a chain's pin: shorten the
span); in `interfere`, a shape whose own field reads t (a belt's teeth).

- **Name the moving part first.** Two shapes (`disjoint`, `clearance`,
  `subset`, `overlaps`, `touches`) are searched in the frame the first
  rests in, so gears or a carriage on its rail are one search and a pin
  turning in its hole does not move: a pin's `subset` of its hole, a
  link's `disjoint` of the sprocket. Two shapes riding one moving
  instance are asked once. Witnesses are given in the world. A part
  turning about a point the motion holds still, though placed through
  another angle (a coupler on its rocker's pin, met by the loop's
  closure), is searched turning about that pivot, its bounded drift
  charged to the claim (a clearance d asked at d plus it; an emptiness at
  its accuracy less it). A turn rounded off in floats (a CAD converter's)
  may leave a pair measured.
- **Speeds** (`speed-limits`, `diff` along t) are the formula's
  derivative, blind to a jump: where a drive may jump (a `select` or
  `floor` on t) the speed is `unknown`. A turn's speed limit reads the
  turn, so its angle wrapped by `mod` at a whole turn (`(mod (* 360 t)
  360)`) is decided. Write profiles with `min`, `max` and `clamp`.
- **Repeating motions** are searched over one period: a claim reading t
  only through a drive (t, or one value all its reads of t pass through:
  a spin-up's angle) and, through it, only through the sine and cosine
  of turns at constant rates of it, over a span where the drive travels
  past their common period, taken from their exact values (`(* 3600 t)`
  and `(/ (* 2 drive) 3)` share a turn; floats such as 0.3 and 0.7 share
  only a period of about 1e16 s). A shape turned about its own z that is
  the same turned by 1/n of a turn (a gear's teeth) repeats whenever the
  drive turns it so: a gear pair at mesh-turn by one pinion tooth, its
  rounding (about 1e-13 mm) charged to the accuracy. Through a drive
  other than t only a hold is an answer; else t is ranged.

The answer's `motion` block gives `how` (`t ranged`, `asked once` or
`measured`), the `undecided` t ranges, and for `measured` the moments
checked (`steps`, and `step` or a package's `basis`). A failure gives its
moment, `t`.

```lisp
(import "std")
(t-span 0 1)
(body crank (shape (move (box 40 6 4) 20 0 0)))
(body post (shape (move (cylinder-z 3 20) 35 35 0)))
(assembly root
  (instance post post)
  (instance c crank (at) (turn z (* 90 t)) (limits -10 100) (speed-limits 0 100)))
(assert clear (clearance (shape-of root/c) (shape-of root/post) 1))
```

`clear` and the joint's `limits:root/c` and `speed-limits:root/c` are
checked over the whole turn.

## Data for other packages

Packages add top-level forms, values, metrics, commands and clauses. The
`design` package owns `solve`, `optimize` and the metrics of
Manufacturing checks. Others (simulation, dynamics) come in distributions
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
naming it (a misspelling would drop its check).

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

- `ink assert` and `ink interfere` prove every claim for every value of
  each ranged param (one with a band) it reads; other commands work at
  the values. `--param bore-d=8.1` moves the band with the value;
  `--param k=lo:hi` gives an absolute range; `--param k=v:v` pins k;
  `--nominal` drops the declared bands. `ink check` lists each param's
  `tol` and `band`, and of one an import sets, its `written` default and
  `set_at`.
- `(ranged lo hi [at])` is a range of the question, not the part (how far
  a tool is withdrawn, every angle from 0 to 90). It is ranged as a band
  is, but `--nominal` keeps it, and `assert` and `interfere` refuse
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
    count, `(bounds f)`, a sampled field, a package form's settings);
  - a param a `solve` sets;
  - a direction not proven nonzero over the bands and the motion;
  - a bound that reads one, except in value and mean claims.
- Ranges are boxes: a round positional zone written as dx, dy in [-r, r]
  is a square. Write round zones in polar params.
- Statistical tolerancing (root-sum-square, Monte Carlo) is not covered.

## Solving constraints

```lisp fragment
(solve loop
  (unknowns (coupler-angle 0.8) (rocker-angle 1.6))        ; names and initial guesses
  (== (+ (* crank (cos crank-angle)) (* coupler (cos coupler-angle)))
      (+ ground (* rocker (cos rocker-angle))))
  (== (+ (* crank (sin crank-angle)) (* coupler (sin coupler-angle)))
      (* rocker (sin rocker-angle))))
```

- Unknowns become params whose values satisfy the equations (least
  squares, Levenberg-Marquardt, exact derivatives). Claims on solved
  params are estimates (`"accuracy": null`), except over ranged params,
  where they are certified roots (below).
- Equations are `(== lhs rhs)` or expressions that must be 0, of params
  and `t`, not x, y or z. One form covers sketch constraints, assembly
  placements and linkage kinematics.
- Systems are solved at load time (t = 0) and again, warm-started, for
  each `--param` value, in the order of what they read, not the file's. A
  solve may read what a later one sets, through a name declared before
  both (`(param b 0)`, then `(unknowns b)` in the later). Two solves that
  read each other's unknowns are refused: give them as one.
- Unknowns may name existing params (`(unknowns tilt-x px)`), closing a
  loop over an assembly's frames: `(== (nth gap 0) 0)`, `gap` the
  difference of two `(get (frame-of ...) "origin")`.
- A solved value cannot shape the model at load time (a `for` list, an
  `if`, a count): use it in fields. Solved params cannot be overridden.
  `ink check` reports solved values and residuals.
- Where a solve does not solve (at the `--t` given), `assert` and
  `interfere` refuse what rests on it (`unknown`, its error as
  `unavailable`; an instance under `unsolved_items`) and decide the rest;
  `check` warns; other commands fail with the error.
- **Derivatives** through a solved value are the implicit function's, so
  a joint driven through a solve has its real speed. At a dead point (no
  regular root) what reads them is refused.
- **Following along t.** A solve that reads t is, at each moment, the
  root followed from its load-time root, not whichever one the solver
  finds there, so a linkage stays on its assembly.
  - Where the followed root meets another or the equations switch (a dead
    point, a `select` on t), every later moment is refused, saying where,
    even where both sides agree. Start the motion where the linkage
    closes, away from a dead point. So, in that call, is every moment
    past where following had got when the call's time limit passed.
  - A solve whose equations repeat in t (t only inside sin and cos, its
    unknowns angles or not repeating) is followed one period; any other
    moment is that period's root, its angles moved by whole turns.
  - A solve with more equations than unknowns has no certified root, so
    is not followed: at each moment it is the solver's answer there, as
    is every solve on t that reads what it sets, and claims reading them
    are measured.

### Over ranged params

A solve whose equations read a ranged param is solved over the ranges:
`ink assert` and `ink interfere` take each unknown to be the one root
near the solution at the ranges' middles, certified for every value of
the params, so claims on it are guaranteed; and so are its rates, a
`diff` through it.

```lisp
(import "std")
(param L 50 (tol -0.1 0.1))   ; a link's length
(param h 30 (tol -0.1 0.1))   ; the rail its end rides on
(solve place (unknowns (theta 37) (sx 40))
  (== (* L (cos (* theta (/ pi 180)))) sx)
  (== (* L (sin (* theta (/ pi 180)))) h))
(assert x-high (value sx <= 40.2))  ; holds (exactly): sx = sqrt(L^2 - h^2) <= 40.1995
```

A root resting on an uncertified solve or a sampled field is an estimate
(`ink check`'s `over_bands_note` names them). Refused, naming the solve:
no root certified (two roots close together, a singular Jacobian, or
ranges too wide: narrow them); more equations than unknowns; an unknown
used at load time. A system that reads `t` is certified at one `--t`;
over a span of t it is `unknown` (its `unavailable` says why).

## Optimization

```lisp fragment
(optimize lighter-bracket
  (vary plate 3 8) (vary web 6 20)                 ; params and their bounds
  (minimize (mass (shape-of bracket) 2.7g/cm3))    ; or (maximize ...)
  (require (min-thickness (shape-of bracket) > 2)) ; any claim
  (require (clearance (shape-of screw) web-region 1))
  (budget 30))                                     ; evaluations
```

- Each param may appear in `(vary ...)` once.
- Objectives, as in claims: `(volume f)`, `(mass f density)`,
  `(value expr)`, a shape's metric (`(min-thickness f)`, `(cavities f)`),
  a package's metric of a form. One with no value at a design fails that
  design, saying why: a volume, mass or shape's metric of a shape not
  proven bounded, a metric with none (a failed solve), one measured on a
  grid that cannot see the part.
- `ink optimize FILE [--name NAME] [--budget N]` runs it (`ns.o` for one
  imported under a namespace). The budget is a soft cap.
- The search steers by each requirement's measured violation at each
  design (`search_violation`; `search_feasible` where all are at most 0),
  which decides nothing. Over a motion it measures at one moment (`note`
  says which).
- At the design returned, each requirement's `answer` is its claim asked
  as `ink assert` asks it, over its file's `(t-span)` unless `--t` selects
  one moment, and over the declared
  tolerances. It is `satisfied` exactly where that answer holds;
  `feasible`: every requirement satisfied. The optimize's answer as a
  whole is an estimate (`accuracy: null`). A run whose answer is its start
  design warns, as does one whose design misses a requirement.
- It reports that design (exact: its values read back as the floats its
  claims were asked of), its history, its `method` and `apply`: the
  options that ask `ink assert` what it was verified for, each param in
  the unit it was given in (in a run stopped early's answer too). Re-check
  the design with them.
- A design the search tried that could not be evaluated (one the model
  refuses there, a study that fails) is in `failed_designs`, with its
  error, and warned of: the search keeps clear of it, so the design
  returned may be limited by it, not by a requirement. A `"gradient"`
  run stops at one instead and answers its error.
- `method` is `"gradient"` where every objective and requirement has
  exact derivatives (volume, mass, `value`, a package metric with an
  adjoint), else `"derivative-free"`, which needs more evaluations.
- Strength, stiffness and heat are not optimized: this release has no
  simulation.

## Standard library

std is written in ink (std/std.ink, std/materials.ink).

- `ink std` lists its sections and definitions' signatures.
- `ink std --section NAME` prints, with comments, the section whose
  banner (heading comment) is NAME, or the definitions it names (`gear`:
  every gear builder); `--all` prints all of it.
- `ink spec --section NAME` finds this reference's sections the same way,
  by title or path (`Solving constraints/Over ranged params`), else looks
  in std and the text of both.
- A call of an unknown function suggests std's name for the CAD word
  (`hole`: subtract a `cylinder-z`; `union`: `min`; `offset`: `-`;
  `extrude`: `extrude-z`) and near spellings.

Sections: scalar helpers (`clamp`, `mix`, `blend`, `atan`,
`acos`, `asin`); combination (`subtract`, smooth blends, fillets,
chamfers, shells); appearance (`paint`, `color-map`); warps (`move`,
`rotate-z`, `scale`, `twist-z`, `extrude-z`, mirrors,
repeats and patterns, thread frames); primitives (2D and 3D, exact
distances); lattices; manufacturing; printing; sheet metal; paths; belts
and chains; power transmission (gears, threads, pulleys).

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
  refused, saying where. Of straight edges only, it is exactly its
  corners: a corner is dropped only exactly on the line between its
  neighbours, and edges that touch are refused. With arcs or rounds, a
  corner or an arc's bulge within 1e-9 of its points' bounding-box
  diagonal of the line that would replace it is dropped, and that line
  used (drops in a row add up); beside a shared corner, a meeting whose
  sliver between the two edges is within 100 times 1e-9 of the diagonal
  (a tangent joint's rounded data) is allowed.
- `(polygon pts drafts z0 z1)`: a moulded or drafted wall. pts is the
  profile at z = 0; each edge's wall leans outward by its draft (degrees;
  one for all, or one per edge) as z rises; rounds keep their radius;
  straight beyond z0 and z1.
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
  toward +z at inside radius r.
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
- Not built: crossing bends, closed corners, hems. A bend runs straight
  across its flange.

```lisp
(import "std")
(def S (sheet (polygon
 (vec (vec 0 0) (vec 60 0) (vec 60 40) (vec 0 40)))
              (vec (move (circle-2d 2.5) 48 20 0)) 2 0.4
              (vec (sheet-bend (vec 30 40) (vec 30 0) 90 2))))   ; x > 30 folds up
(body bracket (shape (sheet-folded S)) (frame upright (sheet-frame S 48 20 -1)))
(body flat (shape (sheet-flat S)) (mark bend (sheet-bend-lines S)))   ; the cut file
(assert flanges (empty (sheet-short-flanges S 4)))
(assert reliefs (empty (sheet-unrelieved S)))
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
  `straight-bevel-gear`. `(mesh-distance za zb module)` and `(mesh-turn
  za zb bearing a)` place a pair. In these, and in `gear-mesh` and
  `rack-mesh`, a gear (za, zb, z) may be `(vec teeth centre-teeth)`, one
  cut to run at another count's centre distance; `bevel-mesh` refuses
  it.
- **Gears in mesh** are relations, A the input:
  - `(gear-mesh za zb module [pa [backlash [efficiency]]])`: spur, each
    axis its gear's z; the efficiency is its coupling's, read by a
    dynamics package. External pairs only: a planet in its ring fails
    `centre-distance`;
  - `(bevel-mesh za zb module face ...)`, `(rack-mesh z module ...)`.
  - pa is the pressure angle (degrees, default 20). Under 10 is an error
    (likely a backlash written in its place).
  - Its claims say the gears sit as a mesh, pitch surfaces parted by no
    more than the backlash allows: `centre-distance`, `parallel-axes`,
    `faces-meet`; a bevel's `pitch-cones`, `square-axes`, `apexes`; a
    rack's `pitch-line`.
  - It has no contact zone: interfere checks the teeth over the pinion's
    angle, so a wrong ratio or phase is a clash. It does not check the
    bodies have the teeth it says, unless they clash.
- **Belt and screw relations**: `belt-drive`, `chain-drive`,
  `belt-span`, `lead-screw`, `held`, `tensile-stiffness`: couplings read
  by a dynamics package; an ink without one refuses them (Data for other
  packages).
- **Threads**: `(thread-iso d p)`, an ISO metric external thread along
  +z. Its field is the distance to the profile, so `(- (thread-iso d p)
  clearance)` is a tapped hole's void. Pulley grooves
  (`gt-pulley-teeth-2d`, `htd-pulley-teeth-2d`) are nominal, with no
  running clearance.

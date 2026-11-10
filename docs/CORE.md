# ink

ink is an implicit geometry kernel for AI agents. A shape is where a
function f(x, y, z) is negative; its surface is where f changes sign.
Operations are arithmetic on these functions, so they cannot fail. Every
spatial question becomes one: is this field ever negative, and where?

## Principles

- **Agent-usable.** Agents write, read and fix models without human CAD
  concepts. Errors say how to fix the file.
- **Fast.** Editing and querying feel interactive, so agents can iterate.
- **Robust.** Operations never fail, and answers are never wrong about what
  they claim.
- **General.** The core has general mechanisms. Domain knowledge (threads,
  gears, materials, parts) lives in libraries.
- **Minimal.** Nothing agents don't need, and no two ways to do one thing.
- **Durable.** Models keep working as ink evolves. A part is checked when it
  changes, and stays checked.

Balance these; don't maximise one at the others' expense.

- **No topology.** No faces, edges or B-rep structure. Agents name the
  pieces they build (`tag`) and the frames parts connect by.
- **No legacy in the core.** Translation to legacy formats is a separate
  layer at the edge, built only when needed.

## Assemblies and motion

- Parts declare named **frames** (a mounting face, a shaft axis).
  Assemblies connect parts, and place sub-assemblies, frame to frame.
- A **joint** is a connection that turns or slides, driven by a
  parameter or by time.
- A connection marks **intended contact**, which clash checks ignore.
- A moving assembly is a field of x, y, z and time t, and time is searched
  like any range. Where the search decides a claim over a span, a hold
  covers every moment, and a failure needs only one moment that violates
  it. Elsewhere the answer is an estimate from sampled moments; a failure
  found at one is guaranteed only where the claim, in exact real numbers,
  fails at that moment.
- A jump in a drive is unbounded speed: speed claims are unknown there,
  never held.

Parts carry properties as data (`attr`: heat output, maximum temperature,
keep-out distance) and requirements, checked in assemblies. Assemblies
need no physics; a physics package would read the same frames and
properties. This release has none.

## Answers

Every answer says whether the claim holds and how far off it may be: its
**accuracy**. An answer with one is checked everywhere to it. One with
none is unknown or an **estimate** (sampled, or solved numerically), with
no guaranteed bound. Parts have tolerances, answers an accuracy, and ink
picks the resolution (how finely it looks).

- **holds**: shown true under the claim's accuracy rule.
- **fails**: shown false under that rule, with a point violating the
  claim where one can show it (none shows shapes failing to touch or
  overlap by being apart).
- **unknown**: undecided at the accuracy asked.

A positive clearance is exact on its distance: a gap under it never
holds. A claim that shapes are apart, that one contains another or lies
within a box, or that they touch, allows a violation no deeper than the
accuracy. A claim is asked once, to its accuracy, so a piece no thicker
than twice the accuracy can pass as absent. A claim stating no accuracy
is asked to the finest tolerance of the parts it is about, a part
declaring none counting at the default.

A param can take a range: a tolerance band, or values to search. Over a
range, a claim holds only where it is shown for every value. An overlap or
a contact found at one value says nothing of the others, so those claims
answer only `fails` (apart at every value, or a contact overlapping too
deep at one) or `unknown`.

Integrals (volume, mass, centre of mass) are estimates, never proven,
sampled over a range. An estimate says what it rests on: resolution,
samples, solved values.

A guaranteed answer is never wrong to its accuracy, given one assumption:
the `libm` crate's transcendentals (sin, cos, tan, atan, exp, ln) are
accurate to a few ulps (steps between adjacent floats), the margin
ink-core's interval arithmetic widens them by
(crates/ink-core/src/interval.rs).

Answers are about each shape as ink builds it: exactly its formula, except
a `polygon` with drafted arcs, or degenerating within its allowance of a
drafted span's end, or with an arc or round it cannot build exactly. Its
edges, arcs, rounds and drafted walls are otherwise exact, a circle read
within bounds on its exact centre and radius (their rounding: near 1e-14
of the circle, more at a sharp corner), a drafted corner's speed within
bounds of its rounding, and the height where a drafted profile's
construction changes within a float step or so: an answer holds for the
drawn shape, and is unknown where those bounds decide it. Each
simplification in building one of the exceptions (merging nearly straight
corners, say) stays within a ten-millionth of its size of the drawn
corners (10 nm on a 100 mm part), but merges add up.

### Range

ink's guarantees hold within 1,000 km of the origin. The finest length a
claim may ask about is a nanometre, growing farther out to a billionth of
the distance from the origin: a micrometre at 1 km, a millimetre at
1,000 km. Distance is judged by the farthest coordinate ink searches,
past the shapes (a box 500 mm out may be judged 1.5 m out), so a model
near either limit can be refused.

Every number a model is built or asked with must be 0 or 1e-100 to 1e100
in size: one written in the file or read from one, a parameter or a
moment, or one ink's own code makes to build the model. std's own
constants are exempt. A value computed from these alone may be any size,
but one that floats lose (an underflow or overflow) is out of range.
Over a tolerance band or time span only the ends are judged; values
between are bounded outward, sound at any scale.

Outside the range a claim answers `unknown`, saying why, as does every
claim of a model built from an out-of-range number. `ink interfere`,
`distance` and `bounds` refuse, saying why.

An answer's numbers are written exactly; a warning's prose rounds them to
10 digits.

A check exists only if a wrong answer would change a decision. A part's
checks can often be cheap (a probe at a point, a coarse envelope), each
proven to its accuracy like any claim. A search over all of space earns
its cost where a defect could hide anywhere, such as a clash in an
assembly.

## Layers

```
core       fields, operations, evaluation, intervals
kernel     search, meshing, measurement, queries
io         formats, meshes for export and seeing, pictures
language   the .ink language: its reader and units, imports, std
packages   physics, design: above the core, never needed by it
clients    the CLI, viewers, the parts catalog, agents
edge       legacy translation, separate, only when needed
```

Each layer and each package is its own crate, so a lower layer cannot
depend on a higher one; the `ink` library gathers them for clients.

A package plugs in only through the language's **registry**, the list of
packages a build has (ink-lang's `ext.rs` and `command.rs`), adding:

- the forms it reads;
- the assembly clauses it owns, read with the assembly's names in view,
  placed with each instance;
- the metrics claims bound;
- CLI and MCP commands.

What no package in a build owns is unknown to it. The CLI holds nothing
of any one package, so a third party's plugs in as ours do. The kernel
builds in one package, design (the `ink` library's `registry()`). A
distribution with more builds on the library and the CLI
(`ink_cli::main`, given its registry). Nothing in the kernel names a
package outside it.

Something belongs in the core only if it is part of what a model is,
every client needs it, and it can be specified exactly.

## Fields

```
(+ a b ...) (* a b ...) (min a b ...) (max a b ...)
(- a b) (/ a b) (pow a b) (atan2 y x)
(neg a) (abs a) (sqrt a) (square a) (exp a) (ln a)
(sin a) (cos a) (tan a) (floor a)
(deg->rad a) (rad->deg a)  degrees in radians and back, pi exact
(emax a b ...)           Euclidean max: exact where the gradients are orthogonal
(select c a b)           a where c < 0, else b
(warp e fx fy fz [ft])   e with x, y, z (and t) replaced
(lipschitz k f)          f, with the promise |grad f| <= k
(tag name f)             f, naming a piece
x y z t, params, numbers, pi
```

Every operation defines:

- a value: an f64, made total (`x / 0` is 0, anything undefined is 0);
- a range, rounded outward so proofs are sound;
- a derivative along x, y, z, t and every parameter.

Growing a shape by d (`offset-shape`, std's `shell`) means everything within d
of it. ink builds that exactly where it shows the field is the distance
there, and warns elsewhere, saying what bound holds.

`lipschitz` is the one promise, and makes clearance proofs fast. A false
promise makes proofs wrong, so the language drops any promise a warp or a
varying size would break, and a proof takes none whose constant, as
computed for the query, is infinite or undefined.

A model is a finite expression graph: `def`, `fn`, `let` and `import`
expand by substitution, and a name means the definition visible where it
is written. Answers are cached by the question as the kernel sees it and
a hash of the ink binary: an unchanged question is answered once, and any
rebuild recomputes everything.

## Not doing now

- Face recovery, edges by name, face fitting.
- Certificates on frames; proofs about rounding.
- Legacy output, a native slicer, a GPU proof engine.

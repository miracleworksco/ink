# Core conformance suite

The tests every implementation of the ink core must pass: each
operation's value, derivatives, and range rule (its bounds over a box of
inputs). The rules are below; docs/IR.md, "Fields", has them in brief.

- `reference.py` is a second implementation, written from the rules
  alone, never from the kernel's outputs. It generates `cases.json`:
  `python3 conformance/reference.py > conformance/cases.json`.
- `cases.json` holds point cases (inputs, the value, the partial
  derivatives) and range cases (an input box, and the hull of values and
  derivatives over a dense sample of it, which a range rule must contain).
  They cover composed expressions (norms, min/max ties, select, pow,
  atan2, nested warps), expressions of time and a param, a sampled field
  (cell interiors, faces, the grid's edge), and extreme floats. The extreme
  cases' inputs are rounded to f32 before the reference evaluates them.
- `crates/ink/tests/facade/conformance.rs` runs the kernel's
  evaluators (points, dual numbers, gradients, intervals, interval
  gradients) against it, and checks that the symbolic, dual and gradient
  derivatives agree bit for bit. Composed and sampled cases' values allow
  64 f64 epsilons of relative error, an operation's at most four; range
  containment has no tolerance. A sampled hull shows a range
  rule unsound where it misses a value, never that it is sound or tight.
- Its `KNOWN` list is where the kernel does not meet the rules, each with
  its repro. One that passes fails the test as `FIXED`, so none goes
  stale.

To change the rules, change them here (and docs/IR.md's brief), then
reference.py; regenerate, and fix the implementations until the suite
passes:

```sh
cargo test --profile fast -p ink-cli --test cli facade::conformance::
```

## The rules

The operations are docs/IR.md's Fields table.

### Values

Every op is total: it has a value for every input, the same on every
implementation, CPU and GPU: its IEEE-754 double result, except:

- There is one zero: -0 reads as 0.
- `x / 0` is 0.
- Anything else undefined is 0: `inf - inf`, `0 * inf`, `ln` of a
  negative, `sin` of infinity, `pow` of a negative base with an exponent
  that is not a whole number.
- `sqrt` of a negative is 0 too: the continuous `sqrt(max(a, 0))`,
  not reported.

Infinities are values, written `inf`. Overflow gives +-inf, `ln 0` is
-inf, and `pow(0, b < 0)` is +inf.

IEEE-754 requires no particular result for `exp`, `ln`, `sin`, `cos`,
`tan`, `atan2`, `pow` and `emax`'s hypot. On the CPU each is one
implementation, the same on every machine: FreeBSD's, through Rust's
`libm` (under an ulp for `exp` and `ln`, nearly rounded for the others),
as are the kernel's packages'. So a model, and every answer about it,
means the same numbers everywhere.

`ink check` warns where its sampled survey finds a body's field passing
through an undefined case, naming the construction: a 0 in a shape means
"on the surface", so one can draw surface that is not there. One between
the samples goes unseen. `ink eval` lists them per point under
`undefined`. One in a number fixed as the model loads (a count, a
condition, a claim's bound) is reported at its line. Only cases the value
depends on count: not a branch not taken, the losing side of a min or max,
or a factor of a constant 0 (a zero-length segment's distance is a
point's, though its projection is 0 / 0). They do not affect proofs: a
range that includes an undefined case covers every value it could take,
so a proof through one is at worst unknown.

In proofs, an expression of constants and params alone means its f64
value (a param is an exact number for a query), except where it reads a
param given a range (a tolerance, a `ranged`). Then it means its value in
real numbers, bounded over the range, even at a single value of it, where a witness
(what shows a claim fails) is confirmed (docs/IR.md, Tolerances).

### Derivatives

Derivatives (`grad`, `diff`, `ink eval`, an optimizer's) are each op's
formula under the same rules, on the branch the value took.

- At an exact tie (`min` or `max` of equal values, `abs` at 0) they are
  the average of both sides, so a sharp edge gets the direction between
  its faces. Variadic `min`/`max` are nested pairs, so at a tie of three
  or more arguments the weights follow the nesting.
- `select`'s derivative is that of the branch its condition takes;
  `floor`'s is 0.
- Through an expression they are the chain rule, op by op: each op's
  partials, from its operands' values (`1 / b` and `-(a / b) / b` for a
  quotient, `0.5 / sqrt a` for a root), times its operands' derivatives,
  summed in order. A tie's average is `da/2 + db/2`, halved first.
- A sampled field's read is an op of its coordinates whose partials are
  the reads one axis further there.

`diff` writes exactly this, and `ink eval`'s gradient and the optimizer
compute it, so all agree bit for bit. A solved value's derivative along t
or a param is the implicit function's (docs/IR.md, Solving constraints):
each adds, to the partial with every param held, each solved value's
partial times its rate, the same rates, so through a solve they agree to
rounding.

Proofs bound a derivative by its formula's range. For `diff` and `grad`
along x, y and z they also use the bounds on the field's own gradient
that clearance proofs use:

- The slope of `sqrt(u1^2 + u2^2 + ...)` along x is at most
  `sqrt(sum (du_i/dx)^2)`; `(dot a a)` counts, as std's `norm` writes it.
- A `min` or `max` takes the hull of its sides' slopes; `emax` a
  combination of them with weights in [0, 1].
- So `(diff (len2 x y) x)` is within [-1, 1] over any box, even across
  the skeleton, where its formula is a quotient of two ranges through 0.
- Along an axis none of an op's operands varies along, its slope is 0,
  even where the op jumps.

These bounds hold the derivative's value in real numbers, with the tie
and quotient conventions above: theorems, taking no promise. A
promise bounds the slope of f where f has one:
`(select (neg (square x)) x (* 1000 x))` is x, but its derivative at
x = 0 is 1000. A warped derivative, no longer the slope of what it is
read on, gets the formula's range only.

`diff` along t of a field that reads t keeps what it differentiates
(through another `diff` and a warp too). So a claim on it over a motion
also checks a joint's drive for jumps, and leaves a def's switches
unchecked (docs/IR.md, Motion).

### Sampled fields

A sampled field is values on a regular grid that travels with the model,
read by trilinear interpolation. The core op is the interpolant inside
the grid box and its eight components: the value, and the derivative
along each subset of x, y, z (an axis twice is 0 inside a cell).
Coordinates clamp into the grid box; a point on a cell face belongs to the
cell above it. Outside the box the continuation is plain ops: the
boundary value plus distance. A package's result that is a quantity, a
stress or a displacement, holds its boundary value instead. Sampled data
enters a model only this way. There is no external-field callback, since
code could lie about its ranges.

### Precision

Proofs do not depend on precision. Any evaluator whose ranges round
outward may prove; lower precision only makes unknown more frequent.
Point values below f64 are approximate. f32 (a GPU's) is for viewing: it
agrees with f64 on the side of the surface wherever f64 is clear of it by
more than a small tolerance of the model's size. In f32, the ops that
would give infinity (division, `exp`, `pow`, `ln`, `tan`) give the largest
finite f32. Sums and products beyond f32's range, `pow` with a small whole
exponent, and subnormals are outside what f32 approximates.

# ink's simulation packages

Packages built on the kernel (../kernel/crates) as a third party's would be,
in a workspace of their own. The kernel names none of them; `ink-ext`
puts them in one `ink`.

| Crate | What it adds | Reference |
|---|---|---|
| ink-sim, on ink-fem | `(study ...)`: strength, stiffness, vibration and heat by the finite cell method on the field, no mesh; `ink solve`; study claims | [SIM.md](docs/SIM.md) |
| ink-topology | `(topology ...)`: a stiff layout within a volume budget; `ink topology` | [TOPOLOGY.md](docs/TOPOLOGY.md) |
| ink-dynamics, with ink-carried | gravity, actuators' effort along a motion, reactions, modes, simulation; `ink effort` | [DYNAMICS.md](docs/DYNAMICS.md), [validation](docs/DYNAMICS-VALIDATION.md) |
| ink-effort | ink-dynamics' effort as a tool of the claims core (../kernel/crates/ink-claims): an actuator's peak torque or force and its speed there, and its torque and speed at each moment the run read, so a rating whose limit moves with speed is held against the motion at every one; facts read off one run, made when one is asked for | its crate's first lines |
| ink-study | ink-sim's solve as a tool of the claims core: a part's peak stress and deflection under a study's loads, each with the band three grids show, as facts read off solves made when one is asked for | its crate's first lines |
| ink-ext | `ink`: the kernel's command line and MCP server with all of the above, and no parts catalog: `ink parts` is a distribution's that holds one | `ink spec` |
| `ink require` | what a file's `(require ...)`s answer (JSON, or words with `--plain`), cheapest source first: what the design states, its parts' tables and formulas, then the two tools above (a motion's efforts at every moment read, a part's stress and deflection under its studies), each started only by a fact nothing before it settled; and every `(require ...)` claim, each saying over how much of what it reads it was shown | ../kernel/crates/ink-cli/src/cli/require.rs's first lines |

## Build

<!-- not run: it builds the binary the tests run -->
```sh
cargo build --release          # target/release/ink
```

Then `target/release/ink spec` is the language reference with each
package's after the kernel's, and `ink mcp` serves it to an agent as the
kernel's does.

## Try

Start with a bar on two supports, pressed at its middle
([examples/beam.ink](examples/beam.ink)): is it stiff enough, and strong
enough? Its comments work both out by hand. It takes seconds on a laptop:

<!-- run -->
```sh
target/release/ink require examples/beam.ink
```

```json
{"requirements": [{"name": "stiff-enough", "verdict": "yes", "value": 2.94},
             {"name": "strong-enough", "verdict": "no", "value": 0.79}]}
```

The bar bends 2.94 mm, under the 4 allowed (by hand, 2.93), and its
safety factor is 0.79, under the 1.5 asked (by hand, 0.80). A claim holds
or fails only where three grids, each twice as coarse as the last, agree on
it; where they don't, it is `unknown`, with a `hint` saying what would
decide it. `ink solve examples/beam.ink` gives the whole solve: stress,
displacement, reactions and its buckling factor.

[examples/](examples/) has a model for each package: a bracket's layout
(layout.ink, topology.ink), an arm's motor torque (arm.ink), a servo
arm's efforts (servo-arm-effort.ink), a simulated arm (simulated-arm.ink), a heat
sink under a fan, its film coefficient given (heatsink.ink), and a
lattice shelf's study (shelf.ink: on a part this
intricate its claims are `unknown`, and say why).
[benchmarks/](benchmarks/) checks each package against hand calculations.

<!-- not run: an example -->
```sh
target/release/ink require examples/slider-loads.ink          # a carriage's loads and its motor's effort: all hold
target/release/ink solve examples/shelf.ink --study books-5kg  # a lattice shelf: minutes, warning what its cells miss
```

## Test

<!-- not run: minutes -->
```sh
./cpu-test.sh                  # every test
./cpu-test.sh --slow           # with the tests too slow for every change
```

## How far to trust it

Every answer says how it was found: a study's claim holds only where
three grids agree, and an estimate says it is one. Each reference says
what its package does not model and where its answers stop holding: read
it before trusting a number for something that matters.

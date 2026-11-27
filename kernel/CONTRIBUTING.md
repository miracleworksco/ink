# Contributing to ink

## Build

```sh
cargo build --release              # target/release/ink
cargo build --profile fast         # target/fast/ink: release speed, quick to rebuild
```

## Layout

Each layer is its own crate, none depending on a layer above it
(docs/CORE.md):

- **ink-core**: what a shape is: the expression graph, interval
  arithmetic, values and derivatives.
- **ink-kernel**: what can be known about it: one branch-and-bound search
  over space, params and time, plus meshing, measurement and queries.
- **ink-io**: STL, 3MF and GLB export, STL import as a field, PNG renders.
- **ink-lang** and **std/**: the `.ink` language and its standard
  library, written in ink.
- **ink-mma**: a gradient optimizer, for packages (sim's topology).
- **ink** bundles them as a library, with one package of its own (the
  printing and machining checks: thinnest wall, enclosed voids); **ink-cli** is the `ink` command and
  its MCP server.

The simulation packages are in `../sim/`, a workspace of their own built on
the kernel as a third party's would be: studies (**ink-sim**, on
**ink-fem**), **ink-topology**, **ink-dynamics** (with **ink-carried**), and **ink-ext**,
an `ink` with them. `../sim/cpu-test.sh` runs their tests.

Also here: `examples/` (small focused models), `benchmarks/` (engineering
tasks with hand calculations) and `conformance/` (an independent reference
implementation of the field rules).

## Test

```sh
cargo test                         # the everyday suite
cargo test -p CRATE NAME           # one crate's tests whose names contain NAME, while you work
cargo fmt                          # the format (rustfmt.toml)
```

`scripts/test.sh` runs the everyday suite with the test binaries in
parallel. The public repository's CI also checks formatting, runs clippy
and builds the release binary on Linux and macOS.

A crate's tests are in that crate: unit tests beside the code, the
library's (source to answer) in crates/ink/tests, and the command line's
in crates/ink-cli/tests/cli, which run the `ink` binary cargo builds.
Their `models` test checks every model in examples/, demos/, benchmarks/
and tests/fixtures except the slow ones and any needing a package the build
lacks. A model passes when `ink check` gives no warnings and its claims
hold. Other expectations go in its file, a line each:

- a claim meant to fail: `; expect-status: NAME fails` (`fails estimate`
  where its answer has no accuracy);
- a warning: `; expect-warning: TEXT`.

A model's inputs (an STL it imports) are committed beside it.

Answers are cached in ~/.cache/ink, keyed by the question and a hash of
the `ink` binary, so the first run after a code change recomputes them.

The slow set (`benchmarks/run.sh`) is what is too slow for every change:
the slow benchmarks, the guide's answers, the heavy models (`SLOW` in
models.rs) and each test marked `#[ignore = "slow: owed for CHANGE"]`,
which names the change that owes it a run. Run it for such a change, and
before a release.

## A good change

- Something that failed now works, with a test that would catch it
  breaking again.
- docs/IR.md (`ink spec`) says what to write and what an answer means,
  and changes with what it describes. A change's history goes in its
  commit message, validation numbers in benchmarks/, method in the code,
  an open problem in an issue.
- It deletes code where it can, and keeps what agents need.
- Its commit message says what now works.

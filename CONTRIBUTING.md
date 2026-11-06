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
- **ink-design** (with **ink-mma**): constraint solving, optimisation
  over params, and the printing and machining checks (thinnest wall,
  enclosed voids). It is a package, plugged in as any other would
  be.
- **ink** bundles them as a library; **ink-cli** is the `ink` command and
  its MCP server.

Also here: `examples/` (small focused models), `benchmarks/` (engineering
tasks with hand calculations) and `conformance/` (an independent reference
implementation of the field rules).

## Test

```sh
cargo test                         # everything, every time: about 1 min on 32 threads once built
cargo test -p CRATE NAME           # one crate's tests whose names contain NAME, while you work
cargo fmt                          # the format (rustfmt.toml)
```

Run `cargo test` and `cargo fmt` before you push. CI also runs `cargo fmt
--check`, `cargo clippy` and the release build, on Linux and macOS
(.github/workflows/ci.yml).

A crate's tests are in that crate: unit tests beside the code, the
library's (source to answer) in crates/ink/tests, and the command line's
in crates/ink-cli/tests/cli, which run the `ink` binary cargo builds.
Their `models` test checks every model in examples/, benchmarks/ and
tests/fixtures except the slow ones and any needing a package the build
lacks. A model passes when `ink check` gives no warnings and its claims
hold. Other expectations go in its file, a line each:

- a claim meant to fail: `; expect-status: NAME fails` (`fails estimate`
  where its answer has no accuracy);
- a warning: `; expect-warning: TEXT`.

A model's inputs (an STL it imports) are committed beside it.

Answers are cached in ~/.cache/ink, keyed by the question and a hash of
the `ink` binary, so the first run after a code change recomputes them.

The slow benchmarks and the guide's fine-resolution answers
(`benchmarks/run.sh`) are validation, not tests: run them when a change
could move what they measure.

## A good change

- Something that failed now works, with a test that would catch it
  breaking again.
- docs/IR.md (`ink spec`) says what to write and what an answer means,
  and changes with what it describes. A change's history goes in its
  commit message, validation numbers in benchmarks/, method in the code,
  an open problem in an issue.
- It deletes code where it can, and keeps what agents need.
- Its commit message says what now works.

## Values

- One general mechanism over special cases, knobs or heuristics.
- Fix the cause, so the whole class of bug can't recur.
- Purity: a shape is one function built from simple fields with min and
  max. What can't be expressed so yet is a stated limitation, not a hack.
- No filler: docs, comments and answers hold facts, or nothing.

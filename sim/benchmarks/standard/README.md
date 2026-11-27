# Public standard verification cases

This is a measurement set for review, not a claim of parity with a vendor.
Each `.ink` header cites its input facts, specifies the asked quantity and
location, and identifies the target's meaning, with the public source it
was taken from. Vendor prose, figures, meshes and decks are not included.

Run from this folder's simulation workspace (`sim/`), using its own freshly
built binary, copied so a rebuild cannot change it mid-run:

```sh
cargo build --profile fast
mkdir -p ~/standard-sets
cp target/fast/ink ~/standard-sets/ink-frozen
python3 benchmarks/standard/run.py \
  --ink ~/standard-sets/ink-frozen --ink-commit "$(git rev-parse HEAD)" \
  --out ~/standard-sets/default
```

The runner executes **one case at a time**, cold (`INK_NO_CACHE=1`), under a
16 GiB memory ceiling and no swap. Each case gets 600 seconds including
loading. The runner rejects grid, support-margin and time-step overrides.
Study durations, materials, loads and supports are problem inputs. Transient
steps are omitted, retaining ink's default 100. There is no refined track in
this measurement. `--resume` reuses only the same binary/model hashes and
deadline. `--reparse` re-extracts retained raw JSON without rerunning a solve;
it preserves the original timing and model hash, and explicitly records that
postprocessing occurred. Do not run two runner processes at once.

Per-case directories retain check/answer stdout, stderr and `result.json`;
`results.json` contains the selected set. The measurement records both the
source commit used to build ink and the benchmark checkout's current head,
plus binary/model SHA256 hashes. Solver source must match `--ink-commit`.

Interpret the fields carefully:

- `status` is ink's measurement-envelope claim status, or `unknown` when a
  point/reaction has no accuracy band. It does **not** mean the target passed.
- `band` is ink's explicit low/high band, or value ± its reported `grid_error`
  for a scalar claim. An explicit buckling band already includes grid error;
  the runner never adds it again. These are numerical **estimates**, not
  guaranteed error bounds. A missing band stays null, even for an exact match.
- `qualified`/`scored` says whether the target and solid problem are comparable.
  An unqualified idealisation cannot become a solver pass or miss. Its raw
  difference and containment, if available, are diagnostic only.
- `band_contains_target` is raw numerical containment; null means no band.
  `band_width_relative` is total width divided by absolute target (or an
  explicit engineering scale for a zero target). Temperature percentages
  depend on the chosen absolute temperature scale; raw Celsius differences
  are retained and are the primary thermal errors here.
- The policy fixed before the batch is: **confidently wrong** means a qualified
  target outside a decided (`holds`/`fails`) estimated band whose total width
  is at most 10% of the target. This operational threshold is not a proof.
  Lower stress/deflection and higher frequency/buckling predictions are
  marked unsafe. Unknowns, missing bands, unqualified targets and broad
  bands cannot satisfy this definition.
- Printed rounding is recorded separately; containment tests the stated
  central target, not an enlarged tolerance chosen from the measured answer.

Analytically equivalent extrusions and boundary substitutions are identified
in model headers. They do not cover the omitted formulation (e.g. TPLA01's
mixed prescribed-flux boundary). Point loads on a 3D solid have singular
loaded-point displacement; a shell's finite point displacement is not a valid
solid target. Problems with conflicting source inputs/targets are quarantined.

Runner regression tests:

```sh
python3 benchmarks/standard/test_runner.py
```

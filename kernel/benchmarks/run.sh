#!/usr/bin/env bash
# The slow set, which `cargo test` leaves out: first every test marked
# `#[ignore = "slow: owed for ..."]` in a crate's tests/ (the slow
# benchmarks' key answers that are not claims, and the tests too slow to
# run per change, each saying the change that owes it a run), then the slow
# models' own claims by the model rule (SLOW in
# crates/ink-cli/tests/cli/models.rs: benchmarks and heavy examples, over
# ~2 s each cold). The quick ones run with `cargo test`.
# usage: benchmarks/run.sh [NN ...]   (default: the whole slow set; NN: those benchmarks, slow or not)
set -euo pipefail
cd "$(dirname "$0")/.."
# A run takes several GB: where systemd runs user scopes it is capped at
# 16 GB (a run over it is killed, not the machine's other work), and where
# flock is, runs that share the lock file (INK_HEAVY_LOCK) go one at a
# time.
run=()
cap=(systemd-run --user --scope -q -p MemoryMax=16G -p MemorySwapMax=0)
if command -v systemd-run >/dev/null && "${cap[@]}" true 2>/dev/null; then
  run=("${cap[@]}")
fi
lock="${INK_HEAVY_LOCK:-${XDG_CACHE_HOME:-$HOME/.cache}/ink-heavy.lock}"
if command -v flock >/dev/null && mkdir -p "$(dirname "$lock")" 2>/dev/null; then
  run=(flock "$lock" "${run[@]}")
fi
failed=0
run_tests() { ${run[@]+"${run[@]}"} cargo test --profile fast "$@"; }
cli() { run_tests -p ink-cli --test cli "$@"; }
if [ $# -eq 0 ]; then
  run_tests --workspace --test '*' -- --ignored || failed=1  # every crate's tests/: unit tests' ignored ones are calibrations, run by hand
  models=slow  # SLOW in crates/ink-cli/tests/cli/models.rs, as the test reads it
else
  tests=(); files=()
  for n in "$@"; do
    match=(benchmarks/"$n"-*.ink)
    if [ ! -e "${match[0]}" ]; then
      echo "benchmarks/run.sh: no benchmark $n (benchmarks/$n-*.ink): the numbers are in benchmarks/README.md" >&2
      exit 2
    fi
    tests+=("benchmarks::b$n"); files+=("${match[@]}")
  done
  cli -- --include-ignored "${tests[@]}" || failed=1
  models="${files[*]}"
fi
INK_MODELS="$models" cli -- models:: || failed=1
exit $failed  # both parts run, so every failure is reported

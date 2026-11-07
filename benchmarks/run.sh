#!/usr/bin/env bash
# The slow models, which `cargo test` leaves out (validation, not tests):
# first the slow benchmarks' key answers that are not claims (benchmarks.rs's
# ignored tests), then the slow models' own claims by the model rule (SLOW in
# crates/ink-cli/tests/cli/models.rs: benchmarks and heavy examples, over
# ~20 s each cold). The quick ones run with `cargo test`.
# usage: benchmarks/run.sh [NN ...]   (default: every slow one; NN: those benchmarks, slow or not)
set -euo pipefail
cd "$(dirname "$0")/.."
# A run takes several GB: where systemd runs user scopes it is capped at
# 16 GB (a run over it is killed, not the machine's other work), inside the
# systemd slice INK_SLICE when given (the catalog watch's,
# scripts/catalog-watch.sh, where no cap is no run), and where flock is, runs
# that share the lock file go one at a time.
run=()
cap=(systemd-run --user --scope -q ${INK_SLICE:+--slice=$INK_SLICE} -p MemoryMax=16G -p MemorySwapMax=0)
if command -v systemd-run >/dev/null && "${cap[@]}" true 2>/dev/null; then
  run=("${cap[@]}")
elif [ -n "${INK_SLICE:-}" ]; then
  echo "benchmarks/run.sh: FAILED: no memory-capped scope in $INK_SLICE" >&2; exit 1
fi
lock="${XDG_CACHE_HOME:-$HOME/.cache}/ink-heavy.lock"
if command -v flock >/dev/null && mkdir -p "$(dirname "$lock")" 2>/dev/null; then
  run=(flock "$lock" "${run[@]}")
fi
cli() { ${run[@]+"${run[@]}"} cargo test -p ink-cli --test cli --profile fast "$@"; }
if [ $# -eq 0 ]; then
  cli -- --ignored benchmarks::
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
  cli -- --include-ignored "${tests[@]}"
  models="${files[*]}"
fi
INK_MODELS="$models" cli -- models::

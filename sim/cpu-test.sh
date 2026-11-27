#!/usr/bin/env bash
# Every test here (none needs a GPU): area.toml's test, and the public CI's
# simulation job. Built once, then the test programs run side by side, each
# printed whole as it ends.
#
# usage: cpu-test.sh [--slow]
#   --slow   the tests too slow for every change too (each marked
#            `cfg_attr(not(feature = "slow"), ignore = ...)`): area.toml's slow
#            set. `--quick`, the default, is taken for the public CI.
set -euo pipefail
cd "$(dirname "$0")"
args=(--workspace --lib --bins
  --test sim --test cli --test topology --test dynamics --test plug)
[ "${1:-}" = --slow ] && args+=(--features ink-ext/slow,ink-dynamics/slow)
out=$(mktemp -d)
trap 'rm -rf "$out"' EXIT
# Each test program, as its package directory and its path.
cargo test --no-run "${args[@]}" --message-format=json-render-diagnostics | python3 -c '
import json, os, sys
for line in sys.stdin:
    m = json.loads(line)
    if m.get("reason") == "compiler-artifact" and m.get("executable") and m["profile"]["test"]:
        print(os.path.dirname(m["manifest_path"]), m["executable"])
' >"$out/programs"
# Tests at a time per program: every program runs at once, and each test
# spreads over the cores itself.
# The slow set two at a time, one program after another: four of its
# tests at once, then its programs side by side two at a time, passed the
# 16 GB a heavy job is given (killed 2026-11-24, 16 GB at its peak).
slow2=4 alone=; [ "${1:-}" = --slow ] && slow2=2 alone=1
threads=${TEST_THREADS:-$slow2}
pids=()
names=()
while read -r dir exe; do
  name=$(basename "$exe")
  (
    start=$SECONDS
    status=0
    (cd "$dir" && "$exe" --test-threads="$threads") >"$out/$name.log" 2>&1 || status=$?
    echo "== ${name%-*}: $((SECONDS - start)) s" >>"$out/$name.log"
    exit $status
  ) &
  pids+=($!)
  [ -z "$alone" ] || wait "$!" || true
  names+=("$name")
done <"$out/programs"
failed=()
for i in "${!pids[@]}"; do
  wait "${pids[$i]}" || failed+=("${names[$i]%-*}")
  cat "$out/${names[$i]}.log"
done
[ ${#failed[@]} -eq 0 ] || { echo "FAILED: ${failed[*]}"; exit 1; }

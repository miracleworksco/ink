#!/usr/bin/env bash
# ext's benchmarks: each file's own
# claims by the model rule (ink-ext/tests/cli/models.rs, given the file in
# INK_MODELS), then
# its key queries against the hand calculations; pass/fail and wall time per
# file, as the kernel's benchmarks/run.sh reports. First the simulation
# benchmarks' key answers that are tests (bNN_ in
# ink-ext/tests/cli/benchmarks.rs, the slow ones ignored by cargo test).
# usage: benchmarks/run.sh [PATTERN]   (e.g. benchmarks/run.sh 13)
#   INK_EXT=/path/to/ink   ext's binary (default: target/fast/ink; INK is the kernel's)
set -u
here="$(cd "$(dirname "$0")" && pwd)"
export INK_EXT="${INK_EXT:-$here/../target/fast/ink}"
export INK="$INK_EXT"
[ -x "$INK" ] || { echo "no ink binary at $INK: cargo build --profile fast in $here/.., or set INK_EXT"; exit 1; }
cd "$here"
pattern="${1:-}"
out="$(mktemp -d)"; trap 'rm -rf "$out"' EXIT
pv() { python3 -c "import json,sys,math; d=json.load(sys.stdin); print($1)"; }
py() { python3 -c "import json,sys,math; d=json.load(sys.stdin); sys.exit(0 if ($1) else 1)"; }
# A solve's {columns, rows} probes as records (as benchmarks/run.sh reads
# them): one probe's field, probe INDEX KEY [COMPONENT].
TABLE='def rec(t):
    out = []
    for row in t["rows"]:
        r = {}
        for c, v in zip(t["columns"], row):
            *up, k = c.split("."); o = r
            for u in up: o = o.setdefault(u, {})
            o[k] = v
        out.append(r)
    return out
d = json.load(sys.stdin)
for k, x in list(d.items()):
    if isinstance(x, dict) and set(x) == {"columns", "rows"}: d[k] = rec(x)'
probe() { python3 -c "import json,sys
$TABLE
p=d['probes'][$1]['$2']; print(p if '${3:-}'=='' else p[${3:-0}])"; }
near() { python3 -c "import sys; v,r,t=map(float,sys.argv[1:]); sys.exit(0 if abs(v-r)<=t*abs(r) else 1)" "$1" "$2" "$3"; }
pass=0; fail=0
tests=("benchmarks::b${pattern}")
if (cd .. && cargo test -q --profile fast -p ink-ext --test cli -- --include-ignored "${tests[@]}"); then
  pass=$((pass + 1)); echo "PASS  ${tests[*]}"
else
  fail=$((fail + 1)); echo "FAIL  ${tests[*]}"
fi
for f in [0-9][0-9]-*.ink; do
  [ -n "$pattern" ] && [[ ! "$f" =~ ^$pattern ]] && continue
  t0=$(date +%s)
  if ! m=$(cd .. && INK_MODELS="benchmarks/$f" cargo test -q --profile fast -p ink-ext --test cli -- --exact models::every_model_passes_the_model_rule 2>&1); then
    echo "FAIL  $f: $(echo "$m" | grep -v '^ *$' | tail -5 | tr '\n' ';')"; fail=$((fail + 1)); continue
  fi
  ok=0; detail=""
  case "$f" in
    *) ok=1 ;;
  esac
  if [ "$ok" = 1 ]; then pass=$((pass + 1)); printf 'PASS  %-28s %5ss  %s\n' "$f" "$(( $(date +%s) - t0 ))" "$detail"
  else fail=$((fail + 1)); printf 'FAIL  %-28s %5ss  %s\n' "$f" "$(( $(date +%s) - t0 ))" "$detail"; fi
done
echo "$pass passed, $fail failed"
[ "$fail" = 0 ] && [ "$pass" -gt 0 ]

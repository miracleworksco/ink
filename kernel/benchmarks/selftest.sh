#!/usr/bin/env bash
# The slow models' runner, checked with a fake cargo and systemd-run:
# benchmarks/run.sh runs uncapped where it cannot cap, and runs both its
# parts whichever fails. Seconds.
#
# usage: benchmarks/selftest.sh   (RUN_SH=PATH: another benchmarks/run.sh to check)
set -u
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
wrong=0
expect() { # pass|fail LABEL COMMAND...: COMMAND, in a subshell of its own, must pass or fail
  local want=$1 label=$2 got; shift 2
  if ( "$@" ) > "$T/out" 2>&1; then got=pass; else got=fail; fi
  if [ "$got" = "$want" ]; then echo "ok    $label"; else echo "WRONG $label: wanted $want, got $got"; sed 's/^/      /' "$T/out"; wrong=$((wrong + 1)); fi
}
mkdir -p "$T/bin"; printf '#!/bin/sh\nexit 0\n' > "$T/bin/cargo"; chmod +x "$T/bin/cargo"  # a cargo that builds and passes at once

# benchmarks/run.sh (RUN_SH, default beside this file)
RUN_SH="${RUN_SH:-$(dirname "$(readlink -f "$0")")/run.sh}"
if [ -f "$RUN_SH" ]; then
  B="$T/bench"; mkdir -p "$B/benchmarks"; cp "$RUN_SH" "$B/benchmarks/run.sh"
  printf '#!/bin/sh\nexit 1\n' > "$T/bin/systemd-run"; chmod +x "$T/bin/systemd-run"
  expect pass "slow models: no cap runs uncapped" env PATH="$T/bin:$PATH" XDG_CACHE_HOME="$T" bash "$B/benchmarks/run.sh"
  mkdir -p "$T/rbin"; cp "$T/bin/systemd-run" "$T/rbin/"
  cat > "$T/rbin/cargo" <<'SH'
#!/bin/sh
echo "$*" >> "$RAN"
case " $* " in *" $FAIL "*) exit 101 ;; esac
SH
  chmod +x "$T/rbin/cargo"
  bothrun() { # FAIL: cargo's runs matching FAIL fail; run.sh still runs both parts, and fails
    : > "$T/ran"
    ! env PATH="$T/rbin:$PATH" XDG_CACHE_HOME="$T" RAN="$T/ran" FAIL="$1" bash "$B/benchmarks/run.sh" &&
      grep -q -- '--ignored' "$T/ran" && grep -q 'models::' "$T/ran"
  }
  expect pass "slow models: failing slow tests still let the models run, and fail the run" bothrun --ignored
  expect pass "slow models: failing models fail the run" bothrun models::
else
  echo "WRONG no benchmarks/run.sh at $RUN_SH"; wrong=$((wrong + 1))
fi

[ $wrong = 0 ] && echo "slow models' self-test: ok" || { echo "slow models' self-test: FAILED ($wrong wrong)"; exit 1; }

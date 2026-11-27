#!/usr/bin/env bash
# Calibrate the view budget's ratios by hand, on the instructions the
# processor counts (Linux perf events): Tape::work's per-instruction costs
# (crates/ink-core/src/tape.rs) and `work`'s TRIANGLE_WORK and READ_WORK
# (crates/ink-cli/src/cli/queries.rs). The tests hold the kernel's own
# counts to what this last printed, on any machine; a change that makes a
# mesh or an evaluation cost otherwise than the kernel counts (a faster
# evaluator, a new pass) is caught only by running this again.
#
# usage: scripts/view-budget-calibrate.sh [MODEL.ink ...]
#   Runs on one thread, uncached, on x86 with perf_event_paranoid at 2 or
#   below (it stops where the system does not let a process count). Prints
#   - ops::COUNTED's rows (crates/ink-cli/src/cli/counted.rs): what a
#     tape of each kind of instruction retires a point, a gradient, a box
#     and a simplification;
#   - VIEW_CORPUS's rows (queries.rs): the instructions each shape's view
#     mesh retires;
#   - for each model (default: examples/curves.ink,
#     examples/gyroid-heat-exchanger.ink; name others, a long extrusion say),
#     the instructions its `ink mesh M.glb` meshes retire and what the
#     budget charged them (`per_triangle`: near TRIANGLE_WORK where the
#     ratios hold), then each of its first meshes;
#   - TRIANGLE_WORK and READ_WORK fitted to every mesh, and each mesh's
#     counted over predicted, at the fit and as in use.
#   Paste the rows and constants that moved; Op::work's table is set by
#   hand from ops::COUNTED (each op's cost beyond its case's arithmetic).
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
[ "$(cat /proc/sys/kernel/perf_event_paranoid)" -le 2 ] || { echo "view-budget-calibrate: perf_event_paranoid is above 2: this system does not let a process count its instructions" >&2; exit 1; }
if [ $# -eq 0 ]; then set -- examples/curves.ink examples/gyroid-heat-exchanger.ink; fi
models=$(for m in "$@"; do realpath "$m"; done | paste -sd:)
export INK_NO_CACHE=1 INK_PATH="$root" INK_WORK_MODELS="$models"
cargo test -p ink-cli --lib -- --ignored --nocapture --test-threads=1 counts_instructions measure_op_work calibrate_view_work

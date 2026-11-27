#!/usr/bin/env bash
# The everyday tests, as `cargo test` runs them, but with the test binaries
# run side by side once built rather than one after another: a binary's
# tail (the cache's answers checked against cold ones; the model rule; the MCP time
# limits, which wait) then overlaps the others' work. Then the doc tests.
# Options go to cargo (`-p CRATE`); test filters are not taken.
# usage: scripts/test.sh [CARGO TEST OPTIONS]
set -euo pipefail
cd "$(dirname "$0")/.."
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
cargo test --no-run --message-format=json-render-diagnostics "$@" > "$T/built"
python3 - "$T/built" <<'PY'
import json, os, subprocess, sys, tempfile
bins = []
for line in open(sys.argv[1]):
    m = json.loads(line)
    if m.get("reason") == "compiler-artifact" and m.get("executable") and m["profile"]["test"]:
        # each in its package's directory, as cargo runs it
        bins.append((m["target"]["src_path"], os.path.dirname(m["manifest_path"]), m["executable"]))
runs = []
for src, cwd, exe in bins:
    out = tempfile.TemporaryFile("w+")
    runs.append((src, exe, out, subprocess.Popen([exe], cwd=cwd, stdout=out, stderr=subprocess.STDOUT)))
failed = []
for src, exe, out, p in runs:
    p.wait()
    out.seek(0)
    print(f"     Running {os.path.relpath(src)} ({exe})\n{out.read()}", end="", flush=True)
    if p.returncode != 0:
        failed.append(src)
if failed:
    sys.exit("failed: " + " ".join(os.path.relpath(f) for f in failed))
PY
cargo test --doc "$@"
# Beside it in the private repository, the viewer builds on the kernel's
# graph (a new node kind breaks it unseen): built, not tested, for every
# kernel change. Not shipped with the kernel, so not there for a release.
if [ -f ../viewer/Cargo.toml ]; then
    cargo check --manifest-path ../viewer/Cargo.toml --locked --bin ink-view
fi

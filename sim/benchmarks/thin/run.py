#!/usr/bin/env python3
"""Thin material: claims on walls, fins, webs, ligaments and split cells
thinner than the cells, each with what is true of it and how that is
known. A design of the solve or its reads is judged by how many of these
it answers wrongly.

    python3 benchmarks/thin/run.py INK [--slow] [--only REGEX] [--no-mirror] [-j N] [--threads N] [--json OUT]

Each case in cases/ is an ink model (`.case`, not `.ink`: the model rule
and the catalog require every claim in a .ink file to hold, and these
hold only where ink answers them wrongly) whose header gives, per claim,
`; truth: NAME holds|fails about|at least|at most VALUE (HOW)`: the claim's
true answer, the metric's value or a bound on it, and how that is known (a
closed form, statics, or the same model on finer cells, itself a case).
`; test: slow` marks a case that runs minutes; --slow includes them.
`; known: NAME (WHY)` marks a claim ink is known to answer wrongly, a
limit its rules state: its wrong answer is reported as KNOWN, apart.

Each claim is also asked the other way, on the same study (so on the same
solves), as `not-NAME`: the comparison reversed at the same bound, whose
truth is the opposite status on the same value. Most claims are written the
safe way (`max-stress <`, `min-safety-factor >`), so a reading that errs
toward failing them shows only in their reverse: a false fail of `>`.
--no-mirror leaves them out.

Per claim the runner reports ink's status and value against the truth:
  FALSE HOLD    ink holds a claim that is false;
  FALSE FAIL    ink fails a claim that is true;
  MISREAD FAIL  ink fails a false claim on a value more than 10 times off
                the true one ("about" truths): a number the part never reaches;
  unknown       ink declines (not wrong);
  ok            ink answers as the truth does;
  KNOWN ...     a wrong answer the case's header says is known.
and a summary line. It exits 1 when any claim not known wrong is a false
hold or fail (or misread), or ink gave no answer.

This is the simulation area's slow set (area.toml): every case, --slow
included, must answer as its header says, its declared `known:` wrong
answers apart.
"""
import argparse, concurrent.futures, json, os, pathlib, re, subprocess, sys, tempfile, time

HERE = pathlib.Path(__file__).resolve().parent
TRUTH = re.compile(r'^; truth: (\S+) (holds|fails) (about|at least|at most) (\S+) \((.*)\)$')
KNOWN = re.compile(r'^; known: (\S+) \((.*)\)$')
ASSERT = re.compile(r'\(require (\S+) \((\S+ \S+) ([<>]) ')


def mirrored(text, truths):
    """The case with each claim that has a truth also asked reversed, on
    the same study, and those claims' truths."""
    extra, more = [], {}
    for m in ASSERT.finditer(text):
        name = m[1]
        if name not in truths:
            continue
        end, depth = m.start(), 0
        while True:
            depth += {'(': 1, ')': -1}.get(text[end], 0)
            end += 1
            if depth == 0:
                break
        claim = text[m.start():end]
        op = '>' if m[3] == '<' else '<'
        extra.append(claim.replace(f'(require {name} ', f'(require not-{name} ', 1).replace(f'{m[2]} {m[3]} ', f'{m[2]} {op} ', 1))
        t = truths[name]
        more['not-' + name] = dict(t, status='fails' if t['status'] == 'holds' else 'holds', how=f"{name} reversed: {t['how']}", mirror=True)
    return text.rstrip('\n') + '\n' + '\n'.join(extra) + '\n', more


def cases(slow, only, mirror, tmp):
    for f in sorted((HERE / 'cases').glob('*.case')):
        text = f.read_text()
        if (' test: slow' in text and not slow) or (only and not re.search(only, f.stem)):
            continue
        known = {m[1] for m in (KNOWN.match(l) for l in text.splitlines()) if m}
        truths = {m[1]: dict(status=m[2], kind=m[3], value=float(m[4]), how=m[5], known=m[1] in known, mirror=False) for m in (TRUTH.match(l) for l in text.splitlines()) if m}
        if mirror:
            text, more = mirrored(text, truths)
            truths.update(more)
            f = pathlib.Path(tmp) / f.name
            f.write_text(text)
        yield f, truths


def run(ink, f, env):
    t = time.monotonic()
    p = subprocess.run([ink, 'require', str(f)], env=env, text=True, capture_output=True)
    try:
        out = {a['name']: a for a in json.loads(p.stdout).get('requirements', [])}
    except json.JSONDecodeError:
        out = None
    return out, p.stderr[-400:], time.monotonic() - t


def verdict(truth, a):
    if a is None:
        return 'error'
    st, v = a.get('status'), a.get('value')
    if st == 'holds':
        return 'FALSE HOLD' if truth['status'] == 'fails' else 'ok'
    if st == 'fails':
        if truth['status'] == 'holds':
            return 'FALSE FAIL'
        if truth['kind'] == 'about' and isinstance(v, (int, float)) and v != 0 and max(v / truth['value'], truth['value'] / v) > 10:
            return 'MISREAD FAIL'
        return 'ok'
    return 'unknown' if st == 'unknown' else 'error'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('ink')
    ap.add_argument('--slow', action='store_true')
    ap.add_argument('--only')
    ap.add_argument('--no-mirror', action='store_true')
    ap.add_argument('-j', type=int, default=3)
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--json')
    a = ap.parse_args()
    env = dict(os.environ, RAYON_NUM_THREADS=str(a.threads))
    env.setdefault('INK_PATH', str(HERE.parents[1]))
    tmp = tempfile.TemporaryDirectory()
    todo = list(cases(a.slow, a.only, not a.no_mirror, tmp.name))
    with concurrent.futures.ThreadPoolExecutor(a.j) as ex:
        done = list(ex.map(lambda c: (c, run(a.ink, c[0], env)), todo))
    rows, count, by = [], {}, ({}, {})
    for (f, truths), (out, err, secs) in done:
        for name, t in truths.items():
            got = out.get(name) if out is not None else None
            v = verdict(t, got)
            if t['known'] and v in ('FALSE HOLD', 'FALSE FAIL', 'MISREAD FAIL'):
                v = 'KNOWN ' + v
            count[v] = count.get(v, 0) + 1
            by[t['mirror']][v] = by[t['mirror']].get(v, 0) + 1
            val = got.get('value') if got else None
            rows.append(dict(case=f.stem, claim=name, truth=f"{t['status']} {t['kind']} {t['value']:.4g}", how=t['how'], status=got.get('status') if got else 'error',
                             value=val, at=got.get('at') if got else None, verdict=v, seconds=round(secs, 1), error=None if out is not None else err))
    w = max(len(r['case']) for r in rows)
    print(f"{'case':{w}}  {'claim':12} {'truth':24} {'ink':8} {'value':>12}  verdict")
    for r in rows:
        val = f"{r['value']:.4g}" if isinstance(r['value'], (int, float)) else ''
        print(f"{r['case']:{w}}  {r['claim']:12} {r['truth']:24} {r['status']:8} {val:>12}  {r['verdict']}")
    wrong = sum(count.get(k, 0) for k in ('FALSE HOLD', 'FALSE FAIL', 'MISREAD FAIL', 'error'))

    def summary(c):
        return (f"{c.get('FALSE HOLD', 0)} false holds, {c.get('FALSE FAIL', 0)} false fails, {c.get('MISREAD FAIL', 0)} misread fails, "
                f"{c.get('unknown', 0)} unknown, {c.get('ok', 0)} ok, {c.get('error', 0)} errors, {sum(n for k, n in c.items() if k.startswith('KNOWN'))} known wrong")

    print(f"{len(rows)} claims in {len(todo)} cases: {summary(count)}")
    if not a.no_mirror:
        print(f"  as written ({sum(by[False].values())}): {summary(by[False])}")
        print(f"  reversed ({sum(by[True].values())}): {summary(by[True])}")
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(rows, indent=1))
    sys.exit(1 if wrong else 0)


if __name__ == '__main__':
    main()

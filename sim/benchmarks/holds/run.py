#!/usr/bin/env python3
"""Holds: claims on parts held by fixes, each against its exact answer, so a
change to how holds are cut, traced or weighed is judged by how many it
answers rightly, wrongly, not at all, or refuses.

    python3 benchmarks/holds/run.py INK [-j N] [--only REGEX] [--json OUT]

Each case in cases/ is an ink model (`.case`: its claims hold only where ink
answers them rightly, so not a `.ink` the model rule requires to hold) with
claims; truths.json gives, per case and claim, the true value of the claim's
metric (or null: read only) and how it is known. Per claim the runner
reports one of:
  right        decided (holds or fails) as the truth says, and the truth
               inside the band: the stated grid error, and the read's error
               as its words state it (a peak stress "under it by X%", "over
               its read by Y%"; a displacement or temperature the solver's
               precision, taken as 2e-6 of the reading), or equal to the
               truth within that precision;
  wrong        a false yes or no, or the truth outside the band (widened
               by the truth's own accuracy, truths.json's `within`, of
               it: a CalculiX truth's last mesh step, solid.py);
  not settled  ink declines;
  refused      the study is refused, or no reading;
  no truth     read only.
and the totals. The families: hs- the statics hunt's (closed forms), hd- the
dynamics hunt's (buckling and modes: the solid's, by CalculiX, in place of
Euler's with shear and Timoshenko's frequencies, which miss it by up to
0.3%), rv- the holds reviews' (slabs, one-sided bars, corner blocks,
vanishing holds, tilted planes), en- the entries' own (a turned cantilever
against CalculiX, a held heat band), pl- a plate's buckling, kt- read only.

Recorded (one run each, 337 cases, 408 claims):
  main 2982806676:                right 236, wrong 69, not settled 96, refused 0, no truth 7
  holds branch (with the peak read
  at surface points and the
  lucky-grid rule's wide band):   right 178, wrong 11, not settled 168, refused 44, no truth 7
  main, its dynamics truths the
  solid's (CalculiX); hd-bsw-fp-r0.5
  read only (ink and CalculiX 0.02%
  apart, unknown):                  right 186, wrong 0, not settled 171, refused 43, no truth 8
  (per row: main.json; its 40 pivot and slide cases' rows read at
  ae7f021767, the rest at c6911fa42a, which no change since reaches)
One run: 51 min on the Spark (arm64, 20 cores), one case at a time at 4
threads; -j 3 here (32 threads) about a third of that.
"""
import argparse, concurrent.futures, json, os, pathlib, re, subprocess, sys, tempfile

HERE = pathlib.Path(__file__).resolve().parent
LIMIT = re.compile(r'\(\S+\s+\S+\s*(<=|>=|<|>)\s*([-+]?[0-9.]+(?:e[-+]?[0-9]+)?)(in|psi|degF|K)?\b')
TO_SI = {'in': lambda x: x * 25.4, 'psi': lambda x: x * 0.006894757293168, 'degF': lambda x: (x - 32) * 5 / 9, 'K': lambda x: x - 273.15}


def classify(r, truth, within=0.0):
    if r is None:
        return 'refused', None
    v, st = r.get('value'), r.get('status')
    if not isinstance(v, (int, float)):
        v = None
    if st not in ('holds', 'fails'):
        return ('refused' if r.get('unavailable') or v is None else 'not settled'), v
    if truth is None:
        return 'no truth', v
    m = LIMIT.search(r.get('claim') or '')
    if m:
        lim = TO_SI.get(m[3], lambda x: x)(float(m[2]))
        true = {'<': truth < lim, '<=': truth <= lim, '>': truth > lim, '>=': truth >= lim}[m[1]]
        if true != (st == 'holds'):
            return 'wrong', v
    if v is not None and abs(v - truth) <= 2e-6 * abs(truth):
        return 'right', v  # equal within the solver's precision
    g = r.get('grid_error')
    if v is not None and g is not None:
        below = above = g
        kind = (r.get('claim') or '').split()[0].lstrip('(')
        if kind in ('max-displacement', 'max-temperature'):
            below = above = g + 2e-6 * abs(v)
        words = r.get('estimate') or ''
        under, over = re.search(r'under it by ([0-9.]+)%', words), re.search(r'over its read by ([0-9.]+)%', words)
        if under:
            below += float(under[1]) / 100 * abs(v)
        if over:
            above += float(over[1]) / 100 * abs(v)
        # (and the truth's own accuracy: a beam formula's, `within`)
        below, above = below + within * abs(truth), above + within * abs(truth)
        if truth < v - below or truth > v + above:
            return 'wrong', v
    return 'right', v


def run(ink, case):
    with tempfile.TemporaryDirectory() as d:
        f = pathlib.Path(d) / (case.stem + '.ink')
        f.write_text(case.read_text())
        p = subprocess.run([ink, 'require', str(f)], capture_output=True, text=True, env={**os.environ, 'INK_NO_CACHE': '1'})
    try:
        return json.loads(p.stdout)
    except ValueError:
        return {'error': (p.stdout + p.stderr)[-300:]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('ink')
    ap.add_argument('-j', type=int, default=1)
    ap.add_argument('--only')
    ap.add_argument('--json')
    a = ap.parse_args()
    truths = json.loads((HERE / 'truths.json').read_text())
    cases = sorted(c for c in (HERE / 'cases').glob('*.case') if not a.only or re.search(a.only, c.stem))
    with concurrent.futures.ThreadPoolExecutor(a.j) as ex:
        out = dict(zip(cases, ex.map(lambda c: run(a.ink, c), cases)))
    tally, rows = {}, []
    for c, d in out.items():
        reqs = {r.get('name'): r for r in d.get('requirements', [])}
        for name, truth in truths[c.stem]['truths'].items():
            cls, v = classify(reqs.get(name), truth, truths[c.stem].get('within', 0.0))
            tally[cls] = tally.get(cls, 0) + 1
            rq = reqs.get(name) or {}
            rows.append({'case': c.stem, 'claim': name, 'class': cls, 'value': v, 'grid_error': rq.get('grid_error'), 'truth': truth, 'status': rq.get('status'), 'estimate': rq.get('estimate') or rq.get('undecided') or rq.get('unavailable')})
            if cls in ('wrong', 'refused'):
                print(f'{cls:12} {c.stem} {name}: {v} against {truth}')
    print(' '.join(f'{k} {tally.get(k, 0)}' for k in ('right', 'wrong', 'not settled', 'refused', 'no truth')))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    sys.exit(main())

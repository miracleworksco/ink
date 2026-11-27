#!/usr/bin/env python3
"""A solve's memory: the estimate it makes before it allocates (ink-fem's
memory.rs, `ink solve`'s grid.memory_estimate_mb) against the peak resident
set it reaches, per study and cell size.

    python3 benchmarks/memory/run.py INK [MODEL...] [--only REGEX] [--json OUT]

The models here at the cell sizes in SIZES, and every study in
benchmarks, examples and each MODEL given (a catalog's analyses,
say, INK_PATH set to find what they import) at their own. It prints each
estimate over the peak, and their range. Linux only.
"""
import argparse, json, os, pathlib, re, subprocess, sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SIZES = {'cube': [2, 1.4, 1], 'plate': [0.6, 0.4, 0.3], 'gyroid': [1, 0.7, 0.5], 'ring': [1, 0.6, 0.4], 'fin': [0.7, 0.8, 1, 1.4]}
OTHERS = ['benchmarks/*.ink', 'examples/*.ink']


def solve(ink, f, study, res):
    env = dict(os.environ, INK_NO_CACHE='1')
    env.setdefault('INK_PATH', str(ROOT))
    args = [ink, 'solve', str(f), '--study', study] + (['--res', str(res)] if res else [])
    p = subprocess.Popen(args, env=env, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    out = p.stdout.read()
    _, _, usage = os.wait4(p.pid, 0)
    try:
        d = json.loads(out)
    except json.JSONDecodeError:
        d = {}
    return d, usage.ru_maxrss * 1024


def studies(models):
    for f in sorted(HERE.glob('*.ink')):
        for res in SIZES[f.stem.split('-')[0]]:
            yield f, 's', res
    for f in [f for pat in OTHERS for f in sorted(ROOT.glob(pat))] + [pathlib.Path(m).resolve() for m in models]:
        for s in re.findall(r'^\(study ([\w-]+)', f.read_text(), re.M):
            yield f, s, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('ink')
    ap.add_argument('models', nargs='*', metavar='MODEL')
    ap.add_argument('--only')
    ap.add_argument('--json')
    a = ap.parse_args()
    rows = []
    for f, study, res in studies(a.models):
        name = f"{os.path.relpath(f, ROOT)}:{study}" + (f"@{res}" if res else '')
        if a.only and not re.search(a.only, name):
            continue
        d, peak = solve(a.ink, f, study, res)
        g = d.get('grid') or {}
        est = (g.get('memory_estimate_mb') or 0) * 2**20
        if not est:
            print(f"{name}: no estimate ({(d.get('error') or {}).get('message') or (d.get('solver') or {}).get('failure')})", flush=True)
            continue
        rows.append(dict(study=name, res=d.get('res'), full=g.get('cells_full', 0), cut=g.get('cells_cut', 0), estimate_gb=est / 2**30, peak_gb=peak / 2**30, ratio=est / peak))
        r = rows[-1]
        print(f"{name:58} {r['full'] + r['cut']:>9} cells  estimate {r['estimate_gb']:6.2f} GB  peak {r['peak_gb']:6.2f} GB  ratio {r['ratio']:.2f}", flush=True)
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(rows, indent=1))
    big = [r for r in rows if r['peak_gb'] >= 0.25]
    for label, rs in (('all', rows), ('over 0.25 GB', big)):
        if rs:
            print(f"{label}: ratio {min(r['ratio'] for r in rs):.2f} to {max(r['ratio'] for r in rs):.2f} over {len(rs)} solves")


if __name__ == '__main__':
    main()

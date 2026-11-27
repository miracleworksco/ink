#!/usr/bin/env python3
"""Cold, serial, bounded measurements. No solver controls are selected per case."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]  # the folder holding kernel/ and sim/ side by side


def metadata(path):
    lines = path.read_text().splitlines()
    rows = [json.loads(s.removeprefix('; scoreboard: ')) for s in lines if s.startswith('; scoreboard: ')]
    if len(rows) != 1:
        raise ValueError(f'{path}: need exactly one scoreboard header')
    code = '\n'.join(s.split(';', 1)[0] for s in lines)
    import re
    if re.search(r'\(\s*(res|support-margin|steps)\b', code):
        raise ValueError(f'{path}: default track forbids resolution/support-margin overrides')
    d = rows[0]
    if 'units_override' in d:
        d['units'] = d.pop('units_override')
    if d['target_kind'] not in ('exact', 'converged reference', 'idealisation'):
        raise ValueError(f'{path}: invalid target kind')
    if not math.isfinite(d['target']):
        raise ValueError(f'{path}: target must be finite')
    return d


def sha(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()


def code_sha(path):
    code = '\n'.join(line.split(';', 1)[0].strip() for line in path.read_text().splitlines())
    return hashlib.sha256(' '.join(code.split()).encode()).hexdigest()


def records(table):
    if isinstance(table, list):
        return table
    result = []
    for row in table['rows']:
        d = {}
        for name, value in zip(table['columns'], row):
            keys = name.split('.')
            p = d
            for key in keys[:-1]:
                p = p.setdefault(key, {})
            p[keys[-1]] = value
        result.append(d)
    return result


def tensor_quantity(t, meta):
    """Read the asked scalar from ink's symmetric Cartesian stress tensor."""
    a = [[t['xx'], t['xy'], t['zx']], [t['xy'], t['yy'], t['yz']], [t['zx'], t['yz'], t['zz']]]
    if 'tensor_direction' in meta:
        n = meta['tensor_direction']
        return sum(n[i] * a[i][j] * n[j] for i in range(3) for j in range(3))
    if meta.get('principal_max'):
        # Orthogonal Jacobi diagonalisation, entirely in the reported tensor.
        for _ in range(32):
            i,j = max(((0,1),(0,2),(1,2)), key=lambda ij: abs(a[ij[0]][ij[1]]))
            if abs(a[i][j]) <= 1e-14 * max(1.0, *(abs(a[k][k]) for k in range(3))):
                break
            angle = 0.5 * math.atan2(2*a[i][j], a[j][j]-a[i][i])
            c, sn = math.cos(angle), math.sin(angle)
            u,v,w = a[i][i],a[j][j],a[i][j]
            a[i][i] = c*c*u-2*c*sn*w+sn*sn*v
            a[j][j] = sn*sn*u+2*c*sn*w+c*c*v
            a[i][j] = a[j][i] = 0.0
            for k in range(3):
                if k not in (i,j):
                    u,v = a[k][i],a[k][j]
                    a[k][i] = a[i][k] = c*u-sn*v
                    a[k][j] = a[j][k] = sn*u+c*v
        return max(a[i][i] for i in range(3))
    raise ValueError('tensor probe needs a direction or principal_max')


def extract(d, meta):
    """Only a reported band, or the reported symmetric grid error; never an invented probe band."""
    if 'probe' in meta:
        if meta['probe'].startswith('reaction:'):
            name = meta['probe'].split(':', 1)[1]
            value = meta.get('probe_scale', 1) * next(r['total'] for r in d['reactions'] if r['name'] == name)
            return value, None, 'unknown', d.get('res'), 'ink supplies no accuracy band for a reaction'
        value = records(d['probes'])[meta.get('probe_index', 0)]
        for key in meta['probe'].split('.'):
            value = value[int(key)] if isinstance(value, list) else value[key]
        if isinstance(value, dict):
            value = tensor_quantity(value, meta)
        return value, None, 'unknown', d.get('res'), 'ink supplies no accuracy band for a point probe'
    row = next(x for x in d.get('requirements', []) if x['name'] == 'measurement')
    value = row.get('value')
    band = row.get('band')
    if band is not None:
        band = [band.get('low'), band.get('high')]
    elif value is not None and row.get('grid_error') is not None:
        e = row['grid_error']
        band = [value - e, value + e]
    reason = row.get('undecided') or row.get('unavailable') or row.get('estimate')
    return value, band, row['status'].split()[0], row.get('res'), reason


def reported_details(data, resolution):
    """Keep warnings and resolution identical in live runs and raw-output replay."""
    return dict(
        resolution_used=resolution if resolution is not None else (
            data.get('grid', {}).get('res_mm') or data.get('solver', {}).get('res_mm')),
        solver=data.get('solver'),
        warnings=data.get('warnings', []) + [
            w for row in data.get('requirements', []) for w in row.get('warnings', [])],
        ink_status=next((row.get('status') for row in data.get('requirements', [])
                         if row.get('name') == 'measurement'), None))


def score(value, band, meta, ink_status):
    target = meta['target']
    valid = value is not None and math.isfinite(value)
    finite_band = band is not None and len(band) == 2 and all(x is not None and math.isfinite(x) for x in band) and band[0] <= band[1]
    contains = band[0] <= target <= band[1] if finite_band else None
    width = band[1] - band[0] if finite_band else None
    scale = abs(target) if target else meta.get('engineering_scale')
    relative = width / scale if width is not None and scale else None
    qualified = meta['qualified']
    # Frozen policy: "narrow" is total width <=10% of the target/declared scale.
    confident_miss = bool(qualified and ink_status in ('holds', 'fails') and contains is False and relative is not None and relative <= 0.1)
    direction = ('low' if value < target else 'high') if valid and value != target else None
    unsafe_direction = bool(valid and ((meta['physics'] == 'static' and meta.get('unsafe_low', True) and abs(value) < abs(target)) or (meta['physics'] in ('modal', 'buckling') and value > target)))
    unsafe = confident_miss and unsafe_direction
    return dict(band_contains_target=contains, band_width=width, band_width_relative=relative,
                value_minus_target=value-target if valid else None,
                relative_error=(value-target)/scale if valid and scale else None,
                scored=qualified, confidently_wrong=confident_miss, unsafe=unsafe, unsafe_direction=unsafe_direction, direction=direction)


def command(argv, seconds, env, raw):
    cmd = ['systemd-run', '--user', '--scope', '-q', '-p', 'MemoryMax=16G', '-p', 'MemorySwapMax=0',
           'timeout', '--signal=TERM', '--kill-after=5s', str(max(0.1, seconds)), *argv]
    t = time.monotonic()
    with raw.with_suffix('.stdout').open('w') as out, raw.with_suffix('.stderr').open('w') as err:
        proc = subprocess.Popen(cmd, stdout=out, stderr=err, env=env, start_new_session=True)
        try:
            rc = proc.wait(timeout=seconds+10)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
            rc = 124
    try:
        data = json.loads(raw.with_suffix('.stdout').read_text())
    except (ValueError, OSError):
        data = None
    return rc, data, time.monotonic()-t, cmd


def run(path, args, binary_hash):
    meta = metadata(path)
    out = args.out / meta['id']
    out.mkdir(parents=True, exist_ok=True)
    signature = dict(binary_sha256=binary_hash, model_sha256=sha(path), timeout_s=args.timeout,
                     ink_commit=args.ink_commit)
    previous = out/'result.json'
    if args.reparse and not previous.exists():
        raise ValueError(f'{path}: no measurement to reparse')
    if args.reparse and previous.exists():
        snapshot = out/'model.ink'
        if not snapshot.exists() or code_sha(snapshot) != code_sha(path):
            raise ValueError(f'{path}: reparse requires identical solver inputs to recorded model')
        original = metadata(snapshot)
        for key in ('target', 'target_kind', 'qualified', 'units', 'points'):
            if original.get(key) != meta.get(key):
                raise ValueError(f'{path}: reparse cannot change {key}')
        prior = json.loads(previous.read_text())
        if prior.get('binary_sha256') != binary_hash or prior.get('ink_commit') != args.ink_commit:
            raise ValueError(f'{path}: reparse binary/build identity differs')
        if sha(snapshot) != prior.get('model_sha256'):
            raise ValueError(f'{path}: recorded model snapshot hash differs')
        prior['solver_input_sha256'] = code_sha(snapshot)
        prior['current_model_sha256'] = sha(path)
        raw = out/'answer.stdout'
        if raw.exists() and prior['status'] not in ('refused', 'timed out'):
            try:
                data = json.loads(raw.read_text())
                if data.get('ok') and not data.get('solver', {}).get('failure'):
                    v,b,st,res,reason = extract(data, meta)
                    prior.update(meta, value=v, band=b, status=st, resolution_used=res, reason=reason)
                    prior.update(reported_details(data, res))
                    prior.update(score(v,b,meta,st))
                    prior['postprocessed'] = 'Re-extracted existing raw output; no solver rerun or timing change'
            except (KeyError, StopIteration, TypeError, ValueError):
                pass
        previous.write_text(json.dumps(prior, indent=2, allow_nan=False)+'\n')
        return prior
    if args.resume and previous.exists():
        prior = json.loads(previous.read_text())
        if all(prior.get(k) == v for k,v in signature.items()):
            print(f"kept {meta['id']}: {prior['status']}", flush=True)
            return prior
    result = dict(meta, **signature, model=str(path.relative_to(REPO)), track='default',
                  status='unknown', value=None, band=None, resolution_used=None,
                  measured_at=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                  repository_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip())
    (out/'model.ink').write_bytes(path.read_bytes())
    result['solver_input_sha256'] = code_sha(path)
    start = time.monotonic()
    env = dict(os.environ, INK_PATH=str(REPO), INK_NO_CACHE='1')
    result['commands'] = []
    rc, check, _, cmd = command([str(args.ink), 'check', str(path), '--load'], args.timeout, env, out/'check')
    result['commands'].append(cmd)
    result['check_warnings'] = check.get('warnings', []) if check else []
    if rc == 124:
        result.update(status='timed out', reason='check exceeded the per-problem deadline')
    elif not check or not check.get('ok'):
        result.update(status='refused', reason=check.get('error') if check else f'check exited {rc}')
    else:
        remain = args.timeout - (time.monotonic()-start)
        if 'probe' in meta:
            argv = [str(args.ink), 'solve', str(path), '--study', meta['study'], '--points', json.dumps(meta['points'])]
        else:
            argv = [str(args.ink), 'require', str(path)]
        rc, data, _, cmd = command(argv, max(0.1, remain), env, out/'answer')
        result['commands'].append(cmd)
        if rc == 124:
            result.update(status='timed out', reason='per-problem deadline')
        elif not data or not data.get('ok'):
            result.update(status='refused', reason=data.get('error') if data else f'process exited {rc} (see stderr)')
        elif data.get('solver', {}).get('failure'):
            result.update(status='unknown', reason=data['solver']['failure'])
        else:
            try:
                v, b, st, res, reason = extract(data, meta)
                result.update(value=v, band=b, status=st, resolution_used=res, reason=reason)
                result.update(reported_details(data, res))
            except (KeyError, StopIteration, TypeError, ValueError) as e:
                result.update(status='unknown', reason=f'quantity unavailable: {type(e).__name__}: {e}')
    result['wall_time_s'] = time.monotonic()-start
    result.update(score(result['value'], result['band'], meta, result['status']))
    previous.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(f"{meta['id']}: {result['status']}, value={result['value']}, band={result['band']}, {result['wall_time_s']:.2f}s", flush=True)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ink', type=Path, default=REPO/'sim/target/fast/ink')
    ap.add_argument('--ink-commit', required=True, help='commit built into this binary (must have identical solver source)')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--case', default='', help='substring of filename; omitted runs all models')
    ap.add_argument('--timeout', type=float, default=600)
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--reparse', action='store_true', help='Re-extract saved output only, retaining original model hash and timing')
    args = ap.parse_args()
    args.ink = args.ink.expanduser().resolve()
    args.out = args.out.expanduser().resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    subprocess.run(['git', 'diff', '--exit-code', args.ink_commit, '--', 'kernel/crates', 'sim/ink-*', 'kernel/std', 'kernel/Cargo.toml', 'kernel/Cargo.lock', 'sim/Cargo.toml', 'sim/Cargo.lock'], cwd=REPO, check=True, stdout=subprocess.DEVNULL)
    paths = sorted(p for p in HERE.glob('*.ink') if args.case in p.name)
    if not paths:
        ap.error('no models selected')
    for path in paths:
        metadata(path)
    binary_hash = sha(args.ink)
    all_results = []
    for path in paths:
        all_results.append(run(path, args, binary_hash))
        (args.out/'results.json').write_text(json.dumps(all_results, indent=2, allow_nan=False)+'\n')
    # A qualified entry is one the scoreboard says can be scored: refused (it did not load or run), that is
    # false, and the run fails saying which. One that cannot be scored yet is marked unqualified, and why.
    refused = [r for r in all_results if r.get('qualified') and r['status'] == 'refused']
    for r in refused:
        print(f"FAILED {r['id']}: qualified, but refused: {r.get('reason')}", flush=True)
    sys.exit(1 if refused else 0)


if __name__ == '__main__':
    main()

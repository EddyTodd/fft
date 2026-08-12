#!/usr/bin/env python3
import argparse
import csv
import json
import os
import platform
import random
import subprocess
import time
from pathlib import Path


def cpu_model():
    try:
        for line in Path('/proc/cpuinfo').read_text().splitlines():
            if line.lower().startswith('model name'):
                return line.split(':', 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or 'unknown'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--binary', default='./build/fft-plan')
    ap.add_argument('--out', required=True)
    ap.add_argument('--sizes', default='64,256,1024,4096,16384,65536')
    ap.add_argument('--sessions', type=int, default=3)
    ap.add_argument('--samples', type=int, default=31)
    ap.add_argument('--warmups', type=int, default=5)
    ap.add_argument('--target-ms', type=float, default=2.0)
    ap.add_argument('--seed', type=int, default=20260812)
    ap.add_argument('--source-commit', default='unknown')
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    sizes = [int(x) for x in args.sizes.split(',')]
    rows = []
    rng = random.Random(args.seed)

    for session in range(args.sessions):
        order = sizes[:]
        rng.shuffle(order)
        for ordinal, n in enumerate(order):
            run_seed = args.seed ^ (session * 0x9E3779B1) ^ n
            cmd = [
                args.binary, '--benchmark', '--size', str(n), '--samples', str(args.samples),
                '--warmups', str(args.warmups), '--target-ms', str(args.target_ms),
                '--seed', str(run_seed), '--raw-csv'
            ]
            text = subprocess.check_output(cmd, text=True)
            for row in csv.DictReader(text.splitlines()):
                rows.append({'session': session, 'order': ordinal, 'seed': run_seed, **row})

    fields = ['session', 'order', 'seed', 'N', 'mode', 'sample', 'iterations_per_sample', 'ns']
    with (out / 'timings.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        'schema': 1,
        'timestamp_utc': time.strftime('%Y%m%dT%H%M%SZ', time.gmtime()),
        'source_commit': args.source_commit,
        'binary': args.binary,
        'platform': platform.platform(),
        'machine': platform.machine(),
        'cpu': cpu_model(),
        'python': platform.python_version(),
        'logical_cpus': os.cpu_count(),
        'sizes': sizes,
        'sessions': args.sessions,
        'samples_per_mode_per_session': args.samples,
        'warmups': args.warmups,
        'target_ms': args.target_ms,
        'randomization_seed': args.seed,
        'timing_semantics': 'forward+inverse pair divided by two; no allocation or trig setup in planned kernels; legacy uses the existing in-place radix-2 kernel; setup is timed separately; mode order randomized inside each sample and size order randomized per session'
    }
    (out / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'wrote {len(rows)} observations to {out}')


if __name__ == '__main__':
    main()

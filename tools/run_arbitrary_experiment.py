#!/usr/bin/env python3
import argparse
import csv
import gzip
import hashlib
import io
import json
import platform
import random
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def cpu_model():
    try:
        for line in Path('/proc/cpuinfo').read_text().splitlines():
            if line.lower().startswith('model name'):
                return line.split(':', 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or 'unknown'


def run(binary, arguments):
    return subprocess.run([str(binary), *arguments], check=True, text=True,
                          capture_output=True).stdout


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--sizes', default='17,31,61,127,257,509,1009,4093')
    parser.add_argument('--sessions', type=int, default=3)
    parser.add_argument('--samples', type=int, default=31)
    parser.add_argument('--setup-samples', type=int, default=1)
    parser.add_argument('--warmups', type=int, default=5)
    parser.add_argument('--target-ms', type=float, default=2.0)
    parser.add_argument('--seed', type=int, default=20260812)
    parser.add_argument('--source-commit', required=True)
    args = parser.parse_args()

    sizes = [int(value) for value in args.sizes.split(',') if value]
    args.out.mkdir(parents=True, exist_ok=True)
    raw_dir = args.out / 'raw'
    raw_dir.mkdir(exist_ok=True)

    info_lines = run(args.binary, ['--info']).strip().splitlines()
    info = dict(line.split('=', 1) for line in info_lines if '=' in line)
    if info.get('fftw_available') != 'yes':
        raise SystemExit('FFTW runtime is required for the formal arbitrary-plan experiment')

    fieldnames = ['session', 'size_order', 'phase', 'backend', 'algorithm', 'N', 'sample',
                  'mode_order', 'iterations', 'ns_per_transform', 'convolution_size', 'direct_cyclic']
    raw_files = []
    total_rows = 0

    for session in range(args.sessions):
        rng = random.Random(args.seed + session)
        ordered_sizes = sizes[:]
        rng.shuffle(ordered_sizes)
        rows = []

        for size_order, n in enumerate(ordered_sizes):
            stdout = run(args.binary, [
                '--raw-csv', '--size', str(n), '--samples', str(args.samples),
                '--setup-samples', str(args.setup_samples), '--warmups', str(args.warmups),
                '--target-ms', str(args.target_ms),
                '--seed', str(args.seed ^ (session << 32) ^ n),
            ])
            for row in csv.DictReader(stdout.splitlines()):
                rows.append({'session': session, 'size_order': size_order, **row})

        path = raw_dir / f'timings-session{session}.csv.gz'
        with open(path, 'wb') as raw_stream:
            with gzip.GzipFile(filename='', mode='wb', fileobj=raw_stream, mtime=0) as compressed:
                text = io.TextIOWrapper(compressed, encoding='utf-8', newline='')
                writer = csv.DictWriter(text, fieldnames=fieldnames, lineterminator='\n')
                writer.writeheader()
                writer.writerows(rows)
                text.flush()
                text.detach()

        raw_files.append(path)
        total_rows += len(rows)
        print(f'{path}: {len(rows)} rows')

    expected = args.sessions * len(sizes) * (6 * args.samples + 4 * args.setup_samples)
    if total_rows != expected:
        raise SystemExit(f'raw row cardinality mismatch: got {total_rows}, expected {expected}')

    metadata = {
        'schema_version': 1,
        'experiment': 'arbitrary-plan-v1',
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'source_commit': args.source_commit,
        'source_provenance': 'Formal run generated from the source tree matching source_commit; later commits add evidence and documentation only.',
        'binary': str(args.binary),
        'binary_sha256': sha256(args.binary),
        'platform': platform.platform(),
        'machine': platform.machine(),
        'cpu_model': cpu_model(),
        'python': platform.python_version(),
        'fftw_library': info.get('fftw_library', ''),
        'fftw_version': info.get('fftw_version', ''),
        'sizes': sizes,
        'sessions': args.sessions,
        'samples_per_execution_mode_session': args.samples,
        'setup_samples_per_mode_session': args.setup_samples,
        'warmups': args.warmups,
        'target_ms': args.target_ms,
        'seed': args.seed,
        'raw_observations': total_rows,
        'timing_semantics': 'Prime-length complex transforms. Planned modes use persistent precomputed convolution kernels and caller-owned output/scratch buffers; legacy modes include per-call allocation/setup. FFTW uses persistent plans and aligned buffers. All execution samples time forward+inverse pairs and divide by two; inverse normalization is included.',
        'raw_files': [
            {'path': str(path.relative_to(args.out)), 'sha256': sha256(path), 'bytes': path.stat().st_size}
            for path in raw_files
        ],
    }
    (args.out / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'total rows: {total_rows}')


if __name__ == '__main__':
    main()

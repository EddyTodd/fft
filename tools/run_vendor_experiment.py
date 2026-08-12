#!/usr/bin/env python3
import argparse, csv, gzip, hashlib, json, os, platform, random, subprocess, time
from pathlib import Path


def parse_info(binary):
    p = subprocess.run([binary, '--info'], text=True, capture_output=True, check=False)
    if p.returncode == 77:
        raise SystemExit('FFTW runtime unavailable')
    if p.returncode:
        raise SystemExit(p.stderr or p.stdout)
    out = {}
    for line in p.stdout.splitlines():
        if ': ' in line:
            k, v = line.split(': ', 1)
            out[k] = v
    return out


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
    ap.add_argument('--binary', default='build/fft-vendor')
    ap.add_argument('--sizes', default='64,256,1024,4096,16384,65536')
    ap.add_argument('--sessions', type=int, default=3)
    ap.add_argument('--samples', type=int, default=31)
    ap.add_argument('--setup-samples', type=int, default=1)
    ap.add_argument('--warmups', type=int, default=5)
    ap.add_argument('--target-ms', type=float, default=2.0)
    ap.add_argument('--seed', type=int, default=20260812)
    ap.add_argument('--source-commit', default='unknown')
    ap.add_argument('--build-description', default='unspecified')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    sizes = [int(x) for x in args.sizes.split(',') if x]
    out = Path(args.output)
    raw = out / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    info = parse_info(args.binary)
    metadata = {
        'schema': 1,
        'timestamp_utc': time.strftime('%Y%m%dT%H%M%SZ', time.gmtime()),
        'source_commit': args.source_commit,
        'build_description': args.build_description,
        'binary': args.binary,
        'platform': platform.platform(),
        'machine': platform.machine(),
        'cpu': cpu_model(),
        'logical_cpus': os.cpu_count(),
        'fftw_library': info.get('library', ''),
        'fftw_version': info.get('version', ''),
        'sizes': sizes,
        'sessions': args.sessions,
        'execution_samples_per_mode_per_session': args.samples,
        'cold_setup_samples_per_mode_per_session': args.setup_samples,
        'warmups': args.warmups,
        'target_ms': args.target_ms,
        'randomization_seed': args.seed,
        'semantics': 'forward+normalized-inverse pair divided by two; FFTW inverse scaling included; caller buffers allocated outside setup/execution timing; FFTW cold setup forgets wisdom before each forward+inverse plan pair',
        'raw_files': [],
    }
    (out / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')

    for session in range(args.sessions):
        rng = random.Random(args.seed ^ (session * 0x9E3779B1))
        order = sizes[:]
        rng.shuffle(order)
        path = raw / f'timings-session{session}.csv.gz'
        with gzip.open(path, 'wt', newline='') as f:
            writer = None
            ordinal = 0
            for n in order:
                cmd = [args.binary, '--raw-csv', '--size', str(n), '--samples', str(args.samples),
                       '--setup-samples', str(args.setup_samples), '--warmups', str(args.warmups),
                       '--target-ms', str(args.target_ms), '--seed', str(args.seed ^ session ^ n)]
                p = subprocess.run(cmd, text=True, capture_output=True, check=True)
                rows = list(csv.DictReader(p.stdout.splitlines()))
                for row in rows:
                    enriched = {'session': session, 'size_order': ordinal, **row}
                    if writer is None:
                        writer = csv.DictWriter(f, fieldnames=list(enriched))
                        writer.writeheader()
                    writer.writerow(enriched)
                ordinal += 1
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        metadata['raw_files'].append({'name': path.name, 'sha256': digest})
        (out / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
        print(path, flush=True)


if __name__ == '__main__':
    main()

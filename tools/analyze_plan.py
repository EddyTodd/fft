#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


def percentile(values, p):
    values = sorted(values)
    pos = p * (len(values) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    return values[lo] if lo == hi else values[lo] * (hi - pos) + values[hi] * (pos - lo)


def ratio_ci(a, b, rng, reps=5000):
    boots = []
    na, nb = len(a), len(b)
    for _ in range(reps):
        ma = statistics.median(a[rng.randrange(na)] for __ in range(na))
        mb = statistics.median(b[rng.randrange(nb)] for __ in range(nb))
        boots.append(ma / mb)
    return percentile(boots, .025), percentile(boots, .975)


def input_paths(source):
    source = Path(source)
    if source.is_dir():
        return sorted([*source.glob('*.csv'), *source.glob('*.csv.gz')])
    return [source]


def load_rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='') as f:
        return list(csv.DictReader(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source', help='raw CSV, .csv.gz, or directory containing either')
    ap.add_argument('--seed', type=int, default=20260812)
    ap.add_argument('--output')
    args = ap.parse_args()

    paths = input_paths(args.source)
    if not paths:
        raise SystemExit('no CSV inputs found')
    rows = []
    for path in paths:
        rows.extend(load_rows(path))

    grouped = defaultdict(list)
    for row in rows:
        grouped[(int(row['N']), row['mode'])].append(float(row['ns']))

    required = {'complex-setup', 'real-setup', 'legacy-complex', 'planned-complex', 'planned-real'}
    sizes = sorted({key[0] for key in grouped})
    for n in sizes:
        missing = required - {mode for size, mode in grouped if size == n}
        if missing:
            raise SystemExit(f'N={n} missing modes: {sorted(missing)}')

    rng = random.Random(args.seed)
    lines = [
        '# Planned and real FFT analysis',
        '',
        f'Raw observations: **{len(rows):,}** from **{len(paths)}** input file(s). Bootstrap seed: `{args.seed}`.',
        '',
        '| N | Legacy complex | Planned complex | Plan speedup (95% CI) | Planned real | Real speedup (95% CI) | Complex setup | Plan break-even | Real setup | Real break-even |',
        '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|'
    ]

    for n in sizes:
        legacy = grouped[n, 'legacy-complex']
        planned = grouped[n, 'planned-complex']
        real = grouped[n, 'planned-real']
        complex_setup = grouped[n, 'complex-setup']
        real_setup = grouped[n, 'real-setup']
        legacy_m = statistics.median(legacy)
        planned_m = statistics.median(planned)
        real_m = statistics.median(real)
        complex_setup_m = statistics.median(complex_setup)
        real_setup_m = statistics.median(real_setup)
        plan_ci = ratio_ci(legacy, planned, rng)
        real_ci = ratio_ci(planned, real, rng)
        plan_break_even = complex_setup_m / (legacy_m - planned_m) if legacy_m > planned_m else math.inf
        real_break_even = max(0.0, real_setup_m - complex_setup_m) / (planned_m - real_m) if planned_m > real_m else math.inf

        def fmt(x):
            return f'{x / 1000:.2f} µs' if x >= 1000 else f'{x:.1f} ns'

        lines.append(
            f'| {n} | {fmt(legacy_m)} | {fmt(planned_m)} | **{legacy_m / planned_m:.3f}×** [{plan_ci[0]:.3f}, {plan_ci[1]:.3f}] | '
            f'{fmt(real_m)} | **{planned_m / real_m:.3f}×** [{real_ci[0]:.3f}, {real_ci[1]:.3f}] | {fmt(complex_setup_m)} | '
            f'{plan_break_even:.2f} transforms | {fmt(real_setup_m)} | {real_break_even:.2f} transforms |'
        )

    lines += [
        '',
        '## Interpretation',
        '',
        '- `legacy-complex` and `planned-complex` execute the same radix-2 decomposition; the planned path moves bit-reversal and trigonometric twiddle construction out of the timed kernel.',
        '- `planned-real` exploits real-input Hermitian symmetry by packing even/odd samples into an N/2 complex FFT and reconstructing only N/2+1 unique frequency bins.',
        '- Setup break-even is descriptive for this implementation and machine: plan construction cost divided by the median per-transform savings.',
        '- The real-path break-even compares the extra real-plan setup cost against planned-complex setup; a value below one means the specialized real representation repays any extra setup within the first transform.',
        '- Confidence intervals use independent nonparametric bootstrap resampling of the pooled raw observations. They characterize this recorded environment, not universal hardware rankings.',
        ''
    ]

    text = '\n'.join(lines)
    if args.output:
        Path(args.output).write_text(text + '\n')
    else:
        print(text)


if __name__ == '__main__':
    main()

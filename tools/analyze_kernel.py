#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path


def percentile(values, p):
    values = sorted(values)
    pos = p * (len(values) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    return values[lo] if lo == hi else values[lo] * (hi - pos) + values[hi] * (pos - lo)


def ratio_ci(numerator, denominator, rng, reps):
    bootstraps = []
    for _ in range(reps):
        a = statistics.median(numerator[rng.randrange(len(numerator))] for __ in range(len(numerator)))
        b = statistics.median(denominator[rng.randrange(len(denominator))] for __ in range(len(denominator)))
        bootstraps.append(a / b)
    return percentile(bootstraps, .025), percentile(bootstraps, .975)


def input_files(path):
    root = Path(path)
    if root.is_dir():
        return sorted(list(root.rglob('*.csv')) + list(root.rglob('*.csv.gz')))
    return [root]


def read_rows(path):
    rows = []
    for file in input_files(path):
        opener = gzip.open if file.suffix == '.gz' else open
        with opener(file, 'rt', newline='') as stream:
            rows.extend(csv.DictReader(stream))
    return rows


def logical_policy(row):
    if row['backend'] == 'kernel' and row['policy'].startswith('auto->'):
        return 'auto'
    return row['policy']


def key(row):
    return row['backend'], logical_policy(row)


def format_time(ns):
    if ns >= 1e6:
        return f'{ns / 1e6:.3f} ms'
    if ns >= 1e3:
        return f'{ns / 1e3:.2f} us'
    return f'{ns:.1f} ns'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('--seed', type=int, default=20260812)
    parser.add_argument('--bootstrap', type=int, default=5000)
    parser.add_argument('--output')
    args = parser.parse_args()

    rows = read_rows(args.input)
    execution = defaultdict(list)
    setup = defaultdict(list)
    auto_sessions = defaultdict(dict)

    for row in rows:
        n = int(row['N'])
        value = float(row['ns_per_transform'])
        destination = execution if row['phase'] == 'execution' else setup
        destination[(n, *key(row))].append(value)
        if row['backend'] == 'kernel' and row['policy'].startswith('auto->'):
            auto_sessions[n][int(row['session'])] = row['policy'].split('->', 1)[1]

    rng = random.Random(args.seed)
    sizes = sorted({group[0] for group in execution})
    labels = {
        'legacy': ('fftlab-plan', 'legacy'),
        'scalar': ('kernel', 'scalar'),
        'avx2': ('kernel', 'avx2'),
        'avx512': ('kernel', 'avx512'),
        'auto': ('kernel', 'auto'),
        'estimate': ('fftw', 'estimate'),
        'measure': ('fftw', 'measure'),
    }

    lines = [
        '# SIMD kernel analysis',
        '',
        f'Raw observations: **{len(rows):,}**. Bootstrap seed: `{args.seed}`; repetitions: **{args.bootstrap:,}**.',
        '',
        '## Steady-state execution',
        '',
        '| N | v3 plan | Scalar codelet | AVX2/FMA | AVX-512/FMA | Auto | FFTW ESTIMATE | FFTW MEASURE | Best SIMD / v3 | Gap to MEASURE closed |',
        '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]

    detail = []
    setup_lines = []
    for n in sizes:
        medians = {}
        for name, backend_policy in labels.items():
            values = execution[(n, *backend_policy)]
            if not values:
                raise SystemExit(f'missing execution mode {name} at N={n}')
            medians[name] = statistics.median(values)

        best_name = 'avx2' if medians['avx2'] <= medians['avx512'] else 'avx512'
        best = medians[best_name]
        legacy_speedup = medians['legacy'] / best
        denominator = medians['legacy'] - medians['measure']
        closure = (medians['legacy'] - best) / denominator if denominator > 0 else math.nan
        lines.append(
            f"| {n} | {format_time(medians['legacy'])} | {format_time(medians['scalar'])} | "
            f"{format_time(medians['avx2'])} | {format_time(medians['avx512'])} | "
            f"{format_time(medians['auto'])} | {format_time(medians['estimate'])} | "
            f"{format_time(medians['measure'])} | **{legacy_speedup:.3f}x** ({best_name}) | "
            f"**{100 * closure:.1f}%** |"
        )

        avx_ci = ratio_ci(execution[(n, 'kernel', 'avx2')], execution[(n, 'kernel', 'avx512')],
                          rng, args.bootstrap)
        avx_ratio = medians['avx2'] / medians['avx512']
        if avx_ci[0] > 1:
            verdict = 'AVX-512 faster'
        elif avx_ci[1] < 1:
            verdict = 'AVX2 faster'
        else:
            verdict = 'unresolved'

        best_ci = ratio_ci(execution[(n, 'fftlab-plan', 'legacy')],
                           execution[(n, 'kernel', best_name)], rng, args.bootstrap)
        fftw_gap = best / medians['measure']
        selections = Counter(auto_sessions[n].values())
        selection_text = ', '.join(
            f'{name}: {count}/{len(auto_sessions[n])}' for name, count in sorted(selections.items()))
        detail.append(
            f'- **N={n}:** AVX2/AVX-512 median ratio {avx_ratio:.3f}x, 95% CI '
            f'[{avx_ci[0]:.3f}, {avx_ci[1]:.3f}] -> {verdict}. Best SIMD speedup over v3 plan '
            f'{legacy_speedup:.3f}x, 95% CI [{best_ci[0]:.3f}, {best_ci[1]:.3f}]. Best SIMD remains '
            f'{fftw_gap:.2f}x slower than FFTW MEASURE. Auto selections: {selection_text}.'
        )

        setup_medians = {
            name: statistics.median(setup[(n, *backend_policy)])
            for name, backend_policy in labels.items()
        }
        explicit_setup = setup_medians[best_name]
        extra_setup = max(0.0, explicit_setup - setup_medians['legacy'])
        saved = medians['legacy'] - best
        break_even = extra_setup / saved if saved > 0 else math.inf
        auto_extra = max(0.0, setup_medians['auto'] - explicit_setup)
        setup_lines.append(
            f"| {n} | {format_time(setup_medians['legacy'])} | {format_time(setup_medians[best_name])} "
            f"({best_name}) | {break_even:.2f} transforms | {format_time(setup_medians['auto'])} | "
            f"{format_time(auto_extra)} | {format_time(setup_medians['estimate'])} | "
            f"{format_time(setup_medians['measure'])} |"
        )

    lines += [
        '',
        '### AVX2 vs AVX-512 and auto-selection',
        '',
        *detail,
        '',
        '## Planning cost',
        '',
        '| N | v3 setup | Best explicit SIMD setup | SIMD break-even vs v3 | Auto-tuned setup | Auto premium vs best explicit | FFTW ESTIMATE setup | FFTW MEASURE setup |',
        '|---:|---:|---:|---:|---:|---:|---:|---:|',
        *setup_lines,
        '',
        '## Interpretation',
        '',
        '- The scalar codelet, AVX2, and AVX-512 paths share the same swap-list permutation and stage-contiguous twiddle layout. Their differences isolate instruction-width/code-generation effects more cleanly than comparing unrelated FFT algorithms.',
        '- `KernelRadix2Plan` stores approximately N-1 stage-local complex twiddles rather than N/2 globally indexed twiddles. The additional persistent memory is part of the setup tradeoff and is not hidden from the analysis.',
        '- The auto policy performs five rotated timing rounds over every supported candidate during plan construction. Auto tuning is therefore a portability/planning feature, not free execution speed; its full construction cost is reported separately.',
        '- AVX2 and AVX-512 are explicit research modes even when auto chooses one of them. Wider vectors are only called faster when the bootstrap speedup interval excludes parity.',
        '- FFTW remains a separately engineered adaptive library. Closing part of the latency gap with one vectorized radix-2 codelet does not imply architectural equivalence or a universal ranking.',
        '- All findings are scoped to the recorded virtualized environment and compiler/runtime configuration.',
        '',
    ]

    text = '\n'.join(lines)
    if args.output:
        Path(args.output).write_text(text + '\n')
    else:
        print(text)


if __name__ == '__main__':
    main()

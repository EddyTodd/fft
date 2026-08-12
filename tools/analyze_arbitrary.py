#!/usr/bin/env python3
import argparse
import base64
import csv
import gzip
import io
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
    root = Path(path)
    if root.is_dir():
        parts = sorted(root.rglob('*.b64part*'))
        groups = {}
        for part in parts:
            prefix = str(part).split('.b64part', 1)[0]
            groups.setdefault(prefix, []).append(part)
        for prefix, group in sorted(groups.items()):
            encoded = ''.join(part.read_text().strip() for part in sorted(group))
            decoded = base64.b64decode(encoded)
            text = gzip.decompress(decoded).decode('utf-8') if prefix.endswith('.gz') else decoded.decode('utf-8')
            rows.extend(csv.DictReader(io.StringIO(text)))
        regular = [file for file in input_files(root) if '.b64part' not in file.name]
    else:
        regular = input_files(root)
    for file in regular:
        opener = gzip.open if file.suffix == '.gz' else open
        with opener(file, 'rt', newline='') as stream:
            rows.extend(csv.DictReader(stream))
    return rows


def key(row):
    return row['backend'], row['algorithm']


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
    structure = {}
    for row in rows:
        n = int(row['N'])
        value = float(row['ns_per_transform'])
        destination = execution if row['phase'] == 'execution' else setup
        destination[(n, *key(row))].append(value)
        if row['backend'] == 'fftlab' and row['algorithm'].startswith('planned-'):
            structure[(n, row['algorithm'])] = (int(row['convolution_size']), row['direct_cyclic'] == '1')

    rng = random.Random(args.seed)
    sizes = sorted({group[0] for group in execution})
    labels = {
        'legacy-blue': ('fftlab', 'legacy-bluestein'),
        'blue': ('fftlab', 'planned-bluestein'),
        'legacy-rader': ('fftlab', 'legacy-rader'),
        'rader': ('fftlab', 'planned-rader'),
        'estimate': ('fftw', 'estimate'),
        'measure': ('fftw', 'measure'),
    }

    lines = [
        '# Arbitrary-length plan analysis',
        '',
        f'Raw observations: **{len(rows):,}**. Bootstrap seed: `{args.seed}`; repetitions: **{args.bootstrap:,}**.',
        '',
        '## Prime-length execution',
        '',
        '| N | Bluestein M | Rader M | Legacy Bluestein | Planned Bluestein | Legacy Rader | Planned Rader | Planned winner | FFTW ESTIMATE | FFTW MEASURE |',
        '|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|',
    ]
    detail = []
    setup_lines = []

    for n in sizes:
        med = {}
        for name, bp in labels.items():
            values = execution[(n, *bp)]
            if not values:
                raise SystemExit(f'missing execution mode {name} at N={n}')
            med[name] = statistics.median(values)

        blue_m, _ = structure[(n, 'planned-bluestein')]
        rader_m, direct = structure[(n, 'planned-rader')]
        ratio = med['blue'] / med['rader']
        ci = ratio_ci(execution[(n, 'fftlab', 'planned-bluestein')],
                      execution[(n, 'fftlab', 'planned-rader')], rng, args.bootstrap)
        if ci[0] > 1:
            winner = 'Rader'
            best = med['rader']
        elif ci[1] < 1:
            winner = 'Bluestein'
            best = med['blue']
        else:
            winner = 'unresolved'
            best = min(med['blue'], med['rader'])
        lines.append(
            f"| {n} | {blue_m} | {rader_m}{' (cyclic)' if direct else ''} | {format_time(med['legacy-blue'])} | "
            f"{format_time(med['blue'])} | {format_time(med['legacy-rader'])} | {format_time(med['rader'])} | "
            f"**{winner}** | {format_time(med['estimate'])} | {format_time(med['measure'])} |"
        )

        blue_plan_speed = med['legacy-blue'] / med['blue']
        rader_plan_speed = med['legacy-rader'] / med['rader']
        fftw_gap = best / med['measure']
        detail.append(
            f"- **N={n}:** planned Bluestein / planned Rader median ratio {ratio:.3f}x, 95% CI "
            f"[{ci[0]:.3f}, {ci[1]:.3f}] -> {winner}. Planning improves Bluestein {blue_plan_speed:.2f}x "
            f"and Rader {rader_plan_speed:.2f}x over the legacy setup-inclusive APIs. Best planned fftlab remains "
            f"{fftw_gap:.2f}x slower than FFTW MEASURE."
        )

        sb = statistics.median(setup[(n, 'fftlab', 'planned-bluestein')])
        sr = statistics.median(setup[(n, 'fftlab', 'planned-rader')])
        se = statistics.median(setup[(n, 'fftw', 'estimate')])
        sm = statistics.median(setup[(n, 'fftw', 'measure')])
        blue_saved = med['legacy-blue'] - med['blue']
        rader_saved = med['legacy-rader'] - med['rader']
        blue_be = sb / blue_saved if blue_saved > 0 else math.inf
        rader_be = sr / rader_saved if rader_saved > 0 else math.inf
        setup_lines.append(
            f"| {n} | {format_time(sb)} | {blue_be:.2f} | {format_time(sr)} | {rader_be:.2f} | "
            f"{format_time(se)} | {format_time(sm)} |"
        )

    lines += [
        '',
        '### Pairwise interpretation',
        '',
        *detail,
        '',
        '## Planning cost and amortization',
        '',
        '| N | Bluestein setup | Bluestein break-even vs legacy | Rader setup | Rader break-even vs legacy | FFTW ESTIMATE setup | FFTW MEASURE setup |',
        '|---:|---:|---:|---:|---:|---:|---:|',
        *setup_lines,
        '',
        '## Interpretation',
        '',
        '- Bluestein and Rader are both reduced to power-of-two convolution kernels so their planned comparison isolates reduction structure and persistent precomputation more cleanly than the legacy allocation-heavy APIs.',
        '- When `N-1` is a power of two, planned Rader performs the cyclic convolution directly at length `N-1`; otherwise it zero-pads a linear convolution and folds it modulo `N-1`. That structural distinction is reported explicitly in the table.',
        '- A lower convolution length is expected to help, but equal convolution lengths do not imply equal machine performance because Rader adds permutations/folding while Bluestein adds chirp multiplications.',
        '- Setup and execution are separate response variables. The legacy functions intentionally remain setup-inclusive historical API baselines; planned-vs-planned and planned-vs-FFTW comparisons use persistent plans.',
        '- FFTW is an adaptive production library and remains a reference point, not a claim that the local reductions have equivalent engineering maturity.',
        '- Results are scoped to the recorded compiler/runtime/virtualized host and should not be converted into a universal prime-size dispatcher without physical-hardware replication.',
        '',
    ]

    text = '\n'.join(lines)
    if args.output:
        Path(args.output).write_text(text + '\n')
    else:
        print(text)


if __name__ == '__main__':
    main()

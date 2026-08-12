#!/usr/bin/env python3
import argparse, csv, gzip, math, random, statistics
from collections import defaultdict
from pathlib import Path


def open_text(path):
    return gzip.open(path, 'rt', newline='') if str(path).endswith('.gz') else open(path, newline='')


def input_files(path):
    p = Path(path)
    if p.is_dir():
        return sorted(list(p.glob('*.csv')) + list(p.glob('*.csv.gz')))
    return [p]


def percentile(values, p):
    values = sorted(values)
    pos = p * (len(values) - 1)
    lo, hi = math.floor(pos), math.ceil(pos)
    return values[lo] if lo == hi else values[lo] * (hi - pos) + values[hi] * (pos - lo)


def ratio_ci(a, b, rng, reps=5000):
    boots = []
    for _ in range(reps):
        ma = statistics.median(a[rng.randrange(len(a))] for __ in range(len(a)))
        mb = statistics.median(b[rng.randrange(len(b))] for __ in range(len(b)))
        boots.append(ma / mb)
    return percentile(boots, .025), percentile(boots, .975)


def common_language_faster(a, b):
    wins = ties = 0
    for x in a:
        for y in b:
            wins += x < y
            ties += x == y
    return (wins + .5 * ties) / (len(a) * len(b))


def fmt(ns):
    if ns >= 1e6:
        return f'{ns / 1e6:.2f} ms'
    if ns >= 1e3:
        return f'{ns / 1e3:.2f} µs'
    return f'{ns:.1f} ns'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input')
    ap.add_argument('--seed', type=int, default=20260812)
    ap.add_argument('--output')
    args = ap.parse_args()

    rows = []
    for path in input_files(args.input):
        with open_text(path) as f:
            rows.extend(csv.DictReader(f))
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row['phase'], row['backend'], row['planner'], row['kind'], int(row['N']))].append(float(row['ns_per_operation']))

    rng = random.Random(args.seed)
    sizes = sorted({key[-1] for key in grouped})
    lines = ['# FFTW comparison analysis', '', f'Raw observations: **{len(rows):,}**. Bootstrap seed: `{args.seed}`.', '']

    for kind in ['complex', 'real']:
        lines += [
            f'## {kind.title()} execution', '',
            '| N | fftlab planned | FFTW ESTIMATE | ESTIMATE speedup (95% CI) | FFTW MEASURE | MEASURE speedup (95% CI) | MEASURE / ESTIMATE |',
            '|---:|---:|---:|---:|---:|---:|---:|'
        ]
        for n in sizes:
            own = grouped['execution', 'fftlab', 'native', kind, n]
            estimate = grouped['execution', 'fftw', 'estimate', kind, n]
            measure = grouped['execution', 'fftw', 'measure', kind, n]
            own_m, estimate_m, measure_m = map(statistics.median, [own, estimate, measure])
            ci_est = ratio_ci(own, estimate, rng)
            ci_measure = ratio_ci(own, measure, rng)
            ci_rigor = ratio_ci(estimate, measure, rng)
            lines.append(
                f'| {n} | {fmt(own_m)} | {fmt(estimate_m)} | **{own_m / estimate_m:.2f}×** [{ci_est[0]:.2f}, {ci_est[1]:.2f}] | '
                f'{fmt(measure_m)} | **{own_m / measure_m:.2f}×** [{ci_measure[0]:.2f}, {ci_measure[1]:.2f}] | '
                f'**{estimate_m / measure_m:.3f}×** [{ci_rigor[0]:.3f}, {ci_rigor[1]:.3f}] |'
            )
        lines.append('')

    lines += [
        '## Planning economics', '',
        '| Kind | N | fftlab setup | FFTW ESTIMATE setup | ESTIMATE break-even vs fftlab | FFTW MEASURE cold setup | MEASURE break-even vs ESTIMATE |',
        '|---|---:|---:|---:|---:|---:|---:|'
    ]
    for kind in ['complex', 'real']:
        for n in sizes:
            own_setup = statistics.median(grouped['setup', 'fftlab', 'native', kind, n])
            estimate_setup = statistics.median(grouped['setup', 'fftw', 'estimate', kind, n])
            measure_setup = statistics.median(grouped['setup', 'fftw', 'measure', kind, n])
            own_exec = statistics.median(grouped['execution', 'fftlab', 'native', kind, n])
            estimate_exec = statistics.median(grouped['execution', 'fftw', 'estimate', kind, n])
            measure_exec = statistics.median(grouped['execution', 'fftw', 'measure', kind, n])
            estimate_break_even = (estimate_setup - own_setup) / (own_exec - estimate_exec) if own_exec > estimate_exec and estimate_setup > own_setup else 0.0
            measure_break_even = (measure_setup - estimate_setup) / (estimate_exec - measure_exec) if estimate_exec > measure_exec and measure_setup > estimate_setup else 0.0
            estimate_text = 'immediate' if own_exec > estimate_exec and estimate_setup <= own_setup else (f'{estimate_break_even:,.1f} transforms' if own_exec > estimate_exec else 'never')
            measure_text = 'immediate' if estimate_exec > measure_exec and measure_setup <= estimate_setup else (f'{measure_break_even:,.0f} transforms' if estimate_exec > measure_exec else 'never')
            lines.append(f'| {kind} | {n} | {fmt(own_setup)} | {fmt(estimate_setup)} | {estimate_text} | {fmt(measure_setup)} | {measure_text} |')

    lines += [
        '', '## Amortized winner', '',
        'Winner minimizes `setup + K × execution` under the normalized forward/inverse contract.', '',
        '| Kind | N | K=1 | K=10 | K=100 | K=10,000 | K=1,000,000 |',
        '|---|---:|---|---|---|---|---|'
    ]
    strategies = [('fftlab', 'native', 'fftlab'), ('fftw', 'estimate', 'FFTW_ESTIMATE'), ('fftw', 'measure', 'FFTW_MEASURE')]
    for kind in ['complex', 'real']:
        for n in sizes:
            winners = []
            for count in [1, 10, 100, 10000, 1000000]:
                costs = []
                for backend, planner, label in strategies:
                    setup = statistics.median(grouped['setup', backend, planner, kind, n])
                    execution = statistics.median(grouped['execution', backend, planner, kind, n])
                    costs.append((setup + count * execution, label))
                winners.append(min(costs)[1])
            lines.append(f'| {kind} | {n} | ' + ' | '.join(winners) + ' |')

    lines += [
        '', '## Effect-size checks', '',
        'Common-language values below are the probability that a random FFTW MEASURE execution observation is faster than a random fftlab planned observation.', '',
        '| Kind | N | P(FFTW MEASURE faster) |', '|---|---:|---:|'
    ]
    for kind in ['complex', 'real']:
        for n in sizes:
            probability = common_language_faster(
                grouped['execution', 'fftw', 'measure', kind, n],
                grouped['execution', 'fftlab', 'native', kind, n])
            lines.append(f'| {kind} | {n} | {probability:.4f} |')

    lines += [
        '', '## Interpretation', '',
        '- FFTW execution is normalized to the same mathematical contract as fftlab: forward plus normalized inverse. FFTW’s required `1/N` inverse scaling is included in its timed path.',
        '- FFTW planning arrays are allocated before setup timing. `FFTW_MEASURE` setup is cold: wisdom is forgotten before each measured forward+inverse plan pair.',
        '- `FFTW_ESTIMATE` and `FFTW_MEASURE` execution use persistent plans and preallocated buffers; planning is excluded from execution timing.',
        '- A MEASURE execution win is not automatically an end-to-end win. The break-even table quantifies how many repeated transforms are required to repay additional planning time.',
        '- Results describe this binary/library/virtualized machine. They do not imply FFTW will lead by the same factor on another architecture, library build, compiler, transform shape, thread count, or planning policy.',
        ''
    ]
    text = '\n'.join(lines) + '\n'
    if args.output:
        Path(args.output).write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()

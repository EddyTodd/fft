#!/usr/bin/env python3
"""Analyze raw fft timing samples without third-party dependencies."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


def percentile(xs: list[float], p: float) -> float:
    ys = sorted(xs)
    pos = p * (len(ys) - 1)
    lo, hi = math.floor(pos), math.ceil(pos)
    if lo == hi:
        return ys[lo]
    return ys[lo] * (hi - pos) + ys[hi] * (pos - lo)


def mad(xs: list[float]) -> float:
    m = statistics.median(xs)
    return statistics.median(abs(x - m) for x in xs)


def bootstrap_median(xs: list[float], rng: random.Random, reps: int) -> tuple[float, float]:
    n = len(xs)
    vals = [statistics.median(xs[rng.randrange(n)] for _ in range(n)) for _ in range(reps)]
    return percentile(vals, 0.025), percentile(vals, 0.975)


def bootstrap_speedup(
    winner: list[float], runner: list[float], rng: random.Random, reps: int
) -> tuple[float, float]:
    nw, nr = len(winner), len(runner)
    ratios = []
    for _ in range(reps):
        wm = statistics.median(winner[rng.randrange(nw)] for _ in range(nw))
        rm = statistics.median(runner[rng.randrange(nr)] for _ in range(nr))
        ratios.append(rm / wm)
    return percentile(ratios, 0.025), percentile(ratios, 0.975)


def common_language_fast(a: list[float], b: list[float]) -> float:
    wins = ties = 0
    for x in a:
        for y in b:
            if x < y:
                wins += 1
            elif x == y:
                ties += 1
    return (wins + 0.5 * ties) / (len(a) * len(b))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", type=Path, help="raw timing CSV/CSV.GZ files or directories containing timings-*.csv")
    ap.add_argument("--bootstrap", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=20260812)
    args = ap.parse_args()

    groups: dict[tuple[int, str], list[float]] = defaultdict(list)
    files: list[Path] = []
    for item in args.inputs:
        files.extend(sorted(item.rglob("timings-*.csv")) if item.is_dir() else [item])
    if not files:
        ap.error("no timing CSV files found")
    for path in files:
        opener = gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")
        with opener as fh:
            for row in csv.DictReader(fh):
                groups[(int(row["N"]), row["algorithm"])].append(float(row["ns_per_transform"]))

    rng = random.Random(args.seed)
    print("# Timing analysis\n")
    print("Medians are reported with a nonparametric 95% bootstrap CI. `P(faster)` is the common-language effect size: the probability that a randomly selected timing sample from the winner is faster than one from the runner-up.\n")
    print("| N | Winner | Median ns | 95% CI ns | MAD ns | Runner-up | Speedup | Speedup 95% CI | P(faster) |")
    print("|---:|---|---:|---:|---:|---|---:|---:|---:|")

    for n in sorted({n for n, _ in groups}):
        algos = [(a, groups[(n, a)]) for nn, a in groups if nn == n]
        algos.sort(key=lambda item: statistics.median(item[1]))
        if len(algos) < 2:
            continue
        (wa, wx), (ra, rx) = algos[:2]
        wm, rm = statistics.median(wx), statistics.median(rx)
        ci = bootstrap_median(wx, rng, args.bootstrap)
        sci = bootstrap_speedup(wx, rx, rng, args.bootstrap)
        prob = common_language_fast(wx, rx)
        print(
            f"| {n} | {wa} | {wm:.3f} | [{ci[0]:.3f}, {ci[1]:.3f}] | {mad(wx):.3f} | "
            f"{ra} | {rm / wm:.3f}x | [{sci[0]:.3f}x, {sci[1]:.3f}x] | {prob:.3f} |"
        )

    print("\n## Full ranking\n")
    print("| N | Rank | Algorithm | Median ns | MAD ns | p05 ns | p95 ns | Samples |")
    print("|---:|---:|---|---:|---:|---:|---:|---:|")
    for n in sorted({n for n, _ in groups}):
        algos = [(a, groups[(n, a)]) for nn, a in groups if nn == n]
        algos.sort(key=lambda item: statistics.median(item[1]))
        for rank, (a, xs) in enumerate(algos, 1):
            print(
                f"| {n} | {rank} | {a} | {statistics.median(xs):.3f} | {mad(xs):.3f} | "
                f"{percentile(xs, .05):.3f} | {percentile(xs, .95):.3f} | {len(xs)} |"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

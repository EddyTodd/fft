#!/usr/bin/env python3
"""Summarize FFT accuracy-suite CSV using robust across-signal statistics."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", type=Path, help="accuracy CSV files or directories containing accuracy-N*.csv")
    args = ap.parse_args()

    groups: dict[tuple[int, str], list[dict[str, str]]] = defaultdict(list)
    files: list[Path] = []
    for item in args.inputs:
        files.extend(sorted(item.rglob("accuracy-N*.csv")) if item.is_dir() else [item])
    if not files:
        ap.error("no accuracy CSV files found")
    for path in files:
        with path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                groups[(int(row["N"]), row["algorithm"])].append(row)

    print("# Accuracy analysis\n")
    print("Ranking uses the worst forward L2 error across the deterministic signal families. Lower is better. This is descriptive of the checked-in long-double reference experiment, not a proof of a global stability ordering.\n")
    print("| N | Rank | Algorithm | Median forward L2 | Worst forward L2 | Worst backward L2 | Worst forward Linf |")
    print("|---:|---:|---|---:|---:|---:|---:|")
    for n in sorted({n for n, _ in groups}):
        rows = []
        for (nn, algo), rs in groups.items():
            if nn != n:
                continue
            f2 = [float(r["forward_l2"]) for r in rs]
            b2 = [float(r["backward_l2"]) for r in rs]
            fi = [float(r["forward_linf"]) for r in rs]
            rows.append((max(f2), algo, statistics.median(f2), max(b2), max(fi)))
        rows.sort()
        for rank, (worst, algo, med, back, linf) in enumerate(rows, 1):
            print(f"| {n} | {rank} | {algo} | {med:.3e} | {worst:.3e} | {back:.3e} | {linf:.3e} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

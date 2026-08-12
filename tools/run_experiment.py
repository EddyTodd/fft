#!/usr/bin/env python3
"""Run a reproducible, randomized multi-session FFT experiment."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import platform
import random
import subprocess
from pathlib import Path


def prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def pow2(n: int) -> bool:
    return n > 0 and n & (n - 1) == 0


def algorithms(n: int) -> list[str]:
    out = ["auto", "mixed-radix", "bluestein"]
    if prime(n):
        out.append("rader")
    if pow2(n):
        out += ["radix2-iterative", "radix2-recursive", "stockham-radix2", "radix4", "split-radix"]
    if n <= 2048:
        out.append("dft")
    return out


def run(cmd: list[str]) -> str:
    return subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE).stdout


def cpu_description() -> str:
    if platform.system() == "Linux":
        try:
            for line in Path("/proc/cpuinfo").read_text().splitlines():
                if line.lower().startswith("model name"):
                    return line.split(":", 1)[1].strip()
        except OSError:
            pass
    return platform.processor() or "unknown"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, default=Path("build/fft"))
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--sizes", default="64,127,256,509,1009,1024,4096")
    ap.add_argument("--sessions", type=int, default=3)
    ap.add_argument("--samples", type=int, default=31)
    ap.add_argument("--target-ms", type=float, default=5.0)
    ap.add_argument("--warmups", type=int, default=5)
    ap.add_argument("--signal", default="tones")
    ap.add_argument("--accuracy-max", type=int, default=512, help="largest N included in the O(N^2) accuracy suite")
    ap.add_argument("--seed", type=int, default=20260812)
    args = ap.parse_args()

    binary = args.binary.resolve()
    sizes = [int(x) for x in args.sizes.split(",")]
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = args.out or Path("results") / f"run-{stamp}"
    out.mkdir(parents=True, exist_ok=False)

    git_commit = "unknown"
    try:
        git_commit = run(["git", "rev-parse", "HEAD"]).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    meta = {
        "schema": 1,
        "timestamp_utc": stamp,
        "git_commit": git_commit,
        "binary": str(binary),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu": cpu_description(),
        "python": platform.python_version(),
        "logical_cpus": os.cpu_count(),
        "sizes": sizes,
        "sessions": args.sessions,
        "samples_per_algorithm_per_session": args.samples,
        "target_ms": args.target_ms,
        "warmups": args.warmups,
        "signal": args.signal,
        "accuracy_max": args.accuracy_max,
        "randomization_seed": args.seed,
        "timing_semantics": "execution includes per-call allocations; setup/planning is not separated in this implementation",
    }
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")

    timing_path = out / "timings.csv"
    with timing_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["session", "order", "algorithm", "N", "sample", "iterations_per_sample", "ns_per_transform"])
        rng = random.Random(args.seed)
        for session in range(args.sessions):
            jobs = [(n, a) for n in sizes for a in algorithms(n)]
            rng.shuffle(jobs)
            for order, (n, algo) in enumerate(jobs):
                text = run([
                    str(binary), "--benchmark", "--algorithm", algo, "--size", str(n),
                    "--samples", str(args.samples), "--target-ms", str(args.target_ms),
                    "--warmups", str(args.warmups), "--signal", args.signal, "--raw-csv",
                ])
                rows = list(csv.DictReader(text.splitlines()))
                for row in rows:
                    writer.writerow([session, order, row["algorithm"], row["N"], row["sample"], row["iterations_per_sample"], row["ns_per_transform"]])
                fh.flush()

    accuracy_sizes = [n for n in sizes if n <= min(2048, args.accuracy_max)]
    if accuracy_sizes:
        (out / "accuracy.csv").write_text(run([
            str(binary), "--accuracy-suite", "--sizes", ",".join(map(str, accuracy_sizes)), "--csv"
        ]))

    with (out / "complexity.csv").open("w", newline="") as fh:
        writer = None
        for n in sizes:
            for algo in algorithms(n):
                text = run([str(binary), "--complexity", "--algorithm", algo, "--size", str(n), "--csv"])
                rows = list(csv.DictReader(text.splitlines()))
                if not rows:
                    continue
                if writer is None:
                    writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
                    writer.writeheader()
                writer.writerows(rows)

    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Benchmark a MembIT command and write machine-readable timing results.

The goal is to measure the practical cost of reading trajectories through each
input path, for example:

  python tools/benchmark_membit.py --label pdb --repeat 3 --out bench.tsv -- \
      python membit.py -f traj.pdb -n index.ndx -thickness 6 1 0 6 25 -deformation -o /tmp/pdb

  python tools/benchmark_membit.py --label xtc --repeat 3 --append --out bench.tsv -- \
      python membit.py -f traj.xtc -s structure.tpr -n index_full.ndx -thickness 6 1 0 6 25 -deformation -o /tmp/xtc

The script records wall time and CPU time.  On Linux it also records the child
process peak resident set size reported by resource.getrusage().
"""
from __future__ import annotations

import argparse
import csv
import json
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path


def run_once(command: list[str]) -> dict[str, float | int]:
    """Run one command and return timing/resource metrics.

    ``resource.getrusage(RUSAGE_CHILDREN)`` is cumulative for all completed
    child processes in the current Python process.  We take before/after
    snapshots and subtract CPU counters.  ``ru_maxrss`` is a maximum rather than
    an additive counter, so the reported value is the current cumulative maximum
    after the command.  That still works well for comparing runs launched from a
    fresh benchmark process.
    """
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()

    completed = subprocess.run(command)

    end = time.perf_counter()
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)

    return {
        "returncode": completed.returncode,
        "wall_seconds": end - start,
        "user_cpu_seconds": usage_after.ru_utime - usage_before.ru_utime,
        "sys_cpu_seconds": usage_after.ru_stime - usage_before.ru_stime,
        "max_rss_kb": usage_after.ru_maxrss,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark a MembIT command. Put the command after -- ."
    )
    parser.add_argument("--label", required=True, help="Short label, e.g. pdb, reduced_xtc, full_xtc")
    parser.add_argument("--repeat", type=int, default=3, help="Number of repeated runs")
    parser.add_argument("--out", type=Path, default=Path("membit_benchmark.tsv"), help="TSV output file")
    parser.add_argument("--json-out", type=Path, default=None, help="Optional JSON output file")
    parser.add_argument("--append", action="store_true", help="Append to an existing TSV instead of replacing it")
    parser.add_argument("command", nargs=argparse.REMAINDER, help="Command to benchmark, preceded by --")
    args = parser.parse_args()

    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise SystemExit("Missing command. Example: benchmark_membit.py --label pdb -- python membit.py ...")

    rows = []
    for i in range(1, args.repeat + 1):
        print(f"[{args.label}] run {i}/{args.repeat}: {' '.join(command)}", file=sys.stderr)
        metrics = run_once(command)
        row = {
            "label": args.label,
            "run": i,
            "command": " ".join(command),
            "python": sys.version.replace("\n", " "),
            "platform": platform.platform(),
            **metrics,
        }
        rows.append(row)
        if metrics["returncode"] != 0:
            # Stop immediately on failure.  A failed benchmark is usually a bad
            # input path or an analysis exception, not a performance data point.
            break

    fieldnames = [
        "label",
        "run",
        "returncode",
        "wall_seconds",
        "user_cpu_seconds",
        "sys_cpu_seconds",
        "max_rss_kb",
        "command",
        "python",
        "platform",
    ]

    write_header = not args.append or not args.out.exists()
    mode = "a" if args.append else "w"
    with args.out.open(mode, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if write_header:
            writer.writeheader()
        writer.writerows(rows)

    if args.json_out:
        args.json_out.write_text(json.dumps(rows, indent=2))

    # Return the last command status so CI/HPC scripts can fail naturally.
    return int(rows[-1]["returncode"])


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Numerically compare two MembIT XVG output files.

This helper is intentionally small and dependency-free so it can be used on
HPC/login nodes where only the MembIT Python environment is available.  It
ignores XVG comments/header lines, parses all numeric columns, treats matching
NaN values as equal, and exits with a non-zero status when a difference exceeds
the requested tolerance.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path


def parse_xvg(path: Path):
    """Return numeric rows from an XVG-like text file.

    MembIT output is mostly whitespace-separated numeric data, but may also
    contain comments, Grace/XVG metadata, or the string "No occurrences".  Those
    non-data lines are skipped so regression tests focus on the calculated
    values.
    """
    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            if "No occurrences" in line:
                continue

            values = []
            for token in line.split():
                if token.lower() == "nan":
                    values.append(math.nan)
                    continue
                try:
                    values.append(float(token))
                except ValueError:
                    # Header tokens such as "time" are ignored.
                    pass

            if values:
                rows.append(values)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare two MembIT XVG output files numerically.")
    parser.add_argument("reference", type=Path, help="Known-good XVG file")
    parser.add_argument("candidate", type=Path, help="Newly generated XVG file")
    parser.add_argument("--atol", type=float, default=1e-3, help="Absolute tolerance per numeric value")
    args = parser.parse_args()

    ref = parse_xvg(args.reference)
    cand = parse_xvg(args.candidate)

    if len(ref) != len(cand):
        raise SystemExit(f"row count differs: reference={len(ref)} candidate={len(cand)}")

    worst = 0.0
    for i, (rrow, crow) in enumerate(zip(ref, cand), start=1):
        if len(rrow) != len(crow):
            raise SystemExit(f"column count differs on row {i}: reference={len(rrow)} candidate={len(crow)}")

        for j, (r, c) in enumerate(zip(rrow, crow), start=1):
            if math.isnan(r) and math.isnan(c):
                continue
            delta = abs(r - c)
            worst = max(worst, delta)
            if delta > args.atol:
                raise SystemExit(
                    f"values differ on row {i}, column {j}: "
                    f"reference={r} candidate={c} delta={delta} > atol={args.atol}"
                )

    print(f"OK: {args.reference} and {args.candidate} match within atol={args.atol}; worst_delta={worst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

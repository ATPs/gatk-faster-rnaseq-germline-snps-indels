#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import DEFAULT_REF, RUST_HC_PREFILTER, existing, resolve_rust_binary, run_command


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build HaplotypeCaller candidate intervals from BAM pileup evidence."
    )
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-interval-list", type=existing, required=True)
    parser.add_argument("--output-interval-list", type=Path, required=True)
    parser.add_argument("--output-summary", type=Path, required=True)
    parser.add_argument("--output-bed", type=Path)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=10)
    parser.add_argument("--min-alt-count", type=int, default=1)
    parser.add_argument("--min-indel-count", type=int, default=1)
    parser.add_argument("--min-alt-fraction", type=float, default=0.0)
    parser.add_argument("--padding", type=int, default=150)
    parser.add_argument("--max-depth", type=int, default=100000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--exclude-supplementary", action="store_true")
    parser.add_argument("--empty-behavior", choices=("error", "input"), default="input")
    parser.add_argument("--rust-bin", type=Path)
    args = parser.parse_args()

    if args.min_mapq < 0 or args.min_baseq < 0:
        raise SystemExit("--min-mapq and --min-baseq must be non-negative.")
    if args.min_alt_count < 1 or args.min_indel_count < 1:
        raise SystemExit("--min-alt-count and --min-indel-count must be at least 1.")
    if not 0.0 <= args.min_alt_fraction <= 1.0:
        raise SystemExit("--min-alt-fraction must be between 0 and 1.")
    if args.padding < 0:
        raise SystemExit("--padding must be non-negative.")
    if args.max_depth < 1 or args.threads < 1:
        raise SystemExit("--max-depth and --threads must be at least 1.")
    rust_bin = resolve_rust_binary(args.rust_bin, RUST_HC_PREFILTER)
    if not rust_bin.exists():
        raise SystemExit(
            f"Rust binary not found: {rust_bin}. "
            "Build it with: cd src/rust && /data/p/sys/rust/1.96.0/bin/cargo build --release --bin rust_hc_prefilter"
        )

    args.output_interval_list.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(rust_bin),
        "candidate-intervals",
        "--input-bam",
        str(args.input_bam),
        "--ref",
        str(args.ref),
        "--input-interval-list",
        str(args.input_interval_list),
        "--output-interval-list",
        str(args.output_interval_list),
        "--output-summary",
        str(args.output_summary),
        "--min-mapq",
        str(args.min_mapq),
        "--min-baseq",
        str(args.min_baseq),
        "--min-alt-count",
        str(args.min_alt_count),
        "--min-indel-count",
        str(args.min_indel_count),
        "--min-alt-fraction",
        str(args.min_alt_fraction),
        "--padding",
        str(args.padding),
        "--max-depth",
        str(args.max_depth),
        "--threads",
        str(args.threads),
        "--empty-behavior",
        args.empty_behavior,
    ]
    if args.exclude_supplementary:
        command.append("--exclude-supplementary")
    if args.output_bed is not None:
        args.output_bed.parent.mkdir(parents=True, exist_ok=True)
        command.extend(["--output-bed", str(args.output_bed)])

    run_command(command)
    print(args.output_interval_list)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from .common import (
    DEFAULT_KNOWN_SITES,
    DEFAULT_REF,
    RUST_BASE_RECALIBRATOR,
    existing,
    gatk_cmd,
    resolve_rust_binary,
    run_command,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run BaseRecalibrator.")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--known-sites", type=existing, action="append")
    parser.add_argument("--output-table", type=Path, required=True)
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--backend", choices=("auto", "gatk", "rust"), default="auto")
    parser.add_argument("--rust-bin", type=Path)
    parser.add_argument("--mismatches-context-size", type=int, default=2)
    parser.add_argument("--low-quality-tail", type=int, default=2)
    parser.add_argument("--maximum-cycle-value", type=int, default=500)
    parser.add_argument("--quantizing-levels", type=int, default=16)
    parser.add_argument("--known-sites-chunk-size", type=int, default=1_000_000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--region-bases", type=int, default=25_000_000)
    args = parser.parse_args()

    known_sites = args.known_sites or [path for path in DEFAULT_KNOWN_SITES if path.exists()]
    if not known_sites:
        raise SystemExit("BaseRecalibrator needs at least one --known-sites VCF.")
    if args.threads < 1:
        raise SystemExit("--threads must be at least 1.")
    if args.region_bases < 1:
        raise SystemExit("--region-bases must be at least 1.")

    args.output_table.parent.mkdir(parents=True, exist_ok=True)
    rust_bin = resolve_rust_binary(args.rust_bin, RUST_BASE_RECALIBRATOR)
    backend = args.backend
    if backend == "auto":
        backend = "rust" if rust_bin.exists() else "gatk"
    if backend == "gatk":
        command = gatk_cmd(
            args.java_mem,
            "BaseRecalibrator",
            "-R",
            args.ref,
            "-I",
            args.input_bam,
            "--use-original-qualities",
            "-O",
            args.output_table,
        )
        command.extend(["--mismatches-context-size", str(args.mismatches_context_size)])
        command.extend(["--low-quality-tail", str(args.low_quality_tail)])
        command.extend(["--maximum-cycle-value", str(args.maximum_cycle_value)])
        command.extend(["--quantizing-levels", str(args.quantizing_levels)])
    else:
        if not rust_bin.exists():
            raise SystemExit(f"rust backend binary not found: {rust_bin}")
        command = [
            str(rust_bin),
            "--ref",
            str(args.ref),
            "--input-bam",
            str(args.input_bam),
            "--use-original-qualities",
            "--output-table",
            str(args.output_table),
            "--mismatches-context-size",
            str(args.mismatches_context_size),
            "--low-quality-tail",
            str(args.low_quality_tail),
            "--maximum-cycle-value",
            str(args.maximum_cycle_value),
            "--quantizing-levels",
            str(args.quantizing_levels),
            "--known-sites-chunk-size",
            str(args.known_sites_chunk_size),
            "--threads",
            str(args.threads),
            "--region-bases",
            str(args.region_bases),
        ]
    for known in known_sites:
        command.extend(["--known-sites", str(known)])
    run_command(command)
    print(args.output_table)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

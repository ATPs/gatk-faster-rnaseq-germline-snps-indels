#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from .common import (
    DEFAULT_REF,
    RUST_SPLIT_N_CIGAR_READS,
    existing,
    gatk_cmd,
    resolve_rust_binary,
    run_command,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run GATK or Rust SplitNCigarReads.")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--output-bam", type=Path, required=True)
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--backend", choices=("auto", "gatk", "rust"), default="auto")
    parser.add_argument("--rust-bin", type=Path)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--rust-mode", choices=("fast", "compatibility"), default="fast")
    parser.add_argument("--skip-mapping-quality-transform", action="store_true")
    parser.add_argument("--process-secondary-alignments", action="store_true")
    args = parser.parse_args()

    args.output_bam.parent.mkdir(parents=True, exist_ok=True)
    rust_bin = resolve_rust_binary(args.rust_bin, RUST_SPLIT_N_CIGAR_READS)
    backend = args.backend
    if backend == "auto":
        backend = "rust" if rust_bin.exists() else "gatk"
    if backend == "gatk":
        command = gatk_cmd(args.java_mem, "SplitNCigarReads", "-R", args.ref, "-I", args.input_bam, "-O", args.output_bam)
    else:
        if args.threads < 1:
            raise SystemExit("--threads must be at least 1.")
        if not rust_bin.exists():
            raise SystemExit(f"rust backend binary not found: {rust_bin}")
        command = [
            str(rust_bin),
            "--input-bam",
            str(args.input_bam),
            "--output-bam",
            str(args.output_bam),
            "--threads",
            str(args.threads),
            "--mode",
            args.rust_mode,
        ]
        if args.rust_mode == "compatibility":
            command.extend(["--reference", str(args.ref)])
        if args.skip_mapping_quality_transform:
            command.append("--skip-mapping-quality-transform")
        if args.process_secondary_alignments:
            command.append("--process-secondary-alignments")
    run_command(command)
    print(args.output_bam)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

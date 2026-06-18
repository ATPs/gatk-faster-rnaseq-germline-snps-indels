#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from .common import (
    RUST_MARK_DUPLICATES,
    SAMBAMBA,
    SAMTOOLS,
    bam_index_path,
    existing,
    gatk_cmd,
    resolve_rust_binary,
    run_command,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run duplicate marking with Picard/GATK, sambamba, or Rust.")
    parser.add_argument("--mode", choices=("auto", "baseline", "sambamba", "rust"), default="auto")
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--output-bam", type=Path, required=True)
    parser.add_argument("--output-metrics", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=40)
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--tmpdir", type=Path, required=True)
    parser.add_argument("--rust-bin", type=Path)
    args = parser.parse_args()

    args.output_bam.parent.mkdir(parents=True, exist_ok=True)
    args.output_metrics.parent.mkdir(parents=True, exist_ok=True)
    args.tmpdir.mkdir(parents=True, exist_ok=True)

    rust_bin = resolve_rust_binary(args.rust_bin, RUST_MARK_DUPLICATES)
    mode = args.mode
    if mode == "auto":
        mode = "rust" if rust_bin.exists() else "baseline"

    if mode == "baseline":
        run_command(
            gatk_cmd(
                args.java_mem,
                "MarkDuplicates",
                "-I",
                args.input_bam,
                "-O",
                args.output_bam,
                "-M",
                args.output_metrics,
                "--CREATE_INDEX",
                "true",
                "--VALIDATION_STRINGENCY",
                "SILENT",
            )
        )
    elif mode == "sambamba":
        run_command(
            [
                str(SAMBAMBA),
                "markdup",
                "-t",
                str(args.threads),
                "--tmpdir",
                str(args.tmpdir),
                str(args.input_bam),
                str(args.output_bam),
            ]
        )
        run_command([str(SAMBAMBA), "index", "-t", str(args.threads), str(args.output_bam)])
    else:
        if not rust_bin.exists():
            raise SystemExit(
                f"Rust binary not found: {rust_bin}. "
                "Build it with: cd src/rust && /data/p/sys/rust/1.96.0/bin/cargo build --release --bin rust_mark_duplicates"
            )
        run_command(
            [
                str(rust_bin),
                "--input-bam",
                args.input_bam,
                "--output-bam",
                args.output_bam,
                "--output-metrics",
                args.output_metrics,
                "--samtools",
                SAMTOOLS,
                "--threads",
                str(args.threads),
            ]
        )

    print(args.output_bam)
    print(bam_index_path(args.output_bam))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

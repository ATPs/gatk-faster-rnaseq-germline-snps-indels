#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import (
    DEFAULT_REF,
    RUST_APPLY_BQSR,
    existing,
    gatk_cmd,
    resolve_rust_binary,
    run_command,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run ApplyBQSR.")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--input-table", type=existing, required=True)
    parser.add_argument("--output-bam", type=Path, required=True)
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--backend", choices=("auto", "gatk", "rust"), default="auto")
    parser.add_argument("--rust-bin", type=Path)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    args.output_bam.parent.mkdir(parents=True, exist_ok=True)
    rust_bin = resolve_rust_binary(args.rust_bin, RUST_APPLY_BQSR)
    backend = args.backend
    if backend == "auto":
        backend = "rust" if rust_bin.exists() else "gatk"
    if backend == "gatk":
        command = gatk_cmd(
            args.java_mem,
            "ApplyBQSR",
            "-R",
            args.ref,
            "-I",
            args.input_bam,
            "--use-original-qualities",
            "--bqsr-recal-file",
            args.input_table,
            "-O",
            args.output_bam,
        )
    else:
        if not rust_bin.exists():
            raise SystemExit(
                f"Rust binary not found: {rust_bin}. "
                "Build it with: cd src/rust && /data/p/sys/rust/1.96.0/bin/cargo build --release --bin rust_apply_bqsr"
            )
        command = [
            str(rust_bin),
            "--ref",
            str(args.ref),
            "--input-bam",
            str(args.input_bam),
            "--input-table",
            str(args.input_table),
            "--output-bam",
            str(args.output_bam),
            "--threads",
            str(args.threads),
            "--use-original-qualities",
        ]
    run_command(command)
    print(args.output_bam)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

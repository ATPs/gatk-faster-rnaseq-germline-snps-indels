#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from .common import (
    DEFAULT_REF,
    RUST_INTERVAL_TOOLS,
    existing,
    gatk_cmd,
    resolve_rust_binary,
    run_command,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Split HaplotypeCaller interval_list into scattered shards.")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-interval-list", type=existing, required=True)
    parser.add_argument("--scatter-count", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--backend", choices=("auto", "gatk", "rust"), default="auto")
    parser.add_argument("--rust-bin", type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rust_bin = resolve_rust_binary(args.rust_bin, RUST_INTERVAL_TOOLS)
    backend = args.backend
    if backend == "auto":
        backend = "rust" if rust_bin.exists() else "gatk"
    if backend == "gatk":
        command = gatk_cmd(
            args.java_mem,
            "SplitIntervals",
            "-R",
            args.ref,
            "-L",
            args.input_interval_list,
            "--scatter-count",
            str(args.scatter_count),
            "--subdivision-mode",
            "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW",
            "-O",
            args.output_dir,
        )
    else:
        if not rust_bin.exists():
            raise SystemExit(f"rust backend binary not found: {rust_bin}")
        command = [
            str(rust_bin),
            "split-intervals",
            "--input-interval-list",
            str(args.input_interval_list),
            "--scatter-count",
            str(args.scatter_count),
            "--output-dir",
            str(args.output_dir),
        ]
    run_command(command)
    print(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

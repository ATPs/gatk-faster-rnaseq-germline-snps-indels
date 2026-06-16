#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import (
    DEFAULT_REF_DICT,
    RUST_INTERVAL_TOOLS,
    existing,
    gatk_cmd,
    resolve_rust_binary,
    run_command,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert BED to Picard/GATK interval_list.")
    parser.add_argument("--input-bed", type=existing, required=True)
    parser.add_argument("--ref-dict", type=existing, default=DEFAULT_REF_DICT)
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--output-interval-list", type=Path, required=True)
    parser.add_argument("--backend", choices=("auto", "gatk", "rust"), default="auto")
    parser.add_argument("--rust-bin", type=Path)
    args = parser.parse_args()

    args.output_interval_list.parent.mkdir(parents=True, exist_ok=True)
    rust_bin = resolve_rust_binary(args.rust_bin, RUST_INTERVAL_TOOLS)
    backend = args.backend
    if backend == "auto":
        backend = "rust" if rust_bin.exists() else "gatk"
    if backend == "gatk":
        command = gatk_cmd(
            args.java_mem,
            "BedToIntervalList",
            "-I",
            args.input_bed,
            "-O",
            args.output_interval_list,
            "-SD",
            args.ref_dict,
        )
    else:
        if not rust_bin.exists():
            raise SystemExit(f"rust backend binary not found: {rust_bin}")
        command = [
            str(rust_bin),
            "bed-to-interval-list",
            "--input-bed",
            str(args.input_bed),
            "--ref-dict",
            str(args.ref_dict),
            "--output-interval-list",
            str(args.output_interval_list),
        ]
    run_command(command)
    print(args.output_interval_list)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

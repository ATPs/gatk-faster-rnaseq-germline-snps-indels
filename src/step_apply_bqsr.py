#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import DEFAULT_REF, existing, gatk_cmd, run_command


def main() -> int:
    parser = argparse.ArgumentParser(description="Run GATK ApplyBQSR.")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--input-table", type=existing, required=True)
    parser.add_argument("--output-bam", type=Path, required=True)
    parser.add_argument("--java-mem", default="24g")
    args = parser.parse_args()

    args.output_bam.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        gatk_cmd(
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
    )
    print(args.output_bam)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

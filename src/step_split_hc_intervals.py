#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import DEFAULT_REF, existing, gatk_cmd, run_command


def main() -> int:
    parser = argparse.ArgumentParser(description="Split HaplotypeCaller interval_list into scattered shards.")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-interval-list", type=existing, required=True)
    parser.add_argument("--scatter-count", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--java-mem", default="24g")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    run_command(
        gatk_cmd(
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
    )
    print(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

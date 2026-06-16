#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import SAMBAMBA, bam_index_path, existing, gatk_cmd, run_command


def main() -> int:
    parser = argparse.ArgumentParser(description="Run duplicate marking with Picard/GATK or sambamba.")
    parser.add_argument("--mode", choices=("baseline", "sambamba"), default="baseline")
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--output-bam", type=Path, required=True)
    parser.add_argument("--output-metrics", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=40)
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--tmpdir", type=Path, required=True)
    args = parser.parse_args()

    args.output_bam.parent.mkdir(parents=True, exist_ok=True)
    args.output_metrics.parent.mkdir(parents=True, exist_ok=True)
    args.tmpdir.mkdir(parents=True, exist_ok=True)

    if args.mode == "baseline":
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
    else:
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

    print(args.output_bam)
    print(bam_index_path(args.output_bam))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

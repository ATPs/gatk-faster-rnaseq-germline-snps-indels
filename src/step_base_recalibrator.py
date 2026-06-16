#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import DEFAULT_KNOWN_SITES, DEFAULT_REF, existing, gatk_cmd, run_command


def main() -> int:
    parser = argparse.ArgumentParser(description="Run GATK BaseRecalibrator.")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--known-sites", type=existing, action="append")
    parser.add_argument("--output-table", type=Path, required=True)
    parser.add_argument("--java-mem", default="24g")
    args = parser.parse_args()

    known_sites = args.known_sites or [path for path in DEFAULT_KNOWN_SITES if path.exists()]
    if not known_sites:
        raise SystemExit("BaseRecalibrator needs at least one --known-sites VCF.")

    args.output_table.parent.mkdir(parents=True, exist_ok=True)
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
    for known in known_sites:
        command.extend(["--known-sites", str(known)])
    run_command(command)
    print(args.output_table)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import existing, gatk_cmd, run_command


def main() -> int:
    parser = argparse.ArgumentParser(description="Merge scattered VCF shards.")
    parser.add_argument("--input-vcf", type=existing, action="append", required=True)
    parser.add_argument("--output-vcf", type=Path, required=True)
    parser.add_argument("--java-mem", default="24g")
    args = parser.parse_args()

    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    command: list[str | Path] = ["MergeVcfs"]
    for vcf in args.input_vcf:
        command.extend(["-I", vcf])
    command.extend(["-O", args.output_vcf])
    run_command(gatk_cmd(args.java_mem, *command))
    print(args.output_vcf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

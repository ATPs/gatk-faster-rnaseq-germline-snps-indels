#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import DEFAULT_REF, existing, gatk_cmd, run_command


def main() -> int:
    parser = argparse.ArgumentParser(description="Run GATK VariantFiltration for the RNA hard filters.")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-vcf", type=existing, required=True)
    parser.add_argument("--output-vcf", type=Path, required=True)
    parser.add_argument("--java-mem", default="24g")
    args = parser.parse_args()

    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        gatk_cmd(
            args.java_mem,
            "VariantFiltration",
            "-R",
            args.ref,
            "-V",
            args.input_vcf,
            "--window",
            "35",
            "--cluster",
            "3",
            "--filter-name",
            "FS",
            "--filter",
            "FS > 30.0",
            "--filter-name",
            "QD",
            "--filter",
            "QD < 2.0",
            "-O",
            args.output_vcf,
        )
    )
    print(args.output_vcf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

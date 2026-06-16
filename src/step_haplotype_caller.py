#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import DEFAULT_DBSNP, DEFAULT_REF, build_haplotype_caller_command, existing, run_command


def main() -> int:
    parser = argparse.ArgumentParser(description="Run GATK HaplotypeCaller on one interval set.")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--input-interval-list", type=existing, required=True)
    parser.add_argument("--output-vcf", type=Path, required=True)
    parser.add_argument("--pair-hmm-threads", type=int, default=8)
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--dbsnp", type=existing, default=DEFAULT_DBSNP)
    args = parser.parse_args()

    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        build_haplotype_caller_command(
            args.ref,
            args.input_bam,
            args.input_interval_list,
            args.output_vcf,
            args.java_mem,
            args.pair_hmm_threads,
            args.dbsnp,
        )
    )
    print(args.output_vcf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

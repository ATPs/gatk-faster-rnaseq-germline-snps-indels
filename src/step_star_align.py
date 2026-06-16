#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import DEFAULT_FASTQ_DIR, DEFAULT_STAR_INDEX, STAR, existing, run_command


def main() -> int:
    parser = argparse.ArgumentParser(description="Run STAR two-pass RNA-seq alignment.")
    parser.add_argument("--sample", default="SRR949115")
    parser.add_argument("--fastq1", type=existing, default=DEFAULT_FASTQ_DIR / "SRR949115_1.fastq.gz")
    parser.add_argument("--fastq2", type=existing, default=DEFAULT_FASTQ_DIR / "SRR949115_2.fastq.gz")
    parser.add_argument("--star-index", type=existing, default=DEFAULT_STAR_INDEX)
    parser.add_argument("--threads", type=int, default=40)
    parser.add_argument("--output-prefix", type=Path, required=True)
    args = parser.parse_args()

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        [
            str(STAR),
            "--runThreadN",
            str(args.threads),
            "--genomeDir",
            str(args.star_index),
            "--readFilesIn",
            str(args.fastq1),
            str(args.fastq2),
            "--readFilesCommand",
            "zcat",
            "--twopassMode",
            "Basic",
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--limitBAMsortRAM",
            "45000000000",
            "--outSAMattrRGline",
            f"ID:{args.sample}.rg1",
            f"SM:{args.sample}",
            f"LB:{args.sample}.lib1",
            "PL:ILLUMINA",
            f"PU:{args.sample}.unit1",
            "--outFileNamePrefix",
            str(args.output_prefix),
        ]
    )
    print(f"{args.output_prefix}Aligned.sortedByCoord.out.bam")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

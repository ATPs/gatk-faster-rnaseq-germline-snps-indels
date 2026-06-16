#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import DEFAULT_GTF, DEFAULT_REF, SAMTOOLS, build_exon_bed, existing, run_command


def main() -> int:
    parser = argparse.ArgumentParser(description="Build exon BED from GTF and reference .fai.")
    parser.add_argument("--gtf", type=existing, default=DEFAULT_GTF)
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--output-bed", type=Path, required=True)
    args = parser.parse_args()

    fai_path = Path(f"{args.ref}.fai")
    if not fai_path.exists():
        run_command([str(SAMTOOLS), "faidx", str(args.ref)])
    build_exon_bed(args.gtf, fai_path, args.output_bed)
    print(args.output_bed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

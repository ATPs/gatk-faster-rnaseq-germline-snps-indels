#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from .common import build_merged_exon_bed, existing


def main() -> int:
    parser = argparse.ArgumentParser(description="Merge a sorted exon BED into non-overlapping intervals.")
    parser.add_argument("--input-bed", type=existing, required=True)
    parser.add_argument("--output-bed", type=Path, required=True)
    args = parser.parse_args()

    build_merged_exon_bed(args.input_bed, args.output_bed)
    print(args.output_bed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

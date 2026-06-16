#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import DEFAULT_REF_DICT, existing, gatk_cmd, run_command


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert BED to Picard/GATK interval_list.")
    parser.add_argument("--input-bed", type=existing, required=True)
    parser.add_argument("--ref-dict", type=existing, default=DEFAULT_REF_DICT)
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--output-interval-list", type=Path, required=True)
    args = parser.parse_args()

    args.output_interval_list.parent.mkdir(parents=True, exist_ok=True)
    run_command(gatk_cmd(args.java_mem, "BedToIntervalList", "-I", args.input_bed, "-O", args.output_interval_list, "-SD", args.ref_dict))
    print(args.output_interval_list)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

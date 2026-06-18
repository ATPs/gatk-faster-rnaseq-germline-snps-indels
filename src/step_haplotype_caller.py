#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

from step_common import (
    DEFAULT_DBSNP,
    DEFAULT_REF,
    RUST_HAPLOTYPE_CALLER,
    build_haplotype_caller_command,
    existing,
    resolve_rust_binary,
    run_command,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run GATK or Rust HaplotypeCaller on one interval set.")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--input-interval-list", type=existing, required=True)
    parser.add_argument("--output-vcf", type=Path, required=True)
    parser.add_argument("--backend", choices=("gatk", "rust"), default="gatk")
    parser.add_argument("--rust-bin", type=Path)
    parser.add_argument("--threads", type=int, default=40)
    parser.add_argument("--memory-gb", type=int, default=128)
    parser.add_argument("--pair-hmm-threads", type=int, default=8)
    parser.add_argument("--pair-hmm-implementation", choices=("rust", "native"), default="native")
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--dbsnp", type=existing, default=DEFAULT_DBSNP)
    parser.add_argument("--exclude-supplementary", action="store_true")
    args = parser.parse_args()

    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    if args.backend == "gatk":
        command = build_haplotype_caller_command(
            args.ref,
            args.input_bam,
            args.input_interval_list,
            args.output_vcf,
            args.java_mem,
            args.pair_hmm_threads,
            args.dbsnp,
        )
    else:
        if args.threads < 1:
            raise SystemExit("--threads must be at least 1.")
        if args.memory_gb < 1:
            raise SystemExit("--memory-gb must be at least 1.")
        if args.pair_hmm_threads < 1:
            raise SystemExit("--pair-hmm-threads must be at least 1.")
        rust_bin = resolve_rust_binary(args.rust_bin, RUST_HAPLOTYPE_CALLER)
        if not rust_bin.exists():
            raise SystemExit(
                f"rust backend binary not found: {rust_bin}. "
                "Build it with: cd src/rust && /data/p/sys/rust/1.96.0/bin/cargo build --release --bin rust_haplotype_caller"
            )
        command = [
            str(rust_bin),
            "call",
            "--input-bam",
            str(args.input_bam),
            "--ref",
            str(args.ref),
            "--input-interval-list",
            str(args.input_interval_list),
            "--output-vcf",
            str(args.output_vcf),
            "--dont-use-soft-clipped-bases",
            "--standard-min-confidence-threshold-for-calling",
            "20",
            "--threads",
            str(args.threads),
            "--memory-gb",
            str(args.memory_gb),
            "--native-pair-hmm-threads",
            str(args.pair_hmm_threads),
            "--pair-hmm-implementation",
            args.pair_hmm_implementation,
        ]
        if args.exclude_supplementary:
            command.append("--exclude-supplementary")
        if args.dbsnp is not None:
            command.extend(["--dbsnp", str(args.dbsnp)])
    run_command(command)
    print(args.output_vcf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

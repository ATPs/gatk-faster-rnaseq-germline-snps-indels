#!/usr/bin/env python3
"""Run and time the RNA-seq germline SNP/indel pipeline."""

from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Sequence
from pathlib import Path

from ..runtime import Step, run_parallel_steps, run_step
from ..steps.common import (
    DEFAULT_DBSNP,
    DEFAULT_FASTQ_DIR,
    DEFAULT_GTF,
    DEFAULT_KNOWN_SITES,
    DEFAULT_REF,
    DEFAULT_REF_DICT,
    DEFAULT_STAR_INDEX,
    SAMTOOLS,
    bam_index_path,
    existing,
    gatk_cmd,
)


REPO_SRC_DIR = Path(__file__).resolve().parents[2]
DEFAULT_RUST_BIN_DIR = REPO_SRC_DIR / "rust" / "bin"
PYTHON_STEP_MODULE_PREFIX = "gatk_faster_rnaseq.steps"


def build_python_step_env() -> dict[str, str]:
    pythonpath = [str(REPO_SRC_DIR)]
    current = os.environ.get("PYTHONPATH")
    if current:
        pythonpath.append(current)
    return {"PYTHONPATH": os.pathsep.join(pythonpath)}


PYTHON_STEP_ENV = build_python_step_env()


def default_rust_binary(bin_dir: Path, name: str) -> Path:
    return bin_dir / name


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a timed GATK RNA-seq germline SNP/indel pipeline."
    )
    parser.add_argument("--sample", default="SRR949115")
    parser.add_argument("--fastq1", type=existing, default=DEFAULT_FASTQ_DIR / "SRR949115_1.fastq.gz")
    parser.add_argument("--fastq2", type=existing, default=DEFAULT_FASTQ_DIR / "SRR949115_2.fastq.gz")
    parser.add_argument("--aligned-bam", type=existing, help="Use an existing STAR coordinate-sorted BAM.")
    parser.add_argument("--dedup-bam", type=existing, help="Use an existing duplicate-marked BAM when skipping mark duplicates.")
    parser.add_argument("--split-bam", type=existing, help="Use an existing SplitNCigarReads BAM when skipping that step.")
    parser.add_argument("--recal-table", type=existing, help="Use an existing recalibration table when skipping BaseRecalibrator.")
    parser.add_argument("--recal-bam", type=existing, help="Use an existing post-BQSR BAM when skipping ApplyBQSR.")
    parser.add_argument("--raw-vcf", type=existing, help="Use an existing raw VCF when skipping HaplotypeCaller.")
    parser.add_argument("--outdir", type=Path, default=Path("/XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115"))
    parser.add_argument("--mode", choices=("auto", "baseline", "sambamba", "rust"), default="auto")
    parser.add_argument("--threads", type=int, default=40)
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--ref-dict", type=existing, default=DEFAULT_REF_DICT)
    parser.add_argument("--gtf", type=existing, default=DEFAULT_GTF)
    parser.add_argument("--star-index", type=existing, default=DEFAULT_STAR_INDEX)
    parser.add_argument("--dbsnp", type=existing, default=DEFAULT_DBSNP)
    parser.add_argument("--known-sites", type=existing, action="append")
    parser.add_argument(
        "--interval-policy",
        choices=("raw", "merged"),
        default="raw",
        help="Use raw exon intervals or merged non-overlapping exon intervals for HaplotypeCaller.",
    )
    parser.add_argument(
        "--rust-bin-dir",
        type=Path,
        default=DEFAULT_RUST_BIN_DIR,
        help="Directory containing release Rust pipeline binaries.",
    )
    parser.add_argument(
        "--no-rust",
        action="store_true",
        help="Disable all Rust backends and force the GATK/non-Rust path.",
    )
    parser.add_argument(
        "--interval-backend",
        choices=("auto", "gatk", "rust"),
        default="auto",
        help="Backend for BED to interval_list conversion and HaplotypeCaller interval splitting.",
    )
    parser.add_argument(
        "--rust-interval-tools",
        type=Path,
        help="Path to the rust_interval_tools binary when --interval-backend=rust.",
    )
    parser.add_argument(
        "--split-n-cigar-backend",
        choices=("auto", "gatk", "rust"),
        default="auto",
        help="Backend for SplitNCigarReads.",
    )
    parser.add_argument(
        "--rust-split-n-cigar",
        type=Path,
        help="Path to the rust_split_n_cigar_reads binary when --split-n-cigar-backend=rust.",
    )
    parser.add_argument(
        "--split-n-cigar-rust-mode",
        choices=("fast", "compatibility"),
        default="fast",
        help="Rust SplitNCigarReads mode. Fast mode is validated for speed-first use; compatibility mode repairs more tags.",
    )
    parser.add_argument(
        "--base-recalibrator-backend",
        choices=("auto", "gatk", "rust"),
        default="auto",
        help="Backend for BaseRecalibrator.",
    )
    parser.add_argument(
        "--rust-base-recalibrator",
        type=Path,
        help="Path to the rust_base_recalibrator binary when --base-recalibrator-backend=rust.",
    )
    parser.add_argument(
        "--base-recalibrator-region-bases",
        type=int,
        default=25_000_000,
        help="Genomic region size used by the Rust BaseRecalibrator work queue.",
    )
    parser.add_argument(
        "--apply-bqsr-backend",
        choices=("auto", "gatk", "rust"),
        default="auto",
        help="Backend for ApplyBQSR.",
    )
    parser.add_argument(
        "--rust-apply-bqsr",
        type=Path,
        help="Path to the rust_apply_bqsr binary when --apply-bqsr-backend=rust.",
    )
    parser.add_argument(
        "--hc-scatter-count",
        type=int,
        default=1,
        help="Split HaplotypeCaller intervals into this many shards and run them in parallel.",
    )
    parser.add_argument(
        "--hc-threads-per-shard",
        type=int,
        help="Native PairHMM threads to assign to each scattered HaplotypeCaller shard.",
    )
    parser.add_argument(
        "--hc-parallel-jobs",
        type=int,
        default=4,
        help="Maximum number of scattered HaplotypeCaller Java processes to run at once.",
    )
    parser.add_argument(
        "--hc-backend",
        choices=("gatk", "rust"),
        default="gatk",
        help="Backend for HaplotypeCaller. The Rust full caller is opt-in while it is under development.",
    )
    parser.add_argument(
        "--rust-haplotype-caller",
        type=Path,
        help="Path to the rust_haplotype_caller binary when --hc-backend=rust.",
    )
    parser.add_argument(
        "--hc-memory-gb",
        type=int,
        default=128,
        help="Memory budget advertised to rust_haplotype_caller.",
    )
    parser.add_argument(
        "--hc-pair-hmm-implementation",
        choices=("rust", "native"),
        default="native",
        help="PairHMM implementation requested from rust_haplotype_caller.",
    )
    parser.add_argument(
        "--hc-prefilter-backend",
        choices=("auto", "none", "rust"),
        default="auto",
        help="Optionally prefilter HaplotypeCaller intervals before calling.",
    )
    parser.add_argument(
        "--rust-hc-prefilter",
        type=Path,
        help="Path to the rust_hc_prefilter binary when --hc-prefilter-backend=rust.",
    )
    parser.add_argument("--hc-prefilter-min-mapq", type=int, default=20)
    parser.add_argument("--hc-prefilter-min-baseq", type=int, default=10)
    parser.add_argument("--hc-prefilter-min-alt-count", type=int, default=1)
    parser.add_argument("--hc-prefilter-min-indel-count", type=int, default=1)
    parser.add_argument("--hc-prefilter-min-alt-fraction", type=float, default=0.0)
    parser.add_argument("--hc-prefilter-padding", type=int, default=150)
    parser.add_argument("--hc-prefilter-max-depth", type=int, default=100000)
    parser.add_argument("--hc-prefilter-exclude-supplementary", action="store_true")
    parser.add_argument("--hc-prefilter-empty-behavior", choices=("error", "input"), default="input")
    parser.add_argument("--skip-star-align", action="store_true", help="Skip STAR and require --aligned-bam.")
    parser.add_argument("--skip-mark-duplicates", action="store_true", help="Skip duplicate marking and require --dedup-bam for downstream steps.")
    parser.add_argument("--skip-split-n-cigar", action="store_true", help="Skip SplitNCigarReads and require --split-bam for downstream steps.")
    parser.add_argument("--skip-base-recalibrator", action="store_true", help="Skip BaseRecalibrator and require --recal-table if ApplyBQSR still runs.")
    parser.add_argument("--skip-apply-bqsr", action="store_true", help="Skip ApplyBQSR and require --recal-bam for downstream steps.")
    parser.add_argument("--skip-hc-prefilter", action="store_true", help="Skip candidate interval prefiltering.")
    parser.add_argument("--skip-haplotype-caller", action="store_true", help="Skip HaplotypeCaller and require --raw-vcf.")
    parser.add_argument("--skip-variant-filtration", action="store_true", help="Skip VariantFiltration and leave the raw VCF as the final output.")
    parser.add_argument("--skip-bqsr", action="store_true", help="Compatibility shortcut for skipping both BaseRecalibrator and ApplyBQSR and using the split BAM directly.")
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--force", action="store_true", help="Rerun steps even if outputs already exist.")
    return parser.parse_args(argv)


def python_step_cmd(module_name: str, *args: str | Path) -> list[str]:
    return [sys.executable, "-m", f"{PYTHON_STEP_MODULE_PREFIX}.{module_name}", *map(str, args)]


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    resolve_rust_defaults(args)
    apply_skip_shortcuts(args)
    normalize_backends(args)
    validate_args(args)

    run_label = build_run_label(args)
    outdir = args.outdir / run_label
    logs = outdir / "logs"
    refwork = outdir / "reference"
    hc_scatter_dir = refwork / "hc_scatter"
    hc_shard_vcf_dir = outdir / "hc_shards"
    logs.mkdir(parents=True, exist_ok=True)
    refwork.mkdir(parents=True, exist_ok=True)
    (outdir / "sambamba_tmp").mkdir(parents=True, exist_ok=True)
    timings = outdir / "timings.tsv"

    sample = args.sample
    aligned_bam_default = outdir / f"{sample}.Aligned.sortedByCoord.out.bam"
    dedup_bam_default = outdir / f"{sample}.dedup.bam"
    split_bam_default = outdir / f"{sample}.split.bam"
    recal_table_default = outdir / f"{sample}.recal.table"
    recal_bam_default = outdir / f"{sample}.recal.bam"
    raw_vcf_default = outdir / f"{sample}.raw.vcf.gz"
    filtered_vcf = outdir / f"{sample}.filtered.vcf.gz"
    markdup_metrics = outdir / f"{sample}.markdup.metrics.txt"
    exon_bed = refwork / "exons.bed"
    merged_exon_bed = refwork / "exons.merged.bed"
    raw_interval_list = refwork / "exons.interval_list"
    merged_interval_list = refwork / "exons.merged.interval_list"
    hc_prefilter_interval_list = refwork / "hc_candidate.interval_list"
    hc_prefilter_bed = refwork / "hc_candidate.bed"
    hc_prefilter_summary = outdir / f"{sample}.hc_prefilter.summary.tsv"
    dict_path = args.ref_dict

    serial_steps: list[Step] = []
    hc_single_step: Step | None = None
    hc_scatter_step: Step | None = None
    hc_merge_step: Step | None = None
    hc_shard_steps: list[Step] = []
    variant_filtration_step: Step | None = None

    need_interval_list = (args.hc_prefilter_backend == "rust" and not args.skip_hc_prefilter) or not args.skip_haplotype_caller
    if need_interval_list:
        dict_path = add_reference_steps(args, serial_steps, refwork, dict_path)
        add_interval_steps(args, serial_steps, exon_bed, merged_exon_bed, raw_interval_list, merged_interval_list, dict_path)

    aligned_bam = resolve_aligned_bam(args, serial_steps, outdir, aligned_bam_default)
    dedup_bam = resolve_dedup_bam(args, serial_steps, aligned_bam, dedup_bam_default, markdup_metrics, outdir)
    split_bam = resolve_split_bam(args, serial_steps, dedup_bam, split_bam_default)
    recal_table = args.recal_table or recal_table_default
    calling_bam = resolve_calling_bam(args, serial_steps, split_bam, recal_table, recal_bam_default)

    interval_list = merged_interval_list if args.interval_policy == "merged" else raw_interval_list
    if args.hc_prefilter_backend == "rust" and not args.skip_hc_prefilter:
        serial_steps.append(build_hc_prefilter_step(args, calling_bam, interval_list, hc_prefilter_interval_list, hc_prefilter_summary, hc_prefilter_bed))
        interval_list = hc_prefilter_interval_list

    pair_hmm_threads = args.hc_threads_per_shard
    if pair_hmm_threads is None:
        pair_hmm_threads = min(8, max(1, args.threads // args.hc_scatter_count))

    raw_vcf = args.raw_vcf or raw_vcf_default
    if not args.skip_haplotype_caller:
        if args.hc_scatter_count == 1:
            hc_single_step = Step(
                "haplotype_caller",
                build_haplotype_caller_step_command(args, calling_bam, interval_list, raw_vcf, pair_hmm_threads),
                (raw_vcf, Path(f"{raw_vcf}.tbi")),
                PYTHON_STEP_ENV,
            )
        else:
            hc_scatter_dir.mkdir(parents=True, exist_ok=True)
            hc_shard_vcf_dir.mkdir(parents=True, exist_ok=True)
            first_scatter_output = hc_scatter_dir / "0000-scattered.interval_list"
            hc_scatter_step = Step(
                "split_haplotype_caller_intervals",
                python_step_cmd(
                    "split_hc_intervals",
                    "--ref",
                    args.ref,
                    "--input-interval-list",
                    interval_list,
                    "--scatter-count",
                    str(args.hc_scatter_count),
                    "--output-dir",
                    hc_scatter_dir,
                    "--java-mem",
                    args.java_mem,
                    "--backend",
                    args.interval_backend,
                    "--rust-bin",
                    args.rust_interval_tools,
                ),
                (first_scatter_output,),
                PYTHON_STEP_ENV,
            )
            hc_merge_step = Step(
                "merge_haplotype_caller_vcfs",
                [],
                (raw_vcf, Path(f"{raw_vcf}.tbi")),
            )

    if not args.skip_variant_filtration:
        variant_filtration_step = Step(
            "variant_filtration",
            python_step_cmd(
                "variant_filtration",
                "--ref",
                args.ref,
                "--input-vcf",
                raw_vcf,
                "--output-vcf",
                filtered_vcf,
                "--java-mem",
                args.java_mem,
            ),
            (filtered_vcf, Path(f"{filtered_vcf}.tbi")),
            PYTHON_STEP_ENV,
        )

    for step in serial_steps:
        print(f"[run] {step.name}", flush=True)
        run_step(step, logs, timings, args.force)

    if hc_single_step is not None:
        print(f"[run] {hc_single_step.name}", flush=True)
        run_step(hc_single_step, logs, timings, args.force)
    elif hc_scatter_step is not None:
        print(f"[run] {hc_scatter_step.name}", flush=True)
        run_step(hc_scatter_step, logs, timings, args.force)

        shard_intervals = sorted(hc_scatter_dir.glob("*-scattered.interval_list"))
        if not shard_intervals:
            raise SystemExit(f"No scattered interval files found under {hc_scatter_dir}")

        for shard_index, shard_interval in enumerate(shard_intervals):
            shard_vcf = hc_shard_vcf_dir / f"{sample}.hc_shard_{shard_index:04d}.vcf.gz"
            shard_step = Step(
                f"haplotype_caller_shard_{shard_index:04d}",
                build_haplotype_caller_step_command(args, calling_bam, shard_interval, shard_vcf, pair_hmm_threads),
                (shard_vcf, Path(f"{shard_vcf}.tbi")),
                PYTHON_STEP_ENV,
            )
            hc_shard_steps.append(shard_step)

        for shard_step in hc_shard_steps:
            print(f"[run] {shard_step.name}", flush=True)
        run_parallel_steps(
            hc_shard_steps,
            logs,
            timings,
            args.force,
            min(args.hc_parallel_jobs, len(hc_shard_steps)),
        )

        if hc_merge_step is None:
            raise AssertionError("merge_haplotype_caller_vcfs step was not initialized")
        shard_vcfs = [step.outputs[0] for step in hc_shard_steps]
        hc_merge_step = Step(
            hc_merge_step.name,
            python_step_cmd(
                "merge_vcfs",
                *[item for vcf in shard_vcfs for item in ("--input-vcf", vcf)],
                "--output-vcf",
                raw_vcf,
                "--java-mem",
                args.java_mem,
            ),
            hc_merge_step.outputs,
            PYTHON_STEP_ENV,
        )
        print(f"[run] {hc_merge_step.name}", flush=True)
        run_step(hc_merge_step, logs, timings, args.force)

    if variant_filtration_step is not None:
        print(f"[run] {variant_filtration_step.name}", flush=True)
        run_step(variant_filtration_step, logs, timings, args.force)

    final_vcf = filtered_vcf if variant_filtration_step is not None else raw_vcf
    print(f"timings: {timings}")
    print(f"vcf: {final_vcf}")
    return 0


def resolve_rust_defaults(args: argparse.Namespace) -> None:
    args.rust_interval_tools = args.rust_interval_tools or default_rust_binary(args.rust_bin_dir, "rust_interval_tools")
    args.rust_split_n_cigar = args.rust_split_n_cigar or default_rust_binary(args.rust_bin_dir, "rust_split_n_cigar_reads")
    args.rust_base_recalibrator = args.rust_base_recalibrator or default_rust_binary(args.rust_bin_dir, "rust_base_recalibrator")
    args.rust_apply_bqsr = args.rust_apply_bqsr or default_rust_binary(args.rust_bin_dir, "rust_apply_bqsr")
    args.rust_hc_prefilter = args.rust_hc_prefilter or default_rust_binary(args.rust_bin_dir, "rust_hc_prefilter")
    args.rust_haplotype_caller = args.rust_haplotype_caller or default_rust_binary(args.rust_bin_dir, "rust_haplotype_caller")
    args.rust_mark_duplicates = default_rust_binary(args.rust_bin_dir, "rust_mark_duplicates")


def apply_skip_shortcuts(args: argparse.Namespace) -> None:
    if args.skip_bqsr:
        if args.recal_table is not None or args.recal_bam is not None:
            raise SystemExit("--skip-bqsr cannot be combined with --recal-table or --recal-bam.")
        args.skip_base_recalibrator = True
        args.skip_apply_bqsr = True


def normalize_backends(args: argparse.Namespace) -> None:
    if args.no_rust:
        args.interval_backend = "gatk"
        args.split_n_cigar_backend = "gatk"
        args.base_recalibrator_backend = "gatk"
        args.apply_bqsr_backend = "gatk"
        args.hc_prefilter_backend = "none"
        args.hc_backend = "gatk"
        if args.mode == "auto":
            args.mode = "baseline"
        return

    if args.interval_backend == "auto":
        args.interval_backend = "rust" if args.rust_interval_tools.exists() else "gatk"
    if args.split_n_cigar_backend == "auto":
        args.split_n_cigar_backend = "rust" if args.rust_split_n_cigar.exists() else "gatk"
    if args.base_recalibrator_backend == "auto":
        args.base_recalibrator_backend = "rust" if args.rust_base_recalibrator.exists() else "gatk"
    if args.apply_bqsr_backend == "auto":
        args.apply_bqsr_backend = "rust" if args.rust_apply_bqsr.exists() else "gatk"
    if args.hc_prefilter_backend == "auto":
        args.hc_prefilter_backend = "rust" if args.rust_hc_prefilter.exists() else "none"
    if args.mode == "auto":
        args.mode = "rust" if args.rust_mark_duplicates.exists() else "baseline"

    if args.skip_hc_prefilter:
        args.hc_prefilter_backend = "none"


def validate_args(args: argparse.Namespace) -> None:
    if args.known_sites is None:
        args.known_sites = [path for path in DEFAULT_KNOWN_SITES if path.exists()]

    if args.dedup_bam is not None and not args.skip_mark_duplicates:
        raise SystemExit("--dedup-bam is only valid together with --skip-mark-duplicates.")
    if args.split_bam is not None and not args.skip_split_n_cigar:
        raise SystemExit("--split-bam is only valid together with --skip-split-n-cigar.")
    if args.recal_table is not None and not args.skip_base_recalibrator:
        raise SystemExit("--recal-table is only valid together with --skip-base-recalibrator.")
    if args.recal_bam is not None and not args.skip_apply_bqsr:
        raise SystemExit("--recal-bam is only valid together with --skip-apply-bqsr.")
    if args.raw_vcf is not None and not args.skip_haplotype_caller:
        raise SystemExit("--raw-vcf is only valid together with --skip-haplotype-caller.")
    if args.skip_star_align and args.aligned_bam is None:
        raise SystemExit("--skip-star-align requires --aligned-bam.")
    if args.skip_mark_duplicates and not args.skip_split_n_cigar and args.dedup_bam is None:
        raise SystemExit("--skip-mark-duplicates requires --dedup-bam unless SplitNCigarReads is also skipped.")
    if need_split_bam(args) and args.skip_split_n_cigar and args.split_bam is None:
        raise SystemExit("--skip-split-n-cigar requires --split-bam for the remaining downstream steps.")
    if args.skip_base_recalibrator and not args.skip_apply_bqsr and args.recal_table is None:
        raise SystemExit("--skip-base-recalibrator requires --recal-table when ApplyBQSR still runs.")
    if need_calling_bam(args) and args.skip_apply_bqsr and not args.skip_bqsr and args.recal_bam is None:
        raise SystemExit("--skip-apply-bqsr requires --recal-bam for the remaining downstream steps.")
    if args.skip_haplotype_caller and args.raw_vcf is None:
        raise SystemExit("--skip-haplotype-caller requires --raw-vcf.")
    if not args.skip_base_recalibrator and not args.known_sites:
        raise SystemExit("BaseRecalibrator needs at least one --known-sites VCF.")
    if args.hc_scatter_count < 1:
        raise SystemExit("--hc-scatter-count must be at least 1.")
    if args.hc_threads_per_shard is not None and args.hc_threads_per_shard < 1:
        raise SystemExit("--hc-threads-per-shard must be at least 1.")
    if args.hc_parallel_jobs < 1:
        raise SystemExit("--hc-parallel-jobs must be at least 1.")
    if args.base_recalibrator_region_bases < 1:
        raise SystemExit("--base-recalibrator-region-bases must be at least 1.")
    if args.hc_prefilter_min_mapq < 0 or args.hc_prefilter_min_baseq < 0:
        raise SystemExit("--hc-prefilter-min-mapq and --hc-prefilter-min-baseq must be non-negative.")
    if args.hc_prefilter_min_alt_count < 1 or args.hc_prefilter_min_indel_count < 1:
        raise SystemExit("--hc-prefilter-min-alt-count and --hc-prefilter-min-indel-count must be at least 1.")
    if not 0.0 <= args.hc_prefilter_min_alt_fraction <= 1.0:
        raise SystemExit("--hc-prefilter-min-alt-fraction must be between 0 and 1.")
    if args.hc_prefilter_padding < 0:
        raise SystemExit("--hc-prefilter-padding must be non-negative.")
    if args.hc_prefilter_max_depth < 1:
        raise SystemExit("--hc-prefilter-max-depth must be at least 1.")
    if args.hc_memory_gb < 1:
        raise SystemExit("--hc-memory-gb must be at least 1.")
    if args.hc_backend == "rust" and not args.rust_haplotype_caller.exists():
        raise SystemExit(f"rust_haplotype_caller binary not found: {args.rust_haplotype_caller}")


def build_run_label(args: argparse.Namespace) -> str:
    run_label = args.mode
    if args.interval_policy != "raw" or (not args.skip_haplotype_caller and args.hc_scatter_count != 1):
        run_label = f"{args.mode}.{args.interval_policy}.hcscatter{args.hc_scatter_count}"
    if args.interval_backend != "gatk":
        run_label = f"{run_label}.intervals-{args.interval_backend}"
    if args.split_n_cigar_backend != "gatk":
        run_label = f"{run_label}.splitncigar-{args.split_n_cigar_backend}"
    if args.base_recalibrator_backend != "gatk":
        run_label = f"{run_label}.baserecal-{args.base_recalibrator_backend}"
    if args.apply_bqsr_backend != "gatk":
        run_label = f"{run_label}.applybqsr-{args.apply_bqsr_backend}"
    if args.hc_prefilter_backend != "none":
        run_label = f"{run_label}.hcpf-{args.hc_prefilter_backend}"
    if args.hc_backend != "gatk":
        run_label = f"{run_label}.hc-{args.hc_backend}"

    skipped = []
    if args.aligned_bam is not None or args.skip_star_align:
        skipped.append("star")
    if args.skip_mark_duplicates:
        skipped.append("markdup")
    if args.skip_split_n_cigar:
        skipped.append("splitncigar")
    if args.skip_base_recalibrator:
        skipped.append("baserecal")
    if args.skip_apply_bqsr:
        skipped.append("applybqsr")
    if args.skip_hc_prefilter:
        skipped.append("hcpf")
    if args.skip_haplotype_caller:
        skipped.append("hc")
    if args.skip_variant_filtration:
        skipped.append("filter")
    if skipped:
        run_label = f"{run_label}.skip-{'-'.join(skipped)}"
    return run_label


def add_reference_steps(args: argparse.Namespace, steps: list[Step], refwork: Path, dict_path: Path) -> Path:
    fai_path = Path(f"{args.ref}.fai")
    if not fai_path.exists():
        steps.append(Step("ref_faidx", [str(SAMTOOLS), "faidx", str(args.ref)], (fai_path,)))

    if dict_path.exists():
        return dict_path

    dict_path = refwork / f"{args.ref.stem}.dict"
    steps.append(
        Step(
            "ref_dict",
            gatk_cmd(args.java_mem, "CreateSequenceDictionary", "-R", args.ref, "-O", dict_path),
            (dict_path,),
        )
    )
    return dict_path


def add_interval_steps(
    args: argparse.Namespace,
    steps: list[Step],
    exon_bed: Path,
    merged_exon_bed: Path,
    raw_interval_list: Path,
    merged_interval_list: Path,
    dict_path: Path,
) -> None:
    steps.append(
        Step(
            "build_exon_bed",
            python_step_cmd(
                "build_exon_bed",
                "--gtf",
                args.gtf,
                "--ref",
                args.ref,
                "--output-bed",
                exon_bed,
            ),
            (exon_bed,),
            PYTHON_STEP_ENV,
        )
    )
    steps.append(
        Step(
            "exon_interval_list",
            python_step_cmd(
                "bed_to_interval_list",
                "--input-bed",
                exon_bed,
                "--ref-dict",
                dict_path,
                "--java-mem",
                args.java_mem,
                "--output-interval-list",
                raw_interval_list,
                "--backend",
                args.interval_backend,
                "--rust-bin",
                args.rust_interval_tools,
            ),
            (raw_interval_list,),
            PYTHON_STEP_ENV,
        )
    )

    if args.interval_policy != "merged":
        return

    steps.append(
        Step(
            "build_merged_exon_bed",
            python_step_cmd(
                "build_merged_exon_bed",
                "--input-bed",
                exon_bed,
                "--output-bed",
                merged_exon_bed,
            ),
            (merged_exon_bed,),
            PYTHON_STEP_ENV,
        )
    )
    steps.append(
        Step(
            "merged_exon_interval_list",
            python_step_cmd(
                "bed_to_interval_list",
                "--input-bed",
                merged_exon_bed,
                "--ref-dict",
                dict_path,
                "--java-mem",
                args.java_mem,
                "--output-interval-list",
                merged_interval_list,
                "--backend",
                args.interval_backend,
                "--rust-bin",
                args.rust_interval_tools,
            ),
            (merged_interval_list,),
            PYTHON_STEP_ENV,
        )
    )


def resolve_aligned_bam(args: argparse.Namespace, steps: list[Step], outdir: Path, aligned_bam_default: Path) -> Path:
    if args.aligned_bam is not None:
        return args.aligned_bam
    if args.skip_star_align:
        raise SystemExit("--skip-star-align requires --aligned-bam.")
    if args.skip_mark_duplicates:
        return aligned_bam_default
    steps.append(
        Step(
            "star_align",
            python_step_cmd(
                "star_align",
                "--sample",
                args.sample,
                "--fastq1",
                args.fastq1,
                "--fastq2",
                args.fastq2,
                "--star-index",
                args.star_index,
                "--threads",
                str(args.threads),
                "--output-prefix",
                outdir / f"{args.sample}.",
            ),
            (aligned_bam_default,),
            PYTHON_STEP_ENV,
        )
    )
    return aligned_bam_default


def resolve_dedup_bam(
    args: argparse.Namespace,
    steps: list[Step],
    aligned_bam: Path,
    dedup_bam_default: Path,
    markdup_metrics: Path,
    outdir: Path,
) -> Path:
    if args.skip_mark_duplicates:
        return args.dedup_bam or dedup_bam_default

    step_name = {
        "baseline": "mark_duplicates_picard",
        "sambamba": "mark_duplicates_sambamba",
    }.get(args.mode, "mark_duplicates_rust")
    steps.append(
        Step(
            step_name,
            python_step_cmd(
                "mark_duplicates",
                "--mode",
                args.mode,
                "--input-bam",
                aligned_bam,
                "--output-bam",
                dedup_bam_default,
                "--output-metrics",
                markdup_metrics,
                "--threads",
                str(args.threads),
                "--java-mem",
                args.java_mem,
                "--tmpdir",
                outdir / "sambamba_tmp",
                "--rust-bin",
                args.rust_mark_duplicates,
            ),
            (dedup_bam_default, bam_index_path(dedup_bam_default)),
            PYTHON_STEP_ENV,
        )
    )
    return dedup_bam_default


def resolve_split_bam(args: argparse.Namespace, steps: list[Step], dedup_bam: Path, split_bam_default: Path) -> Path:
    if args.skip_split_n_cigar:
        return args.split_bam or split_bam_default

    steps.append(
        Step(
            "split_n_cigar_reads",
            python_step_cmd(
                "split_n_cigar",
                "--ref",
                args.ref,
                "--input-bam",
                dedup_bam,
                "--output-bam",
                split_bam_default,
                "--java-mem",
                args.java_mem,
                "--backend",
                args.split_n_cigar_backend,
                "--rust-bin",
                args.rust_split_n_cigar,
                "--threads",
                str(args.threads),
                "--rust-mode",
                args.split_n_cigar_rust_mode,
            ),
            (split_bam_default, bam_index_path(split_bam_default)),
            PYTHON_STEP_ENV,
        )
    )
    return split_bam_default


def resolve_calling_bam(
    args: argparse.Namespace,
    steps: list[Step],
    split_bam: Path,
    recal_table: Path,
    recal_bam_default: Path,
) -> Path:
    if args.skip_bqsr:
        return split_bam

    if not args.skip_base_recalibrator:
        bqsr_cmd = python_step_cmd(
            "base_recalibrator",
            "--ref",
            args.ref,
            "--input-bam",
            split_bam,
            "--output-table",
            recal_table,
            "--java-mem",
            args.java_mem,
            "--backend",
            args.base_recalibrator_backend,
            "--rust-bin",
            args.rust_base_recalibrator,
            "--threads",
            str(args.threads),
            "--region-bases",
            str(args.base_recalibrator_region_bases),
        )
        for known in args.known_sites:
            bqsr_cmd.extend(["--known-sites", str(known)])
        steps.append(Step("base_recalibrator", bqsr_cmd, (recal_table,), PYTHON_STEP_ENV))

    if args.skip_apply_bqsr:
        return args.recal_bam or recal_bam_default

    steps.append(
        Step(
            "apply_bqsr",
            python_step_cmd(
                "apply_bqsr",
                "--ref",
                args.ref,
                "--input-bam",
                split_bam,
                "--input-table",
                recal_table,
                "--output-bam",
                recal_bam_default,
                "--java-mem",
                args.java_mem,
                "--backend",
                args.apply_bqsr_backend,
                "--rust-bin",
                args.rust_apply_bqsr,
                "--threads",
                str(args.threads),
            ),
            (recal_bam_default, bam_index_path(recal_bam_default)),
            PYTHON_STEP_ENV,
        )
    )
    return recal_bam_default


def build_hc_prefilter_step(
    args: argparse.Namespace,
    calling_bam: Path,
    interval_list: Path,
    output_interval_list: Path,
    output_summary: Path,
    output_bed: Path,
) -> Step:
    command = python_step_cmd(
        "hc_prefilter",
        "--ref",
        args.ref,
        "--input-bam",
        calling_bam,
        "--input-interval-list",
        interval_list,
        "--output-interval-list",
        output_interval_list,
        "--output-summary",
        output_summary,
        "--output-bed",
        output_bed,
        "--min-mapq",
        str(args.hc_prefilter_min_mapq),
        "--min-baseq",
        str(args.hc_prefilter_min_baseq),
        "--min-alt-count",
        str(args.hc_prefilter_min_alt_count),
        "--min-indel-count",
        str(args.hc_prefilter_min_indel_count),
        "--min-alt-fraction",
        str(args.hc_prefilter_min_alt_fraction),
        "--padding",
        str(args.hc_prefilter_padding),
        "--max-depth",
        str(args.hc_prefilter_max_depth),
        "--threads",
        str(args.threads),
        "--empty-behavior",
        args.hc_prefilter_empty_behavior,
        "--rust-bin",
        args.rust_hc_prefilter,
    )
    if args.hc_prefilter_exclude_supplementary:
        command.append("--exclude-supplementary")
    return Step(
        "hc_candidate_prefilter",
        command,
        (output_interval_list, output_summary, output_bed),
        PYTHON_STEP_ENV,
    )


def build_haplotype_caller_step_command(
    args: argparse.Namespace,
    calling_bam: Path,
    interval_list: Path,
    output_vcf: Path,
    pair_hmm_threads: int,
) -> list[str]:
    return python_step_cmd(
        "haplotype_caller",
        "--ref",
        args.ref,
        "--input-bam",
        calling_bam,
        "--input-interval-list",
        interval_list,
        "--output-vcf",
        output_vcf,
        "--backend",
        args.hc_backend,
        "--rust-bin",
        args.rust_haplotype_caller,
        "--threads",
        str(args.threads),
        "--memory-gb",
        str(args.hc_memory_gb),
        "--pair-hmm-threads",
        str(pair_hmm_threads),
        "--pair-hmm-implementation",
        args.hc_pair_hmm_implementation,
        "--java-mem",
        args.java_mem,
        "--dbsnp",
        args.dbsnp,
    )


def need_split_bam(args: argparse.Namespace) -> bool:
    if args.skip_bqsr:
        return need_calling_bam(args)
    if not args.skip_base_recalibrator:
        return True
    if not args.skip_apply_bqsr:
        return True
    return False


def need_calling_bam(args: argparse.Namespace) -> bool:
    return (args.hc_prefilter_backend == "rust" and not args.skip_hc_prefilter) or not args.skip_haplotype_caller


if __name__ == "__main__":
    raise SystemExit(main())

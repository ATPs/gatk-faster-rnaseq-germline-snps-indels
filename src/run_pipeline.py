#!/usr/bin/env python3
"""Run and time an RNA-seq germline SNP/indel pipeline.

The runner keeps each pipeline command explicit and writes per-step logs plus a
timings.tsv file, so baseline and optimized runs can be compared directly.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


DEFAULT_FASTQ_DIR = Path("/data1/xlab/researches/AML/Leucegene/raw/PRJNA214592/SRR949115")
DEFAULT_GATK_RESOURCES = Path("/data1/pub/gatk/broad_hg38")
DEFAULT_REF = DEFAULT_GATK_RESOURCES / "Homo_sapiens_assembly38.fasta"
DEFAULT_REF_DICT = DEFAULT_GATK_RESOURCES / "Homo_sapiens_assembly38.dict"
DEFAULT_GTF = Path("/data1/pub/genome/Human/humanGENCODE/gencode.v46.annotation.gtf")
DEFAULT_STAR_INDEX = Path("/XCLabServer002_fastIO/gatk-faster-rnaseq/reference/broad_hg38_gencode_v46_readlen100_STAR.index")
DEFAULT_DBSNP = DEFAULT_GATK_RESOURCES / "Homo_sapiens_assembly38.dbsnp138.vcf.gz"
DEFAULT_KNOWN_SITES = [
    DEFAULT_GATK_RESOURCES / "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
    DEFAULT_GATK_RESOURCES / "Homo_sapiens_assembly38.known_indels.vcf.gz",
    DEFAULT_GATK_RESOURCES / "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
]

GATK = Path("/data/p/gatk/gatk-4.6.2.0/gatk")
SAMTOOLS = Path("/data/p/samtools/samtools-1.22.1_installed/bin/samtools")
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RUST_BIN_DIR = SCRIPT_DIR / "rust" / "bin"


@dataclass(frozen=True)
class Step:
    name: str
    command: list[str]
    outputs: tuple[Path, ...]


def quote_cmd(command: list[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def run_step(step: Step, log_dir: Path, timings_path: Path, force: bool) -> None:
    if step.outputs and not force and all(path.exists() for path in step.outputs):
        append_timing(timings_path, step.name, "skipped", 0.0, quote_cmd(step.command), "")
        return

    log_path = log_dir / f"{step.name}.log"
    start = time.perf_counter()
    status = "ok"
    message = ""
    with log_path.open("w") as log:
        log.write(f"$ {quote_cmd(step.command)}\n\n")
        log.flush()
        try:
            subprocess.run(step.command, check=True, stdout=log, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as exc:
            status = "failed"
            message = f"exit_code={exc.returncode}; log={log_path}"
            raise
        finally:
            elapsed = time.perf_counter() - start
            append_timing(timings_path, step.name, status, elapsed, quote_cmd(step.command), message)


def run_parallel_steps(steps: list[Step], log_dir: Path, timings_path: Path, force: bool, max_parallel: int) -> None:
    if max_parallel < 1:
        raise ValueError("max_parallel must be at least 1")

    for batch_start in range(0, len(steps), max_parallel):
        running: list[tuple[Step, subprocess.Popen[bytes], object, Path, float]] = []
        first_error: subprocess.CalledProcessError | None = None

        for step in steps[batch_start: batch_start + max_parallel]:
            if step.outputs and not force and all(path.exists() for path in step.outputs):
                append_timing(timings_path, step.name, "skipped", 0.0, quote_cmd(step.command), "")
                continue

            log_path = log_dir / f"{step.name}.log"
            log = log_path.open("w")
            log.write(f"$ {quote_cmd(step.command)}\n\n")
            log.flush()
            start = time.perf_counter()
            process = subprocess.Popen(step.command, stdout=log, stderr=subprocess.STDOUT)
            running.append((step, process, log, log_path, start))

        for step, process, log, log_path, start in running:
            status = "ok"
            message = ""
            return_code = process.wait()
            log.close()
            if return_code != 0:
                status = "failed"
                message = f"exit_code={return_code}; log={log_path}"
                if first_error is None:
                    first_error = subprocess.CalledProcessError(return_code, step.command)
            elapsed = time.perf_counter() - start
            append_timing(timings_path, step.name, status, elapsed, quote_cmd(step.command), message)

        if first_error is not None:
            raise first_error


def append_timing(path: Path, step: str, status: str, seconds: float, command: str, message: str) -> None:
    new_file = not path.exists()
    with path.open("a") as handle:
        if new_file:
            handle.write("step\tstatus\tseconds\tminutes\tcommand\tmessage\n")
        handle.write(f"{step}\t{status}\t{seconds:.3f}\t{seconds / 60:.3f}\t{command}\t{message}\n")


def bam_index_path(bam: Path) -> Path:
    if bam.name.endswith(".bam"):
        return bam.with_suffix(".bai")
    return Path(f"{bam}.bai")


def existing(path: Path) -> Path:
    path = Path(path)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"not found: {path}")
    return path


def default_rust_binary(bin_dir: Path, name: str) -> Path:
    return bin_dir / name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a timed GATK RNA-seq germline SNP/indel pipeline."
    )
    parser.add_argument("--sample", default="SRR949115")
    parser.add_argument("--fastq1", type=existing, default=DEFAULT_FASTQ_DIR / "SRR949115_1.fastq.gz")
    parser.add_argument("--fastq2", type=existing, default=DEFAULT_FASTQ_DIR / "SRR949115_2.fastq.gz")
    parser.add_argument("--aligned-bam", type=existing, help="Use an existing STAR coordinate-sorted BAM.")
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
    parser.add_argument("--skip-bqsr", action="store_true", help="Allow a run without known-sites resources.")
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--force", action="store_true", help="Rerun steps even if outputs already exist.")
    return parser.parse_args()


def gatk_cmd(java_mem: str, *args: str | Path) -> list[str]:
    return [str(GATK), "--java-options", f"-Xmx{java_mem}", *map(str, args)]


def python_step_cmd(script_name: str, *args: str | Path) -> list[str]:
    return [sys.executable, str(SCRIPT_DIR / script_name), *map(str, args)]


def main() -> int:
    args = parse_args()
    args.rust_interval_tools = args.rust_interval_tools or default_rust_binary(args.rust_bin_dir, "rust_interval_tools")
    args.rust_split_n_cigar = args.rust_split_n_cigar or default_rust_binary(args.rust_bin_dir, "rust_split_n_cigar_reads")
    args.rust_base_recalibrator = args.rust_base_recalibrator or default_rust_binary(args.rust_bin_dir, "rust_base_recalibrator")
    args.rust_apply_bqsr = args.rust_apply_bqsr or default_rust_binary(args.rust_bin_dir, "rust_apply_bqsr")
    args.rust_hc_prefilter = args.rust_hc_prefilter or default_rust_binary(args.rust_bin_dir, "rust_hc_prefilter")
    rust_mark_duplicates = default_rust_binary(args.rust_bin_dir, "rust_mark_duplicates")

    if args.no_rust:
        args.interval_backend = "gatk"
        args.split_n_cigar_backend = "gatk"
        args.base_recalibrator_backend = "gatk"
        args.apply_bqsr_backend = "gatk"
        args.hc_prefilter_backend = "none"
        if args.mode == "auto":
            args.mode = "baseline"
    else:
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
            args.mode = "rust" if rust_mark_duplicates.exists() else "baseline"

    if args.known_sites is None:
        args.known_sites = [path for path in DEFAULT_KNOWN_SITES if path.exists()]
    if not args.skip_bqsr and not args.known_sites:
        raise SystemExit("BQSR needs at least one --known-sites VCF. Use --skip-bqsr for a timing smoke run.")
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

    run_label = args.mode
    if args.interval_policy != "raw" or args.hc_scatter_count != 1:
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
    outdir = args.outdir / run_label
    logs = outdir / "logs"
    refwork = outdir / "reference"
    hc_scatter_dir = refwork / "hc_scatter"
    hc_shard_vcf_dir = outdir / "hc_shards"
    logs.mkdir(parents=True, exist_ok=True)
    refwork.mkdir(parents=True, exist_ok=True)
    (outdir / "sambamba_tmp").mkdir(parents=True, exist_ok=True)
    hc_scatter_dir.mkdir(parents=True, exist_ok=True)
    hc_shard_vcf_dir.mkdir(parents=True, exist_ok=True)
    timings = outdir / "timings.tsv"

    sample = args.sample
    aligned_bam = outdir / f"{sample}.Aligned.sortedByCoord.out.bam"
    dedup_bam = outdir / f"{sample}.dedup.bam"
    split_bam = outdir / f"{sample}.split.bam"
    recal_table = outdir / f"{sample}.recal.table"
    calling_bam = outdir / f"{sample}.recal.bam"
    raw_vcf = outdir / f"{sample}.raw.vcf.gz"
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

    fai_path = Path(f"{args.ref}.fai")
    steps: list[Step] = []

    if not fai_path.exists():
        steps.append(Step("ref_faidx", [str(SAMTOOLS), "faidx", str(args.ref)], (fai_path,)))

    if not dict_path.exists():
        dict_path = refwork / f"{args.ref.stem}.dict"
        steps.append(
            Step(
                "ref_dict",
                gatk_cmd(args.java_mem, "CreateSequenceDictionary", "-R", args.ref, "-O", dict_path),
                (dict_path,),
            )
        )

    steps.append(
        Step(
            "build_exon_bed",
            python_step_cmd(
                "step_build_exon_bed.py",
                "--gtf",
                args.gtf,
                "--ref",
                args.ref,
                "--output-bed",
                exon_bed,
            ),
            (exon_bed,),
        )
    )

    steps.append(
        Step(
            "exon_interval_list",
            python_step_cmd(
                "step_bed_to_interval_list.py",
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
        )
    )

    if args.interval_policy == "merged":
        steps.append(
            Step(
                "build_merged_exon_bed",
                python_step_cmd(
                    "step_build_merged_exon_bed.py",
                    "--input-bed",
                    exon_bed,
                    "--output-bed",
                    merged_exon_bed,
                ),
                (merged_exon_bed,),
            )
        )
        steps.append(
            Step(
                "merged_exon_interval_list",
                python_step_cmd(
                    "step_bed_to_interval_list.py",
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
            )
        )

    if args.aligned_bam:
        aligned_bam = args.aligned_bam
    else:
        steps.append(
            Step(
                "star_align",
                python_step_cmd(
                    "step_star_align.py",
                    "--sample",
                    sample,
                    "--fastq1",
                    args.fastq1,
                    "--fastq2",
                    args.fastq2,
                    "--star-index",
                    args.star_index,
                    "--threads",
                    str(args.threads),
                    "--output-prefix",
                    outdir / f"{sample}.",
                ),
                (aligned_bam,),
            )
        )

    if args.mode == "baseline":
        steps.append(
            Step(
                "mark_duplicates_picard",
                python_step_cmd(
                    "step_mark_duplicates.py",
                    "--mode",
                    args.mode,
                    "--input-bam",
                    aligned_bam,
                    "--output-bam",
                    dedup_bam,
                    "--output-metrics",
                    markdup_metrics,
                    "--threads",
                    str(args.threads),
                    "--java-mem",
                    args.java_mem,
                    "--tmpdir",
                    outdir / "sambamba_tmp",
                    "--rust-bin",
                    rust_mark_duplicates,
                ),
                (dedup_bam, bam_index_path(dedup_bam)),
            )
        )
    elif args.mode == "sambamba":
        steps.append(
            Step(
                "mark_duplicates_sambamba",
                python_step_cmd(
                    "step_mark_duplicates.py",
                    "--mode",
                    args.mode,
                    "--input-bam",
                    aligned_bam,
                    "--output-bam",
                    dedup_bam,
                    "--output-metrics",
                    markdup_metrics,
                    "--threads",
                    str(args.threads),
                    "--java-mem",
                    args.java_mem,
                    "--tmpdir",
                    outdir / "sambamba_tmp",
                    "--rust-bin",
                    rust_mark_duplicates,
                ),
                (dedup_bam, bam_index_path(dedup_bam)),
            )
        )
    else:
        steps.append(
            Step(
                "mark_duplicates_rust",
                python_step_cmd(
                    "step_mark_duplicates.py",
                    "--mode",
                    args.mode,
                    "--input-bam",
                    aligned_bam,
                    "--output-bam",
                    dedup_bam,
                    "--output-metrics",
                    markdup_metrics,
                    "--threads",
                    str(args.threads),
                    "--java-mem",
                    args.java_mem,
                    "--tmpdir",
                    outdir / "sambamba_tmp",
                    "--rust-bin",
                    rust_mark_duplicates,
                ),
                (dedup_bam, bam_index_path(dedup_bam)),
            )
        )

    steps.append(
        Step(
            "split_n_cigar_reads",
            python_step_cmd(
                "step_split_n_cigar.py",
                "--ref",
                args.ref,
                "--input-bam",
                dedup_bam,
                "--output-bam",
                split_bam,
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
            (split_bam, bam_index_path(split_bam)),
        )
    )

    if args.skip_bqsr:
        calling_bam = split_bam
    else:
        bqsr_cmd = python_step_cmd(
            "step_base_recalibrator.py",
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
        steps.append(Step("base_recalibrator", bqsr_cmd, (recal_table,)))
        steps.append(
            Step(
                "apply_bqsr",
                python_step_cmd(
                    "step_apply_bqsr.py",
                    "--ref",
                    args.ref,
                    "--input-bam",
                    split_bam,
                    "--input-table",
                    recal_table,
                    "--output-bam",
                    calling_bam,
                    "--java-mem",
                    args.java_mem,
                    "--backend",
                    args.apply_bqsr_backend,
                    "--rust-bin",
                    args.rust_apply_bqsr,
                    "--threads",
                    str(args.threads),
                ),
                (calling_bam, bam_index_path(calling_bam)),
            )
        )

    interval_list = merged_interval_list if args.interval_policy == "merged" else raw_interval_list
    if args.hc_prefilter_backend == "rust":
        hc_prefilter_cmd = python_step_cmd(
            "step_hc_prefilter.py",
            "--ref",
            args.ref,
            "--input-bam",
            calling_bam,
            "--input-interval-list",
            interval_list,
            "--output-interval-list",
            hc_prefilter_interval_list,
            "--output-summary",
            hc_prefilter_summary,
            "--output-bed",
            hc_prefilter_bed,
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
            hc_prefilter_cmd.append("--exclude-supplementary")
        steps.append(
            Step(
                "hc_candidate_prefilter",
                hc_prefilter_cmd,
                (hc_prefilter_interval_list, hc_prefilter_summary, hc_prefilter_bed),
            )
        )
        interval_list = hc_prefilter_interval_list

    pair_hmm_threads = args.hc_threads_per_shard
    if pair_hmm_threads is None:
        pair_hmm_threads = min(8, max(1, args.threads // args.hc_scatter_count))

    hc_scatter_step: Step | None = None
    hc_shard_steps: list[Step] = []
    hc_merge_step: Step | None = None

    if args.hc_scatter_count == 1:
        hc_cmd = python_step_cmd(
            "step_haplotype_caller.py",
            "--ref",
            args.ref,
            "--input-bam",
            calling_bam,
            "--input-interval-list",
            interval_list,
            "--output-vcf",
            raw_vcf,
            "--pair-hmm-threads",
            str(pair_hmm_threads),
            "--java-mem",
            args.java_mem,
            "--dbsnp",
            args.dbsnp,
        )
        steps.append(Step("haplotype_caller", hc_cmd, (raw_vcf, Path(f"{raw_vcf}.tbi"))))
    else:
        first_scatter_output = hc_scatter_dir / "0000-scattered.interval_list"
        hc_scatter_step = Step(
            "split_haplotype_caller_intervals",
            python_step_cmd(
                "step_split_hc_intervals.py",
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
        )
        hc_merge_step = Step(
            "merge_haplotype_caller_vcfs",
            [],
            (raw_vcf, Path(f"{raw_vcf}.tbi")),
        )

    steps.append(
        Step(
            "variant_filtration",
            python_step_cmd(
                "step_variant_filtration.py",
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
        )
    )

    for step in steps:
        if step.name == "variant_filtration" and hc_scatter_step is not None:
            break
        print(f"[run] {step.name}", flush=True)
        run_step(step, logs, timings, args.force)

    if hc_scatter_step is not None:
        print(f"[run] {hc_scatter_step.name}", flush=True)
        run_step(hc_scatter_step, logs, timings, args.force)

        shard_intervals = sorted(hc_scatter_dir.glob("*-scattered.interval_list"))
        if not shard_intervals:
            raise SystemExit(f"No scattered interval files found under {hc_scatter_dir}")

        for shard_index, shard_interval in enumerate(shard_intervals):
            shard_vcf = hc_shard_vcf_dir / f"{sample}.hc_shard_{shard_index:04d}.vcf.gz"
            hc_shard_steps.append(
                Step(
                    f"haplotype_caller_shard_{shard_index:04d}",
                    python_step_cmd(
                        "step_haplotype_caller.py",
                        "--ref",
                        args.ref,
                        "--input-bam",
                        calling_bam,
                        "--input-interval-list",
                        shard_interval,
                        "--output-vcf",
                        shard_vcf,
                        "--pair-hmm-threads",
                        str(pair_hmm_threads),
                        "--java-mem",
                        args.java_mem,
                        "--dbsnp",
                        args.dbsnp,
                    ),
                    (shard_vcf, Path(f"{shard_vcf}.tbi")),
                )
            )

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
                "step_merge_vcfs.py",
                *[item for vcf in shard_vcfs for item in ("--input-vcf", vcf)],
                "--output-vcf",
                raw_vcf,
                "--java-mem",
                args.java_mem,
            ),
            hc_merge_step.outputs,
        )
        print(f"[run] {hc_merge_step.name}", flush=True)
        run_step(hc_merge_step, logs, timings, args.force)

        variant_filtration_step = steps[-1]
        print(f"[run] {variant_filtration_step.name}", flush=True)
        run_step(variant_filtration_step, logs, timings, args.force)

    print(f"timings: {timings}")
    print(f"vcf: {filtered_vcf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Shared helpers for independently runnable pipeline step scripts."""

from __future__ import annotations

import argparse
import gzip
import shlex
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
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
STAR = Path("/data/p/star/STAR_2.7.11b/Linux_x86_64_static/STAR")
SAMTOOLS = Path("/data/p/samtools/samtools-1.22.1_installed/bin/samtools")
SAMBAMBA = Path("/data/p/tools/sambamba/bin/sambamba-1.0.1")
RUST_BIN_DIR = SCRIPT_DIR / "rust" / "bin"
RUST_INTERVAL_TOOLS = RUST_BIN_DIR / "rust_interval_tools"
RUST_SPLIT_N_CIGAR_READS = RUST_BIN_DIR / "rust_split_n_cigar_reads"
RUST_APPLY_BQSR = RUST_BIN_DIR / "rust_apply_bqsr"
RUST_BASE_RECALIBRATOR = RUST_BIN_DIR / "rust_base_recalibrator"
RUST_MARK_DUPLICATES = RUST_BIN_DIR / "rust_mark_duplicates"
RUST_HC_PREFILTER = RUST_BIN_DIR / "rust_hc_prefilter"


def quote_cmd(command: list[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def run_command(command: list[str]) -> None:
    print(f"$ {quote_cmd(command)}", flush=True)
    subprocess.run(command, check=True)


def resolve_rust_binary(explicit_path: Path | None, default_path: Path) -> Path:
    if explicit_path is not None:
        return Path(explicit_path)
    return default_path


def existing(path: Path) -> Path:
    path = Path(path)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"not found: {path}")
    return path


def gatk_cmd(java_mem: str, *args: str | Path) -> list[str]:
    return [str(GATK), "--java-options", f"-Xmx{java_mem}", *map(str, args)]


def bam_index_path(bam: Path) -> Path:
    if bam.name.endswith(".bam"):
        return bam.with_suffix(".bai")
    return Path(f"{bam}.bai")


def read_fai_contigs(fai: Path) -> set[str]:
    with fai.open() as handle:
        return {line.split("\t", 1)[0] for line in handle if line.strip()}


def read_fai_contigs_ordered(fai: Path) -> list[str]:
    with fai.open() as handle:
        return [line.split("\t", 1)[0] for line in handle if line.strip()]


def build_exon_bed(gtf: Path, fai: Path, output: Path) -> None:
    contigs = read_fai_contigs(fai)
    opener = gzip.open if gtf.suffix == ".gz" else open
    rows: list[tuple[str, int, int]] = []
    with opener(gtf, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5 or fields[2] != "exon" or fields[0] not in contigs:
                continue
            rows.append((fields[0], int(fields[3]) - 1, int(fields[4])))

    contig_order = {name: idx for idx, name in enumerate(read_fai_contigs_ordered(fai))}
    rows.sort(key=lambda row: (contig_order.get(row[0], sys.maxsize), row[1], row[2]))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for contig, start, end in rows:
            handle.write(f"{contig}\t{start}\t{end}\n")


def build_merged_exon_bed(input_bed: Path, output: Path) -> None:
    current_contig: str | None = None
    current_start: int | None = None
    current_end: int | None = None

    output.parent.mkdir(parents=True, exist_ok=True)
    with input_bed.open() as handle, output.open("w") as out:
        for line in handle:
            if not line.strip():
                continue
            contig, start_str, end_str = line.rstrip("\n").split("\t")
            start = int(start_str)
            end = int(end_str)

            if current_contig is None:
                current_contig = contig
                current_start = start
                current_end = end
                continue

            if contig == current_contig and start <= current_end:
                current_end = max(current_end, end)
                continue

            out.write(f"{current_contig}\t{current_start}\t{current_end}\n")
            current_contig = contig
            current_start = start
            current_end = end

        if current_contig is not None:
            out.write(f"{current_contig}\t{current_start}\t{current_end}\n")


def build_haplotype_caller_command(
    ref: Path,
    input_bam: Path,
    interval_list: Path,
    output_vcf: Path,
    java_mem: str,
    pair_hmm_threads: int,
    dbsnp: Path | None,
) -> list[str]:
    command = gatk_cmd(
        java_mem,
        "HaplotypeCaller",
        "-R",
        ref,
        "-I",
        input_bam,
        "-L",
        interval_list,
        "-O",
        output_vcf,
        "--dont-use-soft-clipped-bases",
        "--standard-min-confidence-threshold-for-calling",
        "20",
        "--native-pair-hmm-threads",
        str(pair_hmm_threads),
    )
    if dbsnp is not None:
        command.extend(["--dbsnp", str(dbsnp)])
    return command

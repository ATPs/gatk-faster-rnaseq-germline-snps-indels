#!/usr/bin/env python3
"""Oracle and comparison utilities for the Rust HaplotypeCaller rewrite."""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from step_common import DEFAULT_DBSNP, DEFAULT_REF, GATK, existing, quote_cmd


LOCAL_GATK = Path("/data/p/gatk/gatk/gatk")
GATK_HC_TEST = (
    Path("/data/p/gatk/gatk")
    / "src/test/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/HaplotypeCallerIntegrationTest.java"
)


@dataclass(frozen=True)
class VcfRecordSet:
    header_ids: dict[str, set[str]]
    info_keys: set[str]
    format_keys: set[str]
    all_records: set[tuple[str, int, str, str]]
    pass_records: set[tuple[str, int, str, str]]
    gvcf_blocks: set[tuple[str, int, int, str]]
    total_rows: int
    pass_rows: int


def main() -> int:
    parser = argparse.ArgumentParser(description="Build and compare HaplotypeCaller oracle fixtures.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    add_inventory_tests_parser(subparsers)
    add_write_interval_parser(subparsers)
    add_run_gatk_parser(subparsers)
    add_compare_vcfs_parser(subparsers)
    add_compare_bam_parser(subparsers)
    add_summarize_vcf_parser(subparsers)

    args = parser.parse_args()
    if args.command == "inventory-tests":
        return inventory_tests(args)
    if args.command == "write-window-interval-list":
        return write_window_interval_list(args)
    if args.command == "run-gatk-hc":
        return run_gatk_hc(args)
    if args.command == "compare-vcfs":
        return compare_vcfs(args)
    if args.command == "compare-bamout":
        return compare_bamout(args)
    if args.command == "summarize-vcf":
        return summarize_vcf(args)
    raise AssertionError(f"unhandled command {args.command}")


def add_inventory_tests_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "inventory-tests",
        help="Inventory HaplotypeCallerIntegrationTest scenarios from the local GATK source.",
    )
    parser.add_argument("--test-java", type=existing, default=GATK_HC_TEST)
    parser.add_argument("--output-tsv", type=Path)
    parser.add_argument("--include-disabled", action="store_true")


def add_write_interval_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "write-window-interval-list",
        help="Write a small interval_list fixture while preserving the input sequence dictionary.",
    )
    parser.add_argument("--input-interval-list", type=existing, required=True)
    parser.add_argument("--output-interval-list", type=Path, required=True)
    parser.add_argument("--region", help="Explicit interval in contig:start-end form.")
    parser.add_argument("--first-records", type=int, default=1)


def add_run_gatk_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "run-gatk-hc",
        help="Run GATK HaplotypeCaller with the pipeline-first Rust-HC acceptance flags.",
    )
    parser.add_argument("--gatk", type=Path, default=LOCAL_GATK)
    parser.add_argument("--fallback-gatk", type=Path, default=GATK)
    parser.add_argument("--use-fallback-gatk", action="store_true")
    parser.add_argument("--ref", type=existing, default=DEFAULT_REF)
    parser.add_argument("--input-bam", type=existing, required=True)
    parser.add_argument("--input-interval-list", type=existing, required=True)
    parser.add_argument("--output-vcf", type=Path, required=True)
    parser.add_argument("--dbsnp", type=existing, default=DEFAULT_DBSNP)
    parser.add_argument("--java-mem", default="24g")
    parser.add_argument("--pair-hmm-threads", type=int, default=8)
    parser.add_argument("--standard-min-confidence-threshold-for-calling", type=float, default=20.0)
    parser.add_argument("--use-soft-clipped-bases", action="store_true")
    parser.add_argument("--dry-run", action="store_true")


def add_compare_vcfs_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "compare-vcfs",
        help="Compare normalized VCF records, headers, annotations, and GVCF block keys.",
    )
    parser.add_argument("--truth-vcf", type=existing, required=True)
    parser.add_argument("--query-vcf", type=existing, required=True)
    parser.add_argument("--ref", type=existing)
    parser.add_argument("--bcftools", default=shutil.which("bcftools") or "bcftools")
    parser.add_argument("--no-normalize", action="store_true")
    parser.add_argument("--output-tsv", type=Path)


def add_compare_bam_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "compare-bamout",
        help="Compare bamout files with samtools flagstat and idxstats summaries.",
    )
    parser.add_argument("--truth-bam", type=existing, required=True)
    parser.add_argument("--query-bam", type=existing, required=True)
    parser.add_argument("--samtools", default=shutil.which("samtools") or "samtools")
    parser.add_argument("--output-tsv", type=Path)


def add_summarize_vcf_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("summarize-vcf", help="Summarize VCF record, PASS, annotation, and GVCF block counts.")
    parser.add_argument("--vcf", type=existing, required=True)
    parser.add_argument("--output-tsv", type=Path)


def inventory_tests(args: argparse.Namespace) -> int:
    rows = []
    pending_test: str | None = None
    with args.test_java.open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("@Test"):
                pending_test = stripped
                continue
            if pending_test is None:
                continue
            match = re.search(r"public\s+void\s+(\w+)\s*\(", stripped)
            if match is None:
                continue
            method = match.group(1)
            enabled = "enabled=false" not in pending_test.replace(" ", "")
            if enabled or args.include_disabled:
                rows.append(
                    {
                        "method": method,
                        "enabled": str(enabled).lower(),
                        "data_provider": annotation_value(pending_test, "dataProvider"),
                        "expected_exceptions": annotation_value(pending_test, "expectedExceptions"),
                    }
                )
            pending_test = None

    write_tsv(args.output_tsv, rows, ["method", "enabled", "data_provider", "expected_exceptions"])
    return 0


def annotation_value(annotation: str, name: str) -> str:
    match = re.search(rf"{name}\s*=\s*([^,\)]+)", annotation)
    if match is None:
        return ""
    return match.group(1).strip().strip('"')


def write_window_interval_list(args: argparse.Namespace) -> int:
    if args.first_records < 1:
        raise SystemExit("--first-records must be at least 1")
    headers: list[str] = []
    body: list[str] = []
    with args.input_interval_list.open() as handle:
        for line in handle:
            if line.startswith("@"):
                headers.append(line)
            elif line.strip() and args.region is None and len(body) < args.first_records:
                body.append(line)

    if args.region is not None:
        contig, start, end = parse_region(args.region)
        body = [f"{contig}\t{start}\t{end}\t+\toracle_window\n"]
    if not body:
        raise SystemExit(f"no interval rows selected from {args.input_interval_list}")

    args.output_interval_list.parent.mkdir(parents=True, exist_ok=True)
    with args.output_interval_list.open("w") as out:
        out.writelines(headers)
        out.writelines(body)
    print(args.output_interval_list)
    return 0


def parse_region(region: str) -> tuple[str, int, int]:
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", region)
    if match is None:
        raise SystemExit("--region must have form contig:start-end")
    contig = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))
    if start < 1 or start > end:
        raise SystemExit("--region coordinates must be 1-based and start <= end")
    return contig, start, end


def run_gatk_hc(args: argparse.Namespace) -> int:
    gatk = args.fallback_gatk if args.use_fallback_gatk else args.gatk
    if not gatk.exists():
        raise SystemExit(f"GATK executable not found: {gatk}")
    if gatk == LOCAL_GATK:
        require_local_gatk_ready(gatk)
    if args.pair_hmm_threads < 1:
        raise SystemExit("--pair-hmm-threads must be at least 1")

    command = [
        str(gatk),
        "--java-options",
        f"-Xmx{args.java_mem}",
        "HaplotypeCaller",
        "-R",
        str(args.ref),
        "-I",
        str(args.input_bam),
        "-L",
        str(args.input_interval_list),
        "-O",
        str(args.output_vcf),
        "--standard-min-confidence-threshold-for-calling",
        format_threshold(args.standard_min_confidence_threshold_for_calling),
        "--native-pair-hmm-threads",
        str(args.pair_hmm_threads),
    ]
    if not args.use_soft_clipped_bases:
        command.append("--dont-use-soft-clipped-bases")
    if args.dbsnp is not None:
        command.extend(["--dbsnp", str(args.dbsnp)])

    print(f"$ {quote_cmd(command)}", flush=True)
    if args.dry_run:
        return 0
    args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(command, check=True)
    print(args.output_vcf)
    return 0


def require_local_gatk_ready(gatk: Path) -> None:
    source_dir = gatk.parent
    if os.environ.get("GATK_LOCAL_JAR"):
        return
    if any((source_dir / "build/libs").glob("*local.jar")):
        return
    raise SystemExit(
        f"{gatk} is the local GATK source wrapper, but no local jar is built. "
        f"Run: {source_dir / 'gradlew'} localJar, set GATK_LOCAL_JAR, or pass --use-fallback-gatk."
    )


def format_threshold(value: float) -> str:
    return str(int(value)) if value.is_integer() else str(value)


def compare_vcfs(args: argparse.Namespace) -> int:
    normalize = not args.no_normalize
    if normalize and args.ref is None:
        raise SystemExit("--ref is required unless --no-normalize is set")
    truth = load_vcf_records(args.truth_vcf, args.ref, args.bcftools, normalize)
    query = load_vcf_records(args.query_vcf, args.ref, args.bcftools, normalize)

    all_shared = truth.all_records & query.all_records
    pass_shared = truth.pass_records & query.pass_records
    rows = [
        metric("truth_total_rows", truth.total_rows),
        metric("query_total_rows", query.total_rows),
        metric("truth_pass_rows", truth.pass_rows),
        metric("query_pass_rows", query.pass_rows),
        metric("shared_records", len(all_shared)),
        metric("truth_private_records", len(truth.all_records - query.all_records)),
        metric("query_private_records", len(query.all_records - truth.all_records)),
        metric("shared_pass_records", len(pass_shared)),
        metric("truth_private_pass_records", len(truth.pass_records - query.pass_records)),
        metric("query_private_pass_records", len(query.pass_records - truth.pass_records)),
        metric("truth_gvcf_blocks", len(truth.gvcf_blocks)),
        metric("query_gvcf_blocks", len(query.gvcf_blocks)),
        metric("shared_gvcf_blocks", len(truth.gvcf_blocks & query.gvcf_blocks)),
    ]
    for section in ("INFO", "FORMAT", "FILTER", "contig"):
        rows.append(metric(f"truth_header_{section}_ids", len(truth.header_ids.get(section, set()))))
        rows.append(metric(f"query_header_{section}_ids", len(query.header_ids.get(section, set()))))
        rows.append(
            metric(
                f"query_missing_header_{section}_ids",
                ",".join(sorted(truth.header_ids.get(section, set()) - query.header_ids.get(section, set()))),
            )
        )
        rows.append(
            metric(
                f"query_extra_header_{section}_ids",
                ",".join(sorted(query.header_ids.get(section, set()) - truth.header_ids.get(section, set()))),
            )
        )
    rows.extend(
        [
            metric("query_missing_info_keys", ",".join(sorted(truth.info_keys - query.info_keys))),
            metric("query_extra_info_keys", ",".join(sorted(query.info_keys - truth.info_keys))),
            metric("query_missing_format_keys", ",".join(sorted(truth.format_keys - query.format_keys))),
            metric("query_extra_format_keys", ",".join(sorted(query.format_keys - truth.format_keys))),
        ]
    )
    write_tsv(args.output_tsv, rows, ["metric", "value"])
    return 0


def summarize_vcf(args: argparse.Namespace) -> int:
    records = load_vcf_records(args.vcf, None, "bcftools", normalize=False)
    rows = [
        metric("total_rows", records.total_rows),
        metric("pass_rows", records.pass_rows),
        metric("nonpass_rows", records.total_rows - records.pass_rows),
        metric("allele_records", len(records.all_records)),
        metric("pass_allele_records", len(records.pass_records)),
        metric("gvcf_blocks", len(records.gvcf_blocks)),
        metric("info_keys", ",".join(sorted(records.info_keys))),
        metric("format_keys", ",".join(sorted(records.format_keys))),
    ]
    write_tsv(args.output_tsv, rows, ["metric", "value"])
    return 0


def load_vcf_records(path: Path, ref: Path | None, bcftools: str, normalize: bool) -> VcfRecordSet:
    if normalize:
        command = [bcftools, "norm", "-f", str(ref), "-m", "-any", str(path)]
        completed = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        lines = completed.stdout.splitlines()
    else:
        with open_text(path) as handle:
            lines = [line.rstrip("\n") for line in handle]
    return parse_vcf_lines(lines)


def parse_vcf_lines(lines: list[str]) -> VcfRecordSet:
    header_ids: dict[str, set[str]] = {"INFO": set(), "FORMAT": set(), "FILTER": set(), "contig": set()}
    info_keys: set[str] = set()
    format_keys: set[str] = set()
    all_records: set[tuple[str, int, str, str]] = set()
    pass_records: set[tuple[str, int, str, str]] = set()
    gvcf_blocks: set[tuple[str, int, int, str]] = set()
    total_rows = 0
    pass_rows = 0

    for line in lines:
        if not line:
            continue
        if line.startswith("##"):
            capture_header_id(line, header_ids)
            continue
        if line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) < 8:
            raise SystemExit(f"malformed VCF row has fewer than 8 columns: {line}")
        chrom, pos_s, _id, ref, alts, _qual, filter_value, info = fields[:8]
        format_value = fields[8] if len(fields) > 8 else ""
        pos = int(pos_s)
        total_rows += 1
        is_pass = filter_value in (".", "PASS")
        if is_pass:
            pass_rows += 1
        for alt in alts.split(","):
            key = (chrom, pos, ref, alt)
            all_records.add(key)
            if is_pass:
                pass_records.add(key)
        info_keys.update(parse_info_keys(info))
        if format_value:
            format_keys.update(format_value.split(":"))
        end = info_end(info)
        if end is not None or any(alt.startswith("<") for alt in alts.split(",")):
            gvcf_blocks.add((chrom, pos, end or pos, alts))

    return VcfRecordSet(
        header_ids=header_ids,
        info_keys=info_keys,
        format_keys=format_keys,
        all_records=all_records,
        pass_records=pass_records,
        gvcf_blocks=gvcf_blocks,
        total_rows=total_rows,
        pass_rows=pass_rows,
    )


def capture_header_id(line: str, header_ids: dict[str, set[str]]) -> None:
    match = re.match(r"##(INFO|FORMAT|FILTER|contig)=<ID=([^,>]+)", line)
    if match is not None:
        header_ids[match.group(1)].add(match.group(2))


def parse_info_keys(info: str) -> set[str]:
    if info in ("", "."):
        return set()
    return {entry.split("=", 1)[0] for entry in info.split(";") if entry}


def info_end(info: str) -> int | None:
    for entry in info.split(";"):
        if entry.startswith("END="):
            return int(entry.split("=", 1)[1])
    return None


def compare_bamout(args: argparse.Namespace) -> int:
    truth_flagstat = samtools_summary(args.samtools, "flagstat", args.truth_bam)
    query_flagstat = samtools_summary(args.samtools, "flagstat", args.query_bam)
    truth_idxstats = samtools_summary(args.samtools, "idxstats", args.truth_bam)
    query_idxstats = samtools_summary(args.samtools, "idxstats", args.query_bam)
    rows = [
        metric("flagstat_identical", str(truth_flagstat == query_flagstat).lower()),
        metric("idxstats_identical", str(truth_idxstats == query_idxstats).lower()),
        metric("truth_flagstat_lines", len(truth_flagstat)),
        metric("query_flagstat_lines", len(query_flagstat)),
        metric("truth_idxstats_lines", len(truth_idxstats)),
        metric("query_idxstats_lines", len(query_idxstats)),
    ]
    write_tsv(args.output_tsv, rows, ["metric", "value"])
    return 0


def samtools_summary(samtools: str, subcommand: str, bam: Path) -> list[str]:
    completed = subprocess.run([samtools, subcommand, str(bam)], check=True, stdout=subprocess.PIPE, text=True)
    return completed.stdout.splitlines()


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def metric(name: str, value: object) -> dict[str, str]:
    return {"metric": name, "value": str(value)}


def write_tsv(path: Path | None, rows: list[dict[str, str]], columns: list[str]) -> None:
    lines = ["\t".join(columns)]
    lines.extend("\t".join(row.get(column, "") for column in columns) for row in rows)
    output = "\n".join(lines) + "\n"
    if path is None:
        print(output, end="")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(output)
    print(path)


if __name__ == "__main__":
    raise SystemExit(main())

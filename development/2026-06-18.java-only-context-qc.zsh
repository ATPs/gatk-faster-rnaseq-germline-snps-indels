#!/usr/bin/env zsh

set -euo pipefail

source /data/p/anaconda3/etc/profile.d/conda.sh
conda activate base
export PATH=/data/p/bin:$PATH

REF=/data1/pub/gatk/broad_hg38/Homo_sapiens_assembly38.fasta
ROUND_DIR=/XCLabServer002_fastIO/gatk-faster-rnaseq/SRR949115_broad_hg38_run2/module6_good_compare_20260618_round26
OUT=${ROUND_DIR}/round26.java_only_selected_context_qc.tsv

export REF ROUND_DIR OUT

python - <<'PY'
import csv
import os
import subprocess
import sys

REF = os.environ["REF"]
ROUND_DIR = os.environ["ROUND_DIR"]
OUT = os.environ["OUT"]
WINDOW = 20
CLOSE_BP = 5

VCFS = {
    "java_filtered": os.path.join(ROUND_DIR, "java.split.vcf.gz"),
    "rust_filtered": os.path.join(ROUND_DIR, "rust.split.vcf.gz"),
    "java_raw": os.path.join(ROUND_DIR, "java.raw.split.vcf.gz"),
    "rust_raw": os.path.join(ROUND_DIR, "rust.raw.split.vcf.gz"),
}

VARIANTS = """\
chr1	1663538	T	C
chr1	1732521	G	A
chr1	37861061	T	C
chr1	85653613	C	A
chr1	109497688	A	G
chr1	153934830	T	C
chr10	78040667	T	A
chr10	100232510	T	C
chr11	73393100	A	C
chr11	119001795	A	G
chr12	53297151	A	G
chr12	53297349	A	G
chr13	40941150	T	A
chr13	49701213	T	A
chr14	101926719	G	T
chr15	40694122	T	C
chr15	67066140	G	A
chr16	18926186	G	A
chr16	29696360	A	G
chr16	88650978	T	C
chr17	27313055	A	G
chr17	40019203	C	T
chr17	42313402	G	A
chr17	57681285	A	G
chr17	82155505	G	C
chr18	32070755	A	G
chr19	2037650	C	T
chr19	18261356	G	A
chr2	88947494	C	A
chr2	232583365	C	T
chr2	234495152	T	A
chr2	234495153	A	T
chr20	2836116	A	G
chr20	45367549	G	A
chr20	62545795	C	G
chr22	22322790	G	A
chr22	22713133	G	A
chr3	113652120	T	A
chr3	113652122	T	A
chr4	24530682	G	A
chr5	160402051	A	T
chr6	18166300	G	A
chr6	137875533	G	A
chr7	7607020	G	T
chr7	87187270	A	G
chr7	127392399	G	A
chr8	673593	G	A
chr8	143925474	A	C
chr8	143927420	C	T
chrX	301507	A	G
chrX	49071643	G	A
"""


def run_text(args):
    return subprocess.check_output(args, text=True)


def load_lengths(fai_path):
    lengths = {}
    with open(fai_path, "r", encoding="ascii") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def fetch_context(chrom, pos, lengths):
    start = max(1, pos - WINDOW)
    end = min(lengths[chrom], pos + WINDOW)
    region = f"{chrom}:{start}-{end}"
    text = run_text(["samtools", "faidx", REF, region])
    seq = "".join(line.strip() for line in text.splitlines() if not line.startswith(">")).upper()
    return start, end, seq


def fmt_offset(offset):
    if offset > 0:
        return f"+{offset}"
    return str(offset)


def n_summary(seq, start, pos):
    offsets = []
    for idx, base in enumerate(seq):
        if base == "N":
            offsets.append(start + idx - pos)
    if not offsets:
        return "0", "0", "0", ".", "."
    min_abs = min(abs(offset) for offset in offsets)
    return (
        "1",
        str(len(offsets)),
        str(min_abs),
        ",".join(fmt_offset(offset) for offset in offsets),
        ",".join(str(pos + offset) for offset in offsets),
    )


def add_repeat(repeats, motif, start_idx, end_idx, kind, repeat_count=None):
    repeats.append(
        {
            "motif": motif,
            "start_idx": start_idx,
            "end_idx": end_idx,
            "kind": kind,
            "repeat_count": repeat_count,
        }
    )


def find_simple_repeats(seq):
    repeats = []
    n = len(seq)

    i = 0
    while i < n:
        j = i + 1
        while j < n and seq[j] == seq[i]:
            j += 1
        if seq[i] != "N" and j - i >= 5:
            add_repeat(repeats, seq[i], i, j - 1, "homopolymer", j - i)
        i = j

    for motif_len in range(2, 7):
        i = 0
        while i <= n - motif_len * 3:
            motif = seq[i : i + motif_len]
            if "N" in motif or len(set(motif)) == 1:
                i += 1
                continue
            count = 1
            j = i + motif_len
            while j + motif_len <= n and seq[j : j + motif_len] == motif:
                count += 1
                j += motif_len
            if count >= 3:
                add_repeat(repeats, motif, i, j - 1, "tandem", count)
                i = j
            else:
                i += 1

    unique = {}
    for repeat in repeats:
        key = (repeat["start_idx"], repeat["end_idx"], repeat["motif"])
        prev = unique.get(key)
        if prev is None or len(repeat["motif"]) < len(prev["motif"]):
            unique[key] = repeat
    return sorted(unique.values(), key=lambda item: (item["start_idx"], item["end_idx"], item["motif"]))


def repeat_summary(seq, start, pos):
    repeats = find_simple_repeats(seq)
    if not repeats:
        return "0", "0", "."

    pos_idx = pos - start
    variant_in_repeat = any(r["start_idx"] <= pos_idx <= r["end_idx"] for r in repeats)
    notes = []
    for repeat in repeats[:8]:
        rel_start = start + repeat["start_idx"] - pos
        rel_end = start + repeat["end_idx"] - pos
        label = repeat["kind"]
        motif = repeat["motif"]
        count = repeat["repeat_count"]
        notes.append(f"{label}:{motif}@{fmt_offset(rel_start)}..{fmt_offset(rel_end)}x{count}")
    if len(repeats) > 8:
        notes.append(f"...{len(repeats) - 8}_more")
    return "1", "1" if variant_in_repeat else "0", ";".join(notes)


def query_close_calls(vcf, chrom, pos, ref, alt):
    start = max(1, pos - CLOSE_BP + 1)
    end = pos + CLOSE_BP - 1
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%QUAL[\\t%GT\\t%AD\\t%DP]\\n"
    text = run_text(["bcftools", "query", "-r", f"{chrom}:{start}-{end}", "-f", fmt, vcf])
    calls = []
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        rec_chrom = fields[0]
        rec_pos = int(fields[1])
        rec_ref = fields[2]
        rec_alt = fields[3]
        if rec_chrom == chrom and rec_pos == pos and rec_ref == ref and rec_alt == alt:
            continue
        delta = rec_pos - pos
        if abs(delta) >= CLOSE_BP:
            continue
        gt = fields[6] if len(fields) > 6 else "."
        ad = fields[7] if len(fields) > 7 else "."
        dp = fields[8] if len(fields) > 8 else "."
        calls.append(f"{fmt_offset(delta)}:{rec_pos}:{rec_ref}>{rec_alt}:{fields[4]}:{gt}:{ad}:{dp}")
    return str(len(calls)), ";".join(calls) if calls else "."


def parse_variants():
    rows = []
    for line in VARIANTS.strip().splitlines():
        chrom, pos, ref, alt = line.split("\t")
        rows.append((chrom, int(pos), ref, alt))
    return rows


def main():
    fai = REF + ".fai"
    if not os.path.exists(fai):
        raise FileNotFoundError(f"missing FASTA index: {fai}")
    for label, path in VCFS.items():
        if not os.path.exists(path):
            raise FileNotFoundError(f"missing {label} VCF: {path}")

    lengths = load_lengths(fai)
    header = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "context_window_bp",
        "context_start",
        "context_end",
        "ref_context",
        "has_N_within_20bp",
        "N_count_within_20bp",
        "N_min_abs_offset",
        "N_offsets",
        "N_positions",
        "simple_repeat_within_20bp",
        "variant_base_in_simple_repeat",
        "simple_repeat_notes",
    ]
    for label in ("java_filtered", "rust_filtered", "java_raw", "rust_raw"):
        header.extend([f"{label}_close_lt5bp_count", f"{label}_close_lt5bp_calls"])

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for chrom, pos, ref, alt in parse_variants():
            if chrom not in lengths:
                raise KeyError(f"{chrom} not present in FASTA index")
            context_start, context_end, seq = fetch_context(chrom, pos, lengths)
            if seq[pos - context_start] != ref:
                print(
                    f"warning: REF mismatch at {chrom}:{pos}: input {ref}, fasta {seq[pos - context_start]}",
                    file=sys.stderr,
                )
            n_cols = n_summary(seq, context_start, pos)
            repeat_cols = repeat_summary(seq, context_start, pos)
            row = [
                chrom,
                str(pos),
                ref,
                alt,
                str(WINDOW),
                str(context_start),
                str(context_end),
                seq,
                *n_cols,
                *repeat_cols,
            ]
            for label in ("java_filtered", "rust_filtered", "java_raw", "rust_raw"):
                row.extend(query_close_calls(VCFS[label], chrom, pos, ref, alt))
            writer.writerow(row)


if __name__ == "__main__":
    main()
PY

echo "${OUT}"

# gatk-faster-rnaseq-germline-snps-indels

An RNA-seq germline SNP/indel pipeline that keeps the Broad/GATK RNA calling structure but replaces selected slow steps with faster implementations, mainly in Rust.

Input can be paired FASTQ files or an already aligned STAR coordinate-sorted BAM. Output is a raw VCF plus the final filtered VCF.

## Status

The repository already supports a practical end-to-end pipeline with:

- STAR alignment from FASTQ
- GATK-compatible RNA preprocessing
- backend selection per step
- `sambamba` or Rust duplicate marking
- Rust interval preparation, SplitNCigarReads, BaseRecalibrator, ApplyBQSR, and HaplotypeCaller candidate prefiltering
- scattered HaplotypeCaller execution
- step skipping and artifact reuse for resume/debug workflows

The full Rust HaplotypeCaller is present but still experimental. GATK remains the default caller backend.

## Repository Layout

- `gatk-faster-rnaseq-germline-snps-indels.py`: top-level pipeline entrypoint
- `src/gatk_faster_rnaseq/`: Python package for orchestration and step wrappers
- `src/rust/`: Rust crate for accelerated tools
- `rust_binary/`: installed release Rust executables used by the pipeline
- `build.md`: build guide
- `design.md`: architecture and tool-by-tool design notes
- `development/`: background notes, benchmarks, and experiment logs

## Requirements

Core external tools used by this project:

- GATK 4
- STAR
- samtools
- sambamba
- Python 3
- Rust and Cargo

See `build.md` for the portable build flow and Rust binary installation steps.

## Important Note About Defaults

The current Python step defaults point to paths on the original development environment, including:

- the SRR949115 test dataset
- Broad hg38 resources
- a prebuilt STAR index
- the default output directory under `/XCLabServer002_fastIO`

On another machine, pass your own paths explicitly.

## Quick Start

### 1. Build and install the Rust binaries

Follow `build.md`.

### 2. Run the full pipeline from FASTQ

Example:

```bash
python gatk-faster-rnaseq-germline-snps-indels.py \
  --sample SAMPLE1 \
  --fastq1 /path/to/SAMPLE1_1.fastq.gz \
  --fastq2 /path/to/SAMPLE1_2.fastq.gz \
  --ref /path/to/Homo_sapiens_assembly38.fasta \
  --ref-dict /path/to/Homo_sapiens_assembly38.dict \
  --gtf /path/to/gencode.annotation.gtf \
  --star-index /path/to/star.index \
  --dbsnp /path/to/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  --known-sites /path/to/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  --known-sites /path/to/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites /path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --outdir /path/to/output \
  --threads 40
```

### 3. Start from an already aligned STAR BAM

Example:

```bash
python gatk-faster-rnaseq-germline-snps-indels.py \
  --sample SAMPLE1 \
  --aligned-bam /path/to/SAMPLE1.Aligned.sortedByCoord.out.bam \
  --skip-star-align \
  --ref /path/to/Homo_sapiens_assembly38.fasta \
  --ref-dict /path/to/Homo_sapiens_assembly38.dict \
  --gtf /path/to/gencode.annotation.gtf \
  --dbsnp /path/to/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  --known-sites /path/to/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
  --known-sites /path/to/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites /path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --outdir /path/to/output \
  --threads 40
```

## Backend Control

Main backend switches:

- `--mode baseline|sambamba|rust|auto` controls duplicate marking
- `--interval-backend gatk|rust|auto`
- `--split-n-cigar-backend gatk|rust|auto`
- `--base-recalibrator-backend gatk|rust|auto`
- `--apply-bqsr-backend gatk|rust|auto`
- `--hc-prefilter-backend none|rust|auto`
- `--hc-backend gatk|rust`
- `--no-rust` forces the non-Rust path

The usual safe default is `auto`, which uses installed Rust binaries when present and otherwise falls back to GATK.

## Resume And Debug Workflows

The runner can skip completed or intentionally omitted steps and continue from supplied artifacts:

- `--skip-star-align` with `--aligned-bam`
- `--skip-mark-duplicates` with `--dedup-bam`
- `--skip-split-n-cigar` with `--split-bam`
- `--skip-base-recalibrator` with `--recal-table`
- `--skip-apply-bqsr` with `--recal-bam`
- `--skip-haplotype-caller` with `--raw-vcf`
- `--skip-bqsr` as a shortcut for skipping both BQSR steps

This is useful for benchmarking one stage at a time and for reusing large intermediate files.

## Outputs

A run creates a labeled subdirectory under `--outdir` containing:

- final BAM and VCF artifacts
- `logs/` with one log per step
- `timings.tsv` with step timing and status
- `reference/` with generated BED and `interval_list` files
- HaplotypeCaller shard outputs when scatter mode is enabled

The final output is normally `SAMPLE.filtered.vcf.gz`. If `--skip-variant-filtration` is used, the raw VCF is left as the final output.

## Running Individual Steps

Each major step is also callable directly, for example:

```bash
export PYTHONPATH="$PWD/src"
python -m gatk_faster_rnaseq.steps.haplotype_caller --help
python -m gatk_faster_rnaseq.steps.apply_bqsr --help
python -m gatk_faster_rnaseq.steps.split_n_cigar --help
```

## Documentation

- `design.md`: current architecture and every tool in the repository
- `build.md`: build guide
- `development/background.md`: background context
- `development/*.md`: experiment records and validation notes

## License

MIT License. See `LICENSE`.

# Project Design

## Purpose

`gatk-faster-rnaseq-germline-snps-indels` is a single-sample RNA-seq germline SNP and indel pipeline. It accepts either:

- paired FASTQ files, or
- an already aligned STAR coordinate-sorted BAM

and produces a filtered VCF.

The design goal is not to invent a new calling workflow from scratch. The project keeps the Broad/GATK RNA-seq germline calling structure, then replaces slow stages with faster implementations when that can be done with acceptable compatibility.

## Current State

The repository already supports:

- a baseline GATK-oriented path
- faster duplicate marking with `sambamba`
- Rust replacements for several preprocessing steps
- optional interval prefiltering before calling
- scattered HaplotypeCaller execution
- artifact reuse and step skipping for resume/debug workflows

The full Rust HaplotypeCaller exists, but it is still opt-in and under active validation. The default caller backend remains GATK.

## High-Level Architecture

### 1. Thin top-level entrypoint

- `gatk-faster-rnaseq-germline-snps-indels.py`

This script only adds `src/` to `sys.path` and dispatches into the Python package.

### 2. Python package as the control plane

- `src/gatk_faster_rnaseq/cli.py`
- `src/gatk_faster_rnaseq/pipeline/runner.py`
- `src/gatk_faster_rnaseq/runtime.py`

This layer owns orchestration, argument parsing, output layout, logging, timing, skip logic, and backend selection.

### 3. Independently runnable step wrappers

- `src/gatk_faster_rnaseq/steps/*.py`

Each step can be run on its own. These wrappers normalize arguments and choose between GATK, `sambamba`, or Rust backends where available.

### 4. Rust crate as the acceleration layer

- `src/rust/src/`
- `src/rust/cli/`

The Rust code is split into library modules plus thin executable entrypoints. Release binaries are installed into `rust_binary/` and used by the Python pipeline.

## Pipeline Data Flow

The main pipeline in `src/gatk_faster_rnaseq/pipeline/runner.py` works as follows:

1. Prepare reference-side interval assets from the GTF.
2. Run STAR alignment when starting from FASTQ.
3. Mark duplicates with GATK, `sambamba`, or Rust.
4. Run SplitNCigarReads with GATK or Rust.
5. Run BaseRecalibrator with GATK or Rust.
6. Run ApplyBQSR with GATK or Rust.
7. Optionally shrink calling intervals with the Rust HaplotypeCaller prefilter.
8. Run HaplotypeCaller on one interval set or on scattered shards.
9. Merge shard VCFs if scattering was used.
10. Apply standard RNA hard filters with GATK VariantFiltration.

## Execution Model

### Step abstraction

`src/gatk_faster_rnaseq/runtime.py` defines a `Step` dataclass with:

- step name
- command
- expected outputs
- optional environment overrides

The runtime runner executes steps, writes one log per step, and appends timing records to `timings.tsv`.

### Resumability

If all declared outputs for a step already exist, the step is skipped unless `--force` is used.

### Artifact reuse

The runner can resume from intermediate files using:

- `--aligned-bam`
- `--dedup-bam`
- `--split-bam`
- `--recal-table`
- `--recal-bam`
- `--raw-vcf`

### Parallelism

The heavy parallel path today is HaplotypeCaller scattering:

- `split_hc_intervals` creates shards
- shard callers run in parallel
- `merge_vcfs` combines results

## Output Layout

For each run label under `--outdir`, the pipeline writes:

- final BAM and VCF artifacts
- `logs/` with one log file per step
- `timings.tsv` with wall-clock timings and skip/fail status
- `reference/` with generated BED and `interval_list` files
- `hc_shards/` and `reference/hc_scatter/` when scattered calling is enabled

## Tool Catalog

This section explains every tool currently present in the project.

### External tools used by the pipeline

- `GATK`: baseline implementation for dictionary creation, interval conversion, MarkDuplicates, SplitNCigarReads, BaseRecalibrator, ApplyBQSR, HaplotypeCaller, MergeVcfs, and VariantFiltration.
- `STAR`: two-pass RNA-seq aligner used when the input is FASTQ.
- `samtools`: FASTA indexing and BAM indexing support.
- `sambamba`: multi-threaded duplicate marking alternative to Picard/GATK.
- `Python`: pipeline orchestration and step wrappers.
- `Rust/Cargo`: build system and runtime for accelerated internal tools.

### Python user-facing tools

- `gatk-faster-rnaseq-germline-snps-indels.py`: top-level end-to-end pipeline launcher.
- `src/gatk_faster_rnaseq/steps/star_align.py`: STAR two-pass alignment wrapper.
- `src/gatk_faster_rnaseq/steps/build_exon_bed.py`: extract exon intervals from a GTF, filtered and sorted by reference contigs.
- `src/gatk_faster_rnaseq/steps/build_merged_exon_bed.py`: merge overlapping exon BED rows into non-overlapping regions.
- `src/gatk_faster_rnaseq/steps/bed_to_interval_list.py`: convert BED intervals into Picard/GATK `interval_list` format using GATK or Rust.
- `src/gatk_faster_rnaseq/steps/mark_duplicates.py`: choose duplicate marking backend: baseline GATK, `sambamba`, or Rust.
- `src/gatk_faster_rnaseq/steps/split_n_cigar.py`: choose GATK or Rust SplitNCigarReads.
- `src/gatk_faster_rnaseq/steps/base_recalibrator.py`: choose GATK or Rust BaseRecalibrator.
- `src/gatk_faster_rnaseq/steps/apply_bqsr.py`: choose GATK or Rust ApplyBQSR.
- `src/gatk_faster_rnaseq/steps/hc_prefilter.py`: call the Rust prefilter that keeps only intervals with pileup evidence for variation.
- `src/gatk_faster_rnaseq/steps/split_hc_intervals.py`: scatter an `interval_list` into balanced shards with GATK or Rust.
- `src/gatk_faster_rnaseq/steps/haplotype_caller.py`: run GATK HaplotypeCaller or the opt-in Rust HaplotypeCaller on one interval set.
- `src/gatk_faster_rnaseq/steps/merge_vcfs.py`: merge scattered VCF shards with GATK MergeVcfs.
- `src/gatk_faster_rnaseq/steps/variant_filtration.py`: apply the standard RNA hard filters (`FS > 30.0`, `QD < 2.0`).

### Python internal support modules

- `src/gatk_faster_rnaseq/runtime.py`: step execution, log writing, timing, and limited batch parallel execution.
- `src/gatk_faster_rnaseq/steps/common.py`: shared defaults, tool paths, command builders, BED helpers, and path validation.
- `src/gatk_faster_rnaseq/pipeline/runner.py`: main orchestration logic, backend normalization, skip validation, run labeling, and shard wiring.

### Rust pipeline binaries

These are the Rust tools the Python pipeline can call directly.

- `rust_interval_tools`: interval utility binary with three subcommands:
  - `bed-to-interval-list`: convert BED to sorted, merged `interval_list`
  - `split-intervals`: shard an `interval_list` into balanced pieces
  - `prepare`: do BED merge, `interval_list` creation, and sharding in one tool
- `rust_split_n_cigar_reads`: Rust replacement for SplitNCigarReads, with `fast` and `compatibility` modes.
- `rust_base_recalibrator`: mismatch-only Rust BaseRecalibrator replacement for RNA-seq preprocessing.
- `rust_apply_bqsr`: apply recalibration tables to a BAM in Rust.
- `rust_mark_duplicates`: Rust duplicate marking for coordinate-sorted BAMs.
- `rust_hc_prefilter`: scan pileups and emit candidate calling intervals before HaplotypeCaller.
- `rust_haplotype_caller`: assembly-based Rust HaplotypeCaller workbench and pipeline-callable variant caller.

### Rust HaplotypeCaller validation and debug tools

These tools are not core preprocessing steps. They exist to compare Rust behavior with Java GATK and accelerate development.

- `rust_hc_vcf_compare`: exact-allele-key VCF comparison between two callers.
- `rust_hc_select_regions`: select focused debug regions from two VCFs.
- `rust_hc_run_gatk_oracle`: run Java GATK HaplotypeCaller with debug/oracle outputs enabled.
- `rust_hc_gatk_debug_extract`: parse GATK debug output into structured TSV tables.
- `rust_hc_region_replay`: replay chosen regions through the Rust caller and emit stage tables.
- `rust_hc_vcf_to_genotype_table`: convert a VCF into genotype-stage TSV rows.
- `rust_hc_active_region_diff`: diff active-region stage tables between Java and Rust.
- `rust_hc_read_finalize_diff`: diff finalized-read stage tables between Java and Rust.
- `rust_hc_assembly_diff`: diff assembly and haplotype stage tables.
- `rust_hc_pairhmm_diff`: diff PairHMM likelihood tables.
- `rust_hc_genotype_diff`: diff genotype and annotation tables.
- `rust_hc_acceptance_report`: combine comparison outputs into one acceptance report.

### Rust library modules

These modules implement the real acceleration logic behind the binaries.

- `interval_tools`: interval parsing, sorting, merging, and scattering.
- `split_n_cigar`: record transformation logic for RNA splice-aware BAM cleanup.
- `base_recalibrator`: recalibration-table construction.
- `apply_bqsr`: quality-score adjustment using recal tables.
- `mark_duplicates`: duplicate detection and metrics generation.
- `hc_prefilter`: pileup-driven candidate interval selection.
- `haplotype_caller`: Rust caller orchestration and stage-level logic.
- `assembly`: local haplotype assembly helpers.
- `smith_waterman`: haplotype-to-reference alignment utilities.
- `pair_hmm`: read-to-haplotype likelihood scoring.
- `hc_tools`: comparison and reporting helpers for caller validation.
- `gatk_oracle`: wrapper for running Java GATK with debug outputs.

## Why The Design Looks Like This

The current design follows a few deliberate rules:

- Replace slow stages incrementally instead of attempting a one-shot rewrite.
- Keep the Python layer thin and operational, not algorithm-heavy.
- Keep Rust binaries thin and move real logic into library modules.
- Make every major step runnable in isolation for benchmarking and debugging.
- Prefer explicit artifacts and skip flags over hidden state.
- Preserve a clear fallback path to GATK for correctness checks.

## Main Limitations

- The repository is still tuned around one sample at a time.
- Several default resource paths are specific to the original development server.
- The Rust HaplotypeCaller is not yet a drop-in replacement for GATK in all cases.
- Validation tooling around HaplotypeCaller is extensive because exact parity is still an active engineering target.

## Related Documents

- `README.md`: user-facing overview and quick start
- `build.md`: build and installation details
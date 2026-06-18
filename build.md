# Build Guide

This file is the portable build guide for contributors who are not on the original development server.

For the archived server-specific commands and historical benchmark notes, see [development/2026-06-18.server-specific-build.md](/data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/development/2026-06-18.server-specific-build.md).

## Repository layout

- Python entrypoint: `gatk-faster-rnaseq-germline-snps-indels.py`
- Python package: `src/gatk_faster_rnaseq`
- Rust crate root: `src/rust`
- Rust library modules: `src/rust/src`
- Rust CLI entrypoints: `src/rust/cli`
- Installed Rust binaries used by the Python pipeline: `rust_binary`

## Prerequisites

You need:

- Python 3
- Rust and Cargo
- a C/C++ build toolchain
- `libclang` and Clang headers for `bindgen`
- system libraries required by `rust-htslib` / `hts-sys`

Typical Linux package names vary by distribution, but you usually need packages equivalent to:

- `clang`
- `libclang-dev`
- `pkg-config`
- `zlib` development headers
- `bz2` development headers
- `xz` / `lzma` development headers
- `libcurl` development headers

## Important local override note

This repository currently includes [src/rust/.cargo/config.toml](/data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/src/rust/.cargo/config.toml:1), and that file may contain environment variables with paths from the original development machine.

Before building on a different machine, review that file:

- if the paths are valid for your environment, keep them
- if they are not valid, update or remove them before building

If `bindgen` cannot find `libclang`, set `LIBCLANG_PATH` for your local installation before running Cargo.

## Build Rust binaries

Use a normal release build first. It is the most portable default.

```bash
cd /path/to/gatk-faster-rnaseq-germline-snps-indels/src/rust
cargo build --release --bins
```

Build one binary only:

```bash
cd /path/to/gatk-faster-rnaseq-germline-snps-indels/src/rust
cargo build --release --bin rust_haplotype_caller
```

## Install built binaries for the pipeline

The Python pipeline looks for installed binaries in `rust_binary`.

```bash
cd /path/to/gatk-faster-rnaseq-germline-snps-indels
mkdir -p rust_binary
for bin in \
  rust_interval_tools \
  rust_split_n_cigar_reads \
  rust_base_recalibrator \
  rust_apply_bqsr \
  rust_mark_duplicates \
  rust_hc_prefilter \
  rust_haplotype_caller \
  rust_hc_vcf_compare \
  rust_hc_select_regions \
  rust_hc_run_gatk_oracle \
  rust_hc_gatk_debug_extract \
  rust_hc_region_replay \
  rust_hc_vcf_to_genotype_table \
  rust_hc_active_region_diff \
  rust_hc_read_finalize_diff \
  rust_hc_assembly_diff \
  rust_hc_pairhmm_diff \
  rust_hc_genotype_diff \
  rust_hc_acceptance_report
do
  install -m 755 "src/rust/target/release/$bin" rust_binary/
done
```

## Verify the build

Run formatting and tests:

```bash
cd /path/to/gatk-faster-rnaseq-germline-snps-indels/src/rust
cargo fmt --check
cargo test --release
```

Check one installed binary:

```bash
cd /path/to/gatk-faster-rnaseq-germline-snps-indels
rust_binary/rust_haplotype_caller --help
```

## Verify the Python entrypoints

The repository entry script works directly:

```bash
cd /path/to/gatk-faster-rnaseq-germline-snps-indels
python gatk-faster-rnaseq-germline-snps-indels.py --help
```

The package step modules need `src` on `PYTHONPATH`:

```bash
cd /path/to/gatk-faster-rnaseq-germline-snps-indels
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python -m gatk_faster_rnaseq.steps.haplotype_caller --help
python -m gatk_faster_rnaseq.steps.apply_bqsr --help
```

## Platform-specific builds

Static builds, custom linker flags, and server-specific environment exports depend on the target machine and libc/toolchain setup. Those commands are intentionally not the default guide here.

If you need the original server-specific static build flow, use [development/2026-06-18.server-specific-build.md](/data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/development/2026-06-18.server-specific-build.md) as a reference and adapt it to your environment.

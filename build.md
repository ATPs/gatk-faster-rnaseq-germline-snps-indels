# Build Rust code

This repository keeps all Rust code for the RNA-seq pipeline in a single crate:

- crate root: `src/rust`
- shared Rust modules: `src/rust/src`
- binary entrypoints: `src/rust/src/bin`
- installed release binaries: `src/rust/bin`

The Python pipeline looks for the installed binaries in `src/rust/bin` by default.

## Binaries

Current Rust binaries:

- `rust_interval_tools`
- `rust_split_n_cigar_reads`
- `rust_base_recalibrator`
- `rust_apply_bqsr`
- `rust_mark_duplicates`
- `rust_hc_prefilter`

## Environment

Use the same shell environment as normal local runs:

```bash
source /data/p/anaconda3/bin/activate base
export PATH=/data/p/bin:$PATH
cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels
```

The crate already includes [src/rust/.cargo/config.toml](/data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/src/rust/.cargo/config.toml:1), which sets:

- `LIBCLANG_PATH`
- `LD_LIBRARY_PATH`
- `BINDGEN_EXTRA_CLANG_ARGS`

Those settings are needed so `rust-htslib` / `hts-sys` bindgen can find the conda `libclang` and sysroot headers in this environment.

## Build release binaries

Standard release build:

```bash
source /data/p/anaconda3/bin/activate base
export PATH=/data/p/bin:$PATH
cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/src/rust
cargo build --release --bins
```

Build a single binary:

```bash
cargo build --release --bin rust_interval_tools
```

## Build with static CRT

In this environment, the fully static `musl` target is not installed, so this repo currently uses the GNU target with static CRT enabled.

From `src/rust`:

```bash
source /data/p/anaconda3/bin/activate base
export PATH=/data/p/bin:$PATH
cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/src/rust
RUSTFLAGS='-C target-feature=+crt-static' \
  cargo build --release --bins --target x86_64-unknown-linux-gnu
```

That produces `static-pie` binaries under:

```text
src/rust/target/x86_64-unknown-linux-gnu/release/
```

## Install binaries to `src/rust/bin`

After a successful static-CRT build, copy the binaries into the location used by the Python pipeline:

```bash
source /data/p/anaconda3/bin/activate base
export PATH=/data/p/bin:$PATH
cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels
install -m 755 src/rust/target/x86_64-unknown-linux-gnu/release/rust_interval_tools src/rust/bin/
install -m 755 src/rust/target/x86_64-unknown-linux-gnu/release/rust_split_n_cigar_reads src/rust/bin/
install -m 755 src/rust/target/x86_64-unknown-linux-gnu/release/rust_base_recalibrator src/rust/bin/
install -m 755 src/rust/target/x86_64-unknown-linux-gnu/release/rust_apply_bqsr src/rust/bin/
install -m 755 src/rust/target/x86_64-unknown-linux-gnu/release/rust_mark_duplicates src/rust/bin/
install -m 755 src/rust/target/x86_64-unknown-linux-gnu/release/rust_hc_prefilter src/rust/bin/
```

Or in one loop:

```bash
source /data/p/anaconda3/bin/activate base
export PATH=/data/p/bin:$PATH
cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels
for bin in \
  rust_interval_tools \
  rust_split_n_cigar_reads \
  rust_base_recalibrator \
  rust_apply_bqsr \
  rust_mark_duplicates \
  rust_hc_prefilter
do
  install -m 755 "src/rust/target/x86_64-unknown-linux-gnu/release/$bin" src/rust/bin/
done
```

## Verify the build

Run tests:

```bash
source /data/p/anaconda3/bin/activate base
export PATH=/data/p/bin:$PATH
cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/src/rust
cargo test --release --bins
```

Check that the installed binaries are present:

```bash
cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels
ls -lah src/rust/bin
```

Check that a binary is statically linked:

```bash
file src/rust/bin/rust_interval_tools
ldd src/rust/bin/rust_interval_tools
```

Expected result:

- `file` shows `static-pie linked`
- `ldd` shows `statically linked`

Check wrapper and pipeline syntax:

```bash
source /data/p/anaconda3/bin/activate base
export PATH=/data/p/bin:$PATH
cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels
python -m py_compile src/*.py
python src/run_pipeline.py --help
python src/step_bed_to_interval_list.py --help
python src/step_split_n_cigar.py --help
python src/step_base_recalibrator.py --help
python src/step_apply_bqsr.py --help
python src/step_mark_duplicates.py --help
python src/step_split_hc_intervals.py --help
python src/step_hc_prefilter.py --help
```

## How the pipeline uses the binaries

By default, the Python code resolves Rust tools from `src/rust/bin`.

Relevant files:

- [src/step_common.py](/data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/src/step_common.py:32)
- [src/run_pipeline.py](/data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/src/run_pipeline.py:159)
- [src/run_pipeline.py](/data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/src/run_pipeline.py:278)

Useful runtime options:

- `--rust-bin-dir /path/to/bin_dir`
- `--no-rust`
- per-wrapper `--rust-bin /path/to/binary`

## Troubleshooting

If bindgen or `hts-sys` fails during build:

1. Make sure conda `base` is activated.
2. Build from `src/rust`, so the local `.cargo/config.toml` is picked up.
3. Confirm these files exist:
   - `/data/p/anaconda3/lib/libclang.so`
   - `/data/p/anaconda3/x86_64-conda-linux-gnu/sysroot`
4. Retry with a clean build:

```bash
source /data/p/anaconda3/bin/activate base
export PATH=/data/p/bin:$PATH
cd /data/p/gatk/gatk-faster-rnaseq-germline-snps-indels/src/rust
cargo clean
cargo build --release --bins
```

If the GNU static-CRT build fails:

1. check that Rust supports `x86_64-unknown-linux-gnu` on this machine:

```bash
rustc --print target-list | rg x86_64-unknown-linux-gnu
```

2. retry without static CRT to confirm the crate itself still builds:

```bash
cargo build --release --bins
```

If you need to place build outputs elsewhere, use Cargo's target directory support:

```bash
CARGO_TARGET_DIR=/some/fast/storage/path cargo build --release --bins
```

# Aim

Develop a fast pipeline for RNA-seq germline SNP and indel calling.

- Use Python for pipeline orchestration.
- Use `argparse`.
- Accept paired FASTQ input or an aligned STAR BAM.
- Produce a VCF.

Read `.local/AGENTS.md` if this file exists.

# Project Rules

- Keep the pipeline structure close to the Broad/GATK RNA-seq germline workflow, then replace slow steps incrementally where that improves runtime.
- Python package code should live under `src/gatk_faster_rnaseq/`.
- The preferred main entrypoint is the repo-root script `gatk-faster-rnaseq-germline-snps-indels.py`.
- Do not add top-level pipeline wrapper files back under `src/` such as `src/run_pipeline.py`, `src/step_*.py`, or `src/step_common.py`.
- Rust implementation code should live in `src/rust/src` library modules.
- Rust executable entrypoints should live in `src/rust/cli/*.rs` and stay thin.
- Installed release Rust executables used by the Python pipeline should live in repo-root `rust_binary/`.
- Keep new Python and Rust source files below 1000 lines when practical. If a file is approaching that size, split it into submodules before adding more logic.


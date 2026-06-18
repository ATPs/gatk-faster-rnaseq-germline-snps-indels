# Aim
develop a ultra fast pipeline for using RNAseq data to call germline snps and indels. Use python to run the pipeline. use argparse. input is fastq or aligned bam from STAR. output is vcf.

# how to do
To do this, we need to change slow steps in GATK pipeline to faster ones. we may write some new tools in rust to replace slow tools. rust codes should store in folder `gatk-faster-rnaseq-germline-snps-indels/src`.

# software
gatk: /data/p/gatk/gatk-4.6.2.0
STAR: /data/p/star/STAR_2.7.11b/Linux_x86_64_static/STAR
samtools: /data/p/samtools/samtools-1.22.1_installed/bin/samtools
sambamba: /data/p/tools/sambamba/bin/sambamba-1.0.1
python: /data/p/anaconda3/bin/python
bwa: /data/p/tools/bwa/bwa
rustc: /data/p/sys/rust/1.96.0/bin/rustc
cargo: /data/p/sys/rust/1.96.0/bin/cargo

the source code of gatk is in folder `/data/p/gatk/gatk`. we may need to understand the code to write some rust codes.

# development
- `gatk-faster-rnaseq-germline-snps-indels/development/background.md` stores the background knowledge. need to check this file to understand the pipeline.
- developement notes: each time after changing the code, write a markdown file in folder  `gatk-faster-rnaseq-germline-snps-indels/development` with name like `2023-01-01.{title}.md` to record the aim, changes and conclusion. here "title" means the title of the aim.
- during testing running the code, if the job run for a long time, check status every 10 minutes.
- sambamba markdup 比 Picard/GATK MarkDuplicates 更适合多线程
- python package code should live under `src/gatk_faster_rnaseq/`. the preferred main entrypoint is repo root `gatk-faster-rnaseq-germline-snps-indels.py`.
- do not add top-level pipeline wrapper files back under `src/` such as `src/run_pipeline.py`, `src/step_*.py`, or `src/step_common.py`. use the package modules directly.
- rust implementation code should live in `src/rust/src` library modules. rust executable entrypoints should live in `src/rust/cli/*.rs` and stay thin.
- keep new python and rust source files below 1000 lines when practical. if a file is approaching that size, split it into submodules before adding more logic.

# test the code
use the `/XCLabServer002_fastIO` folder to save large files. 
use this RNA-seq data: /data1/xlab/researches/AML/Leucegene/raw/PRJNA214592/SRR949115
use this reference genome: /data1/pub/gatk/broad_hg38/Homo_sapiens_assembly38.fasta
gatk files in folder: /data1/pub/gatk/broad_hg38
use up to 40 threads

# rust
build follow the build.md file. update this file if needed.

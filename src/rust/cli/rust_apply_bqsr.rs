use anyhow::Result;
use clap::Parser;
use gatk_faster_rnaseq_rust::{apply_bqsr, ApplyBqsrConfig};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Apply GATK BQSR recalibration tables to a BAM"
)]
struct Args {
    #[arg(short = 'R', long = "ref")]
    _reference: Option<PathBuf>,

    #[arg(short = 'I', long = "input-bam")]
    input_bam: PathBuf,

    #[arg(long = "input-table", alias = "bqsr-recal-file")]
    recal_table: PathBuf,

    #[arg(short = 'O', long = "output-bam")]
    output_bam: PathBuf,

    #[arg(long = "output-index")]
    output_index: Option<PathBuf>,

    #[arg(long, default_value_t = 4)]
    threads: usize,

    #[arg(long = "use-original-qualities")]
    use_original_qualities: bool,

    #[arg(long = "allow-missing-read-group")]
    allow_missing_read_groups: bool,

    #[arg(long = "use-report-quantization")]
    use_report_quantization: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let stats = apply_bqsr(&ApplyBqsrConfig {
        input_bam: args.input_bam,
        recal_table: args.recal_table,
        output_bam: args.output_bam.clone(),
        output_index: args.output_index,
        threads: args.threads,
        use_original_qualities: args.use_original_qualities,
        allow_missing_read_groups: args.allow_missing_read_groups,
        use_report_quantization: args.use_report_quantization,
    })?;
    eprintln!(
        "processed {} records and {} bases",
        stats.records, stats.bases
    );
    println!("{}", args.output_bam.display());
    Ok(())
}

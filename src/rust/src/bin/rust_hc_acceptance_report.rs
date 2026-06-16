use anyhow::Result;
use clap::Parser;
use gatk_faster_rnaseq_rust::hc_tools::write_acceptance_report;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Combine HaplotypeCaller matching reports into one acceptance report"
)]
struct Args {
    #[arg(long)]
    output_md: PathBuf,

    #[arg(long, default_value = "HaplotypeCaller Rust acceptance report")]
    title: String,

    #[arg(required = true)]
    input_reports: Vec<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    write_acceptance_report(&args.output_md, &args.input_reports, &args.title)?;
    println!("{}", args.output_md.display());
    Ok(())
}

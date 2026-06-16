use anyhow::Result;
use clap::Parser;
use gatk_faster_rnaseq_rust::hc_tools::{compare_vcfs, write_vcf_comparison};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Compare two HaplotypeCaller VCFs by exact allele key"
)]
struct Args {
    #[arg(long)]
    a_vcf: PathBuf,

    #[arg(long)]
    b_vcf: PathBuf,

    #[arg(long, default_value = "A")]
    a_label: String,

    #[arg(long, default_value = "B")]
    b_label: String,

    #[arg(long)]
    output_prefix: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let comparison = compare_vcfs(&args.a_vcf, &args.b_vcf, &args.a_label, &args.b_label)?;
    write_vcf_comparison(&args.output_prefix, &comparison)?;
    println!(
        "{}\t{}\tshared_pass={}\ta_private_pass={}\tb_private_pass={}\ta_sensitivity={:.3}\tb_precision={:.3}",
        args.a_label,
        args.b_label,
        comparison.pass_records.shared,
        comparison.pass_records.a_private,
        comparison.pass_records.b_private,
        comparison.pass_records.a_sensitivity,
        comparison.pass_records.b_precision_vs_a
    );
    Ok(())
}

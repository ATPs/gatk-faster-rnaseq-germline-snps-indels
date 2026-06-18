use anyhow::Result;
use clap::Parser;
use gatk_faster_rnaseq_rust::hc_tools::{select_regions, write_selected_regions};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Select focused HaplotypeCaller debug regions from two VCFs"
)]
struct Args {
    #[arg(long)]
    a_vcf: PathBuf,

    #[arg(long)]
    b_vcf: PathBuf,

    #[arg(long)]
    output_prefix: PathBuf,

    #[arg(long, default_value_t = 100)]
    padding: u64,

    #[arg(long, default_value_t = 100)]
    max_per_category: usize,

    #[arg(long, default_value_t = true)]
    pass_only: bool,

    #[arg(long)]
    interval_list_template: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let rows = select_regions(
        &args.a_vcf,
        &args.b_vcf,
        args.padding,
        args.max_per_category,
        args.pass_only,
    )?;
    write_selected_regions(
        &args.output_prefix,
        &rows,
        args.interval_list_template.as_deref(),
    )?;
    println!("selected_regions={}", rows.len());
    Ok(())
}

use anyhow::Result;
use clap::Parser;
use gatk_faster_rnaseq_rust::hc_tools::write_vcf_genotype_table;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(author, version, about = "Convert a VCF to a genotype stage TSV")]
struct Args {
    #[arg(long)]
    input_vcf: PathBuf,

    #[arg(long)]
    output_tsv: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let rows = write_vcf_genotype_table(&args.input_vcf, &args.output_tsv)?;
    println!("rows={rows}");
    Ok(())
}

use anyhow::Result;
use clap::Parser;
use gatk_faster_rnaseq_rust::hc_tools::{extract_gatk_debug_tables, GatkDebugExtractConfig};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Extract structured TSV tables from GATK HaplotypeCaller debug output"
)]
struct Args {
    #[arg(long)]
    genotyper_debug: Option<PathBuf>,

    #[arg(long)]
    assembly_state: Option<PathBuf>,

    #[arg(long)]
    output_prefix: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let stats = extract_gatk_debug_tables(&GatkDebugExtractConfig {
        genotyper_debug: args.genotyper_debug,
        assembly_state: args.assembly_state,
        output_prefix: args.output_prefix,
    })?;
    println!(
        "genotyper_haplotypes={} pairhmm_scores={} event_allele_links={} allele_likelihoods={} read_quality_rows={} assembly_reads={} assembly_haplotypes={}",
        stats.genotyper_haplotypes,
        stats.pairhmm_scores,
        stats.event_allele_links,
        stats.allele_likelihoods,
        stats.read_quality_rows,
        stats.assembly_reads,
        stats.assembly_haplotypes
    );
    Ok(())
}

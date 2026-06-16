use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use gatk_faster_rnaseq_rust::hc::activity::discover_active_regions;
use gatk_faster_rnaseq_rust::hc::args::{
    ActiveRegionDiscoveryConfig, HaplotypeCallerConfig,
};
use gatk_faster_rnaseq_rust::hc::call_variants;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(author, version, about = "Rust HaplotypeCaller rewrite workbench")]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Run the full HaplotypeCaller rewrite. The assembly/genotyping core is under development.
    Call(CallArgs),

    /// Discover HaplotypeCaller-style candidate active regions from a BAM.
    DiscoverActiveRegions(DiscoverActiveRegionsArgs),
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
enum PairHmmImplementationArg {
    Rust,
    Native,
}

impl PairHmmImplementationArg {
    fn as_str(self) -> &'static str {
        match self {
            PairHmmImplementationArg::Rust => "rust",
            PairHmmImplementationArg::Native => "native",
        }
    }
}

#[derive(Debug, Parser)]
struct CallArgs {
    #[arg(short = 'I', long)]
    input_bam: PathBuf,

    #[arg(short = 'R', long = "ref")]
    reference: PathBuf,

    #[arg(short = 'L', long)]
    input_interval_list: PathBuf,

    #[arg(short = 'O', long)]
    output_vcf: PathBuf,

    #[arg(long)]
    dbsnp: Option<PathBuf>,

    #[arg(long, default_value_t = false)]
    dont_use_soft_clipped_bases: bool,

    #[arg(long, default_value_t = 20.0)]
    standard_min_confidence_threshold_for_calling: f64,

    #[arg(long, default_value_t = 40)]
    threads: usize,

    #[arg(long, default_value_t = 128)]
    memory_gb: usize,

    #[arg(long, default_value_t = 8)]
    native_pair_hmm_threads: usize,

    #[arg(long, value_enum, default_value_t = PairHmmImplementationArg::Native)]
    pair_hmm_implementation: PairHmmImplementationArg,
}

#[derive(Debug, Parser)]
struct DiscoverActiveRegionsArgs {
    #[arg(long)]
    input_bam: PathBuf,

    #[arg(long = "ref")]
    reference: PathBuf,

    #[arg(long)]
    input_interval_list: PathBuf,

    #[arg(long)]
    output_active_bed: PathBuf,

    #[arg(long)]
    output_summary: PathBuf,

    #[arg(long, default_value_t = 20)]
    min_mapq: u8,

    #[arg(long, default_value_t = 10)]
    min_baseq: u8,

    #[arg(long, default_value_t = 1)]
    min_alt_count: u32,

    #[arg(long, default_value_t = 1)]
    min_indel_count: u32,

    #[arg(long, default_value_t = 0.0)]
    min_alt_fraction: f64,

    #[arg(long, default_value_t = 150)]
    active_region_padding: u64,

    #[arg(long, default_value_t = 100_000)]
    max_depth: u32,

    #[arg(long, default_value_t = 1)]
    threads: usize,

    #[arg(long, default_value_t = false)]
    exclude_supplementary: bool,
}

fn main() -> Result<()> {
    match Args::parse().command {
        Command::Call(args) => call_variants(&HaplotypeCallerConfig {
            input_bam: args.input_bam,
            reference: args.reference,
            input_interval_list: args.input_interval_list,
            output_vcf: args.output_vcf,
            dbsnp: args.dbsnp,
            dont_use_soft_clipped_bases: args.dont_use_soft_clipped_bases,
            standard_min_confidence_threshold_for_calling: args
                .standard_min_confidence_threshold_for_calling,
            threads: args.threads,
            memory_gb: args.memory_gb,
            native_pair_hmm_threads: args.native_pair_hmm_threads,
            pair_hmm_implementation: args.pair_hmm_implementation.as_str().to_string(),
        }),
        Command::DiscoverActiveRegions(args) => {
            let summary = discover_active_regions(&ActiveRegionDiscoveryConfig {
                input_bam: args.input_bam,
                reference: args.reference,
                input_interval_list: args.input_interval_list,
                output_active_bed: args.output_active_bed,
                output_summary: args.output_summary,
                min_mapq: args.min_mapq,
                min_baseq: args.min_baseq,
                min_alt_count: args.min_alt_count,
                min_indel_count: args.min_indel_count,
                min_alt_fraction: args.min_alt_fraction,
                active_region_padding: args.active_region_padding,
                max_depth: args.max_depth,
                threads: args.threads,
                exclude_supplementary: args.exclude_supplementary,
            })?;
            println!(
                "active_regions\t{}\nactive_bases\t{}",
                summary.active_regions, summary.active_bases
            );
            Ok(())
        }
    }
}

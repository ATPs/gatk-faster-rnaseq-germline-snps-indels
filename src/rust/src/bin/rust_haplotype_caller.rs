use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use gatk_faster_rnaseq_rust::haplotype_caller::{call_variants, HaplotypeCallerConfig};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(author, version, about = "Rust HaplotypeCaller rewrite workbench")]
struct Args {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Run the pipeline-only Rust HaplotypeCaller replacement.
    Call(CallArgs),
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
    exclude_supplementary: bool,

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

fn main() -> Result<()> {
    match Args::parse().command {
        Command::Call(args) => call_variants(&HaplotypeCallerConfig {
            input_bam: args.input_bam,
            reference: args.reference,
            input_interval_list: args.input_interval_list,
            output_vcf: args.output_vcf,
            dbsnp: args.dbsnp,
            exclude_supplementary: args.exclude_supplementary,
            dont_use_soft_clipped_bases: args.dont_use_soft_clipped_bases,
            standard_min_confidence_threshold_for_calling: args
                .standard_min_confidence_threshold_for_calling,
            threads: args.threads,
            memory_gb: args.memory_gb,
            native_pair_hmm_threads: args.native_pair_hmm_threads,
            pair_hmm_implementation: args.pair_hmm_implementation.as_str().to_string(),
        }),
    }
}

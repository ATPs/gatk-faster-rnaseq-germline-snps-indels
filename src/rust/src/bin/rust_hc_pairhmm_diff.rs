use anyhow::Result;
use clap::Parser;
use gatk_faster_rnaseq_rust::hc_tools::{run_stage_diff, StageDiffConfig};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(author, version, about = "Diff Java/Rust PairHMM likelihood tables")]
struct Args {
    #[arg(long)]
    java: PathBuf,
    #[arg(long)]
    rust: PathBuf,
    #[arg(long)]
    output_prefix: PathBuf,
    #[arg(long, default_value = "region,read,haplotype")]
    key_columns: String,
    #[arg(long, default_value_t = 1e-3)]
    numeric_tolerance: f64,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let summary = run_stage_diff(&StageDiffConfig {
        java_path: args.java,
        rust_path: args.rust,
        key_columns: split_key_columns(&args.key_columns),
        numeric_tolerance: args.numeric_tolerance,
        output_prefix: args.output_prefix,
        stage_name: "pairhmm".to_string(),
    })?;
    println!(
        "shared_rows={} java_private={} rust_private={} field_diffs={}",
        summary.shared_rows,
        summary.java_private_rows,
        summary.rust_private_rows,
        summary.field_diffs
    );
    Ok(())
}

fn split_key_columns(value: &str) -> Vec<String> {
    value
        .split(',')
        .map(|field| field.trim().to_string())
        .filter(|field| !field.is_empty())
        .collect()
}

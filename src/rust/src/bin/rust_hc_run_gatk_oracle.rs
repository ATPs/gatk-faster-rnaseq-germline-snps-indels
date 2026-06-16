use anyhow::{Context, Result};
use clap::Parser;
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Run Java GATK HaplotypeCaller with stage debug/oracle outputs"
)]
struct Args {
    #[arg(long, default_value = "/data/p/gatk/gatk-4.6.2.0/gatk")]
    gatk: PathBuf,

    #[arg(long, default_value = "24g")]
    java_mem: String,

    #[arg(short = 'I', long)]
    input_bam: PathBuf,

    #[arg(short = 'R', long = "ref")]
    reference: PathBuf,

    #[arg(short = 'L', long)]
    input_interval_list: PathBuf,

    #[arg(long)]
    dbsnp: Option<PathBuf>,

    #[arg(short = 'O', long)]
    output_vcf: PathBuf,

    #[arg(long)]
    output_dir: PathBuf,

    #[arg(long, default_value_t = 8)]
    pair_hmm_threads: usize,

    #[arg(long, default_value_t = true)]
    debug_assembly: bool,

    #[arg(long, default_value_t = false)]
    emit_assembly_events_vcf: bool,

    #[arg(long = "extra-gatk-arg")]
    extra_gatk_args: Vec<String>,

    #[arg(long, default_value_t = false)]
    dry_run: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();
    fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("creating {}", args.output_dir.display()))?;
    if let Some(parent) = args.output_vcf.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).with_context(|| format!("creating {}", parent.display()))?;
        }
    }

    let assembly_state = args.output_dir.join("gatk_hc_oracle.assembly_state.txt");
    let pairhmm = args.output_dir.join("gatk_hc_oracle.pairhmm.txt");
    let genotyper = args.output_dir.join("gatk_hc_oracle.genotyper.txt");

    let mut command = vec![
        args.gatk.display().to_string(),
        "--java-options".to_string(),
        format!("-Xmx{}", args.java_mem),
        "HaplotypeCaller".to_string(),
        "-R".to_string(),
        args.reference.display().to_string(),
        "-I".to_string(),
        args.input_bam.display().to_string(),
        "-L".to_string(),
        args.input_interval_list.display().to_string(),
        "-O".to_string(),
        args.output_vcf.display().to_string(),
        "--dont-use-soft-clipped-bases".to_string(),
        "--standard-min-confidence-threshold-for-calling".to_string(),
        "20".to_string(),
        "--native-pair-hmm-threads".to_string(),
        args.pair_hmm_threads.to_string(),
        "--debug-assembly-region-state".to_string(),
        assembly_state.display().to_string(),
        "--pair-hmm-results-file".to_string(),
        pairhmm.display().to_string(),
        "--debug-genotyper-output".to_string(),
        genotyper.display().to_string(),
    ];
    if args.debug_assembly {
        command.push("--debug-assembly".to_string());
    }
    if args.emit_assembly_events_vcf {
        let assembly_variants = args
            .output_dir
            .join("gatk_hc_oracle.assembly_events.vcf.gz");
        command.push("--debug-assembly-variants-out".to_string());
        command.push(assembly_variants.display().to_string());
    }
    if let Some(dbsnp) = args.dbsnp {
        command.push("--dbsnp".to_string());
        command.push(dbsnp.display().to_string());
    }
    command.extend(args.extra_gatk_args);

    let command_path = args.output_dir.join("gatk_hc_oracle.command.txt");
    let mut command_file = File::create(&command_path)
        .with_context(|| format!("creating {}", command_path.display()))?;
    writeln!(command_file, "{}", shell_quote_command(&command))?;
    println!("{}", shell_quote_command(&command));

    if args.dry_run {
        return Ok(());
    }

    let status = Command::new(&command[0])
        .args(&command[1..])
        .status()
        .context("running GATK HaplotypeCaller oracle command")?;
    if !status.success() {
        anyhow::bail!("GATK oracle command failed with status {status}");
    }
    Ok(())
}

fn shell_quote_command(command: &[String]) -> String {
    command
        .iter()
        .map(|part| {
            if part
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || "-_./:=+".contains(c))
            {
                part.clone()
            } else {
                format!("'{}'", part.replace('\'', "'\\''"))
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

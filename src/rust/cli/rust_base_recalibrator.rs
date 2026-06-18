use anyhow::Result;

fn main() -> Result<()> {
    gatk_faster_rnaseq_rust::base_recalibrator::run_cli()
}

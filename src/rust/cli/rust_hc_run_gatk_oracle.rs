use anyhow::Result;

fn main() -> Result<()> {
    gatk_faster_rnaseq_rust::gatk_oracle::run_cli()
}

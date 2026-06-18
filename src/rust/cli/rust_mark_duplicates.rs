use anyhow::Result;

fn main() -> Result<()> {
    gatk_faster_rnaseq_rust::mark_duplicates::run_cli()
}

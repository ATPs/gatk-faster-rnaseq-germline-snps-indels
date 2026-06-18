use anyhow::Result;

fn main() -> Result<()> {
    gatk_faster_rnaseq_rust::hc_prefilter::run_cli()
}

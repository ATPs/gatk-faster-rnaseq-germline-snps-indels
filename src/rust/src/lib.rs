pub mod apply_bqsr;
pub mod assembly;
pub mod base_recalibrator;
pub mod gatk_oracle;
pub mod haplotype_caller;
pub mod hc_prefilter;
pub mod hc_tools;
pub mod interval_tools;
pub mod mark_duplicates;
pub mod pair_hmm;
pub mod smith_waterman;
pub mod split_n_cigar;

pub use apply_bqsr::{apply_bqsr, bam_index_path, ApplyBqsrConfig, ApplyBqsrStats};

pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

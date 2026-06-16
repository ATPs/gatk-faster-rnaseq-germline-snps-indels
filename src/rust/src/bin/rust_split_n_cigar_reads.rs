use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use gatk_faster_rnaseq_rust::bam_index_path;
use gatk_faster_rnaseq_rust::split_n_cigar::{
    transform_record_for_compatibility_manager, transform_record_with_contig_names,
    OverhangFixingManager, OverhangOptions, PendingRecord, SortKey, SortQueue, SplitMode,
    SplitOptions, SplitStats,
};
use rust_htslib::bam::index::{self, Type};
use rust_htslib::bam::{self, Format, Header, Read, Writer};
use std::cmp::Reverse;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Rust SplitNCigarReads replacement for RNA-seq BAM preprocessing"
)]
struct Args {
    #[arg(short = 'I', long = "input-bam")]
    input_bam: PathBuf,

    #[arg(short = 'O', long = "output-bam")]
    output_bam: PathBuf,

    #[arg(short = 'R', long = "reference")]
    reference: Option<PathBuf>,

    #[arg(long = "threads", default_value_t = 1)]
    threads: usize,

    #[arg(long = "skip-mapping-quality-transform")]
    skip_mapping_quality_transform: bool,

    #[arg(long = "process-secondary-alignments")]
    process_secondary_alignments: bool,

    #[arg(long = "mode", value_enum, default_value_t = ModeArg::Fast)]
    mode: ModeArg,

    #[arg(long = "max-reads-in-memory", default_value_t = 150_000)]
    max_reads_in_memory: usize,

    #[arg(long = "max-mismatches-in-overhang", default_value_t = 1)]
    max_mismatches_in_overhang: usize,

    #[arg(long = "max-bases-in-overhang", default_value_t = 40)]
    max_bases_in_overhang: usize,

    #[arg(long = "do-not-fix-overhangs")]
    do_not_fix_overhangs: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ModeArg {
    Fast,
    Compatibility,
}

fn main() -> Result<()> {
    let args = Args::parse();
    run(args)
}

fn run(args: Args) -> Result<()> {
    if args.threads == 0 {
        anyhow::bail!("--threads must be at least 1");
    }
    if args.max_reads_in_memory == 0 {
        anyhow::bail!("--max-reads-in-memory must be at least 1");
    }
    if matches!(args.mode, ModeArg::Compatibility)
        && !args.do_not_fix_overhangs
        && args.reference.is_none()
    {
        anyhow::bail!(
            "--reference is required for --mode compatibility unless --do-not-fix-overhangs is set"
        );
    }

    let mut reader = bam::Reader::from_path(&args.input_bam)
        .with_context(|| format!("failed to open input BAM {}", args.input_bam.display()))?;
    if args.threads > 1 {
        reader
            .set_threads(args.threads)
            .context("failed to configure BAM reader threads")?;
    }

    let header = Header::from_template(reader.header());
    let mut writer = Writer::from_path(&args.output_bam, &header, Format::Bam)
        .with_context(|| format!("failed to create output BAM {}", args.output_bam.display()))?;
    if args.threads > 1 {
        writer
            .set_threads(args.threads)
            .context("failed to configure BAM writer threads")?;
    }

    let options = SplitOptions {
        skip_mq_transform: args.skip_mapping_quality_transform,
        process_secondary_alignments: args.process_secondary_alignments,
        mode: match args.mode {
            ModeArg::Fast => SplitMode::Fast,
            ModeArg::Compatibility => SplitMode::Compatibility,
        },
    };
    let contig_names = reader
        .header()
        .target_names()
        .iter()
        .map(|name| String::from_utf8_lossy(name).to_string())
        .collect::<Vec<_>>();
    let stats = if options.mode == SplitMode::Compatibility && !args.do_not_fix_overhangs {
        let mut manager = OverhangFixingManager::from_reference_path(
            args.reference
                .as_ref()
                .expect("validated compatibility reference"),
            &contig_names,
            OverhangOptions {
                max_records_in_memory: args.max_reads_in_memory,
                max_mismatches_in_overhang: args.max_mismatches_in_overhang,
                max_bases_in_overhang: args.max_bases_in_overhang,
                do_not_fix_overhangs: args.do_not_fix_overhangs,
                process_secondary_alignments: args.process_secondary_alignments,
            },
        )?;
        split_bam_compatibility(
            &args.input_bam,
            &mut reader,
            &mut writer,
            options,
            &mut manager,
            args.threads,
        )?
    } else {
        split_bam(&mut reader, &mut writer, options, &contig_names)?
    };
    drop(writer);

    let output_index = bam_index_path(&args.output_bam);
    index::build(
        &args.output_bam,
        Some(&output_index),
        Type::Bai,
        args.threads as u32,
    )
    .with_context(|| format!("failed to build BAI for {}", output_index.display()))?;

    eprintln!(
        "rust_split_n_cigar_reads\tinput_records={}\toutput_records={}\tsplit_records={}\tmq_transformed_records={}",
        stats.input_records,
        stats.output_records,
        stats.split_records,
        stats.mq_transformed_records,
    );
    println!("{}", args.output_bam.display());
    Ok(())
}

fn split_bam_compatibility<R: Read>(
    input_bam: &PathBuf,
    reader: &mut R,
    writer: &mut Writer,
    options: SplitOptions,
    manager: &mut OverhangFixingManager,
    threads: usize,
) -> Result<SplitStats> {
    let mut stats = SplitStats::default();
    let mut output = Vec::new();

    for result in reader.records() {
        let record = result.context("failed to read BAM record")?;
        stats.input_records += 1;
        let mapq_will_change = !options.skip_mq_transform && record.mapq() == 255;
        let input_will_split = !record.is_unmapped()
            && !(record.is_secondary() && !options.process_secondary_alignments)
            && record
                .cigar()
                .iter()
                .any(|element| matches!(element, bam::record::Cigar::RefSkip(_)));
        let ready = transform_record_for_compatibility_manager(&record, options, manager)?;

        if mapq_will_change {
            stats.mq_transformed_records += 1;
        }
        if input_will_split {
            stats.split_records += 1;
        }
        for record in ready {
            output.push(record);
        }
    }

    for record in manager.activate_writing() {
        output.push(record);
    }

    let mut second_reader = bam::Reader::from_path(input_bam)
        .with_context(|| format!("failed to reopen input BAM {}", input_bam.display()))?;
    if threads > 1 {
        second_reader
            .set_threads(threads)
            .context("failed to configure second-pass BAM reader threads")?;
    }

    for result in second_reader.records() {
        let record = result.context("failed to read BAM record")?;
        let ready = transform_record_for_compatibility_manager(&record, options, manager)?;
        for record in ready {
            output.push(record);
        }
    }

    for record in manager.flush() {
        output.push(record);
    }

    output.sort_by_key(SortKey::from_record);
    for record in output {
        stats.output_records += 1;
        writer
            .write(&record)
            .context("failed to write BAM record")?;
    }

    Ok(stats)
}

fn split_bam<R: Read>(
    reader: &mut R,
    writer: &mut Writer,
    options: SplitOptions,
    contig_names: &[String],
) -> Result<SplitStats> {
    let mut stats = SplitStats::default();
    let mut queue = SortQueue::new();
    let mut ordinal = 0u64;

    for result in reader.records() {
        let record = result.context("failed to read BAM record")?;
        stats.input_records += 1;
        let current_key = SortKey::from_record(&record);
        let mapq_will_change = !options.skip_mq_transform && record.mapq() == 255;
        let transformed = transform_record_with_contig_names(&record, options, Some(contig_names))?;

        if transformed.len() > 1 {
            stats.split_records += 1;
        }
        if mapq_will_change {
            stats.mq_transformed_records += 1;
        }

        for split in transformed {
            stats.output_records += 1;
            queue.push(Reverse(PendingRecord::new(split, ordinal)));
            ordinal += 1;
        }

        flush_ready(&mut queue, writer, current_key)?;
    }

    while let Some(Reverse(pending)) = queue.pop() {
        writer
            .write(&pending.into_record())
            .context("failed to write BAM record")?;
    }

    Ok(stats)
}

fn flush_ready(queue: &mut SortQueue, writer: &mut Writer, current_key: SortKey) -> Result<()> {
    while queue
        .peek()
        .map(|Reverse(pending)| pending.should_flush_before_or_at(current_key))
        .unwrap_or(false)
    {
        let pending = queue.pop().expect("queue peeked Some").0;
        writer
            .write(&pending.into_record())
            .context("failed to write BAM record")?;
    }
    Ok(())
}

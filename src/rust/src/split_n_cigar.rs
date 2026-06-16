use anyhow::{bail, Context, Result};
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::faidx;
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeSet, BinaryHeap, HashMap, VecDeque};
use std::convert::TryFrom;
use std::path::Path;

const SUPPLEMENTARY_FLAG: u16 = 0x800;
const TAGS_TO_REMOVE_FAST_MODE: [&[u8]; 5] = [b"NM", b"MD", b"NH", b"MC", b"SA"];
const TAGS_TO_REMOVE_COMPATIBILITY_MODE: [&[u8]; 3] = [b"NM", b"MD", b"NH"];
const DEFAULT_MAX_RECORDS_IN_MEMORY: usize = 150_000;
const DEFAULT_MAX_MISMATCHES_IN_OVERHANG: usize = 1;
const DEFAULT_MAX_BASES_IN_OVERHANG: usize = 40;
const MAX_SPLICES_TO_KEEP: usize = 1_000;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplitPlan {
    pub ref_offset: i64,
    pub cigar: CigarString,
}

#[derive(Debug, Clone, Copy)]
pub struct SplitOptions {
    pub skip_mq_transform: bool,
    pub process_secondary_alignments: bool,
    pub mode: SplitMode,
}

impl Default for SplitOptions {
    fn default() -> Self {
        Self {
            skip_mq_transform: false,
            process_secondary_alignments: false,
            mode: SplitMode::Fast,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitMode {
    Fast,
    Compatibility,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SplitStats {
    pub input_records: u64,
    pub output_records: u64,
    pub split_records: u64,
    pub mq_transformed_records: u64,
}

#[derive(Debug)]
pub struct PendingRecord {
    key: SortKey,
    ordinal: u64,
    record: bam::Record,
}

impl PendingRecord {
    pub fn new(record: bam::Record, ordinal: u64) -> Self {
        Self {
            key: SortKey::from_record(&record),
            ordinal,
            record,
        }
    }

    pub fn should_flush_before_or_at(&self, current: SortKey) -> bool {
        self.key <= current
    }

    pub fn into_record(self) -> bam::Record {
        self.record
    }
}

impl PartialEq for PendingRecord {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key && self.ordinal == other.ordinal
    }
}

impl Eq for PendingRecord {}

impl PartialOrd for PendingRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PendingRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key
            .cmp(&other.key)
            .then_with(|| self.ordinal.cmp(&other.ordinal))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SortKey {
    tid: i32,
    pos: i64,
}

impl SortKey {
    pub fn from_record(record: &bam::Record) -> Self {
        let tid = if record.tid() < 0 {
            i32::MAX
        } else {
            record.tid()
        };
        let pos = if record.pos() < 0 {
            i64::MAX
        } else {
            record.pos()
        };
        Self { tid, pos }
    }
}

pub type SortQueue = BinaryHeap<Reverse<PendingRecord>>;

#[derive(Debug, Clone, Copy)]
pub struct OverhangOptions {
    pub max_records_in_memory: usize,
    pub max_mismatches_in_overhang: usize,
    pub max_bases_in_overhang: usize,
    pub do_not_fix_overhangs: bool,
    pub process_secondary_alignments: bool,
}

impl Default for OverhangOptions {
    fn default() -> Self {
        Self {
            max_records_in_memory: DEFAULT_MAX_RECORDS_IN_MEMORY,
            max_mismatches_in_overhang: DEFAULT_MAX_MISMATCHES_IN_OVERHANG,
            max_bases_in_overhang: DEFAULT_MAX_BASES_IN_OVERHANG,
            do_not_fix_overhangs: false,
            process_secondary_alignments: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Splice {
    tid: i32,
    start: i64,
    end: i64,
    reference: Vec<u8>,
}

#[derive(Debug, Clone)]
struct ManagedRead {
    old_pos: i64,
    old_cigar: CigarString,
    record: bam::Record,
}

impl ManagedRead {
    fn new(record: bam::Record) -> Self {
        let old_pos = record.pos();
        let old_cigar = CigarString(record.cigar().iter().copied().collect());
        Self {
            old_pos,
            old_cigar,
            record,
        }
    }

    fn has_been_overhang_clipped(&self) -> bool {
        self.old_pos != self.record.pos()
            || self.old_cigar.to_string() != self.record.cigar().to_string()
    }
}

#[derive(Debug, Clone)]
struct ReadGroup {
    records: Vec<ManagedRead>,
}

impl ReadGroup {
    fn new(records: Vec<bam::Record>) -> Self {
        Self {
            records: records.into_iter().map(ManagedRead::new).collect(),
        }
    }

    fn first_key(&self) -> SortKey {
        self.records
            .first()
            .map(|read| SortKey::from_record(&read.record))
            .unwrap_or(SortKey {
                tid: i32::MAX,
                pos: i64::MAX,
            })
    }

    fn len(&self) -> usize {
        self.records.len()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct MateKey {
    name: Vec<u8>,
    first_of_pair: bool,
    mate_start: i64,
}

#[derive(Debug, Clone)]
struct MatePrediction {
    start: i64,
    cigar: String,
}

#[derive(Debug)]
pub struct OverhangFixingManager {
    reference: faidx::Reader,
    contig_names: Vec<String>,
    options: OverhangOptions,
    output_to_file: bool,
    waiting_groups: VecDeque<ReadGroup>,
    waiting_reads: usize,
    splices: BTreeSet<Splice>,
    mate_changed_reads: HashMap<MateKey, MatePrediction>,
}

impl OverhangFixingManager {
    pub fn from_reference_path<P: AsRef<Path>>(
        reference: P,
        contig_names: &[String],
        options: OverhangOptions,
    ) -> Result<Self> {
        let reference = faidx::Reader::from_path(reference.as_ref()).with_context(|| {
            format!(
                "failed to open FASTA reference {}",
                reference.as_ref().display()
            )
        })?;
        Ok(Self {
            reference,
            contig_names: contig_names.to_vec(),
            options,
            output_to_file: false,
            waiting_groups: VecDeque::new(),
            waiting_reads: 0,
            splices: BTreeSet::new(),
            mate_changed_reads: HashMap::new(),
        })
    }

    pub fn activate_writing(&mut self) -> Vec<bam::Record> {
        let records = self.flush();
        self.splices.clear();
        self.output_to_file = true;
        records
    }

    pub fn add_read_group(&mut self, records: Vec<bam::Record>) -> Result<Vec<bam::Record>> {
        if records.is_empty() {
            return Ok(Vec::new());
        }

        let too_many_reads = self.options.max_records_in_memory > 0
            && self.waiting_reads >= self.options.max_records_in_memory;
        let encountered_new_contig = self.encountered_new_contig(&records[0]);
        let target_queue_size = if encountered_new_contig {
            0
        } else {
            self.options.max_records_in_memory / 2
        };
        let mut ready = if too_many_reads || encountered_new_contig {
            self.write_reads(target_queue_size)?
        } else {
            Vec::new()
        };

        let mut group = ReadGroup::new(records);
        let splices = self.splices.iter().cloned().collect::<Vec<_>>();
        for splice in &splices {
            for read in &mut group.records {
                self.fix_split(read, splice)?;
            }
        }

        self.waiting_reads += group.len();
        self.waiting_groups.push_back(group);
        Ok(ready.drain(..).collect())
    }

    pub fn add_splice_position(&mut self, tid: i32, start: i64, end: i64) -> Result<()> {
        if self.options.do_not_fix_overhangs || tid < 0 || start > end {
            return Ok(());
        }
        if self
            .splices
            .iter()
            .any(|splice| splice.tid == tid && splice.start == start && splice.end == end)
        {
            return Ok(());
        }

        if let Some(first) = self.splices.iter().next() {
            if first.tid != tid {
                self.splices.clear();
            }
        }

        let contig = self
            .contig_names
            .get(tid as usize)
            .with_context(|| format!("missing contig name for tid {tid}"))?;
        let reference = self
            .reference
            .fetch_seq(contig, start as usize, end as usize)
            .with_context(|| {
                format!(
                    "failed to fetch reference {contig}:{}-{}",
                    start + 1,
                    end + 1
                )
            })?;
        let splice = Splice {
            tid,
            start,
            end,
            reference,
        };

        for group in &mut self.waiting_groups {
            for read in &mut group.records {
                fix_split_with_options(read, &splice, self.options)?;
            }
        }

        self.splices.insert(splice);
        if self.splices.len() > MAX_SPLICES_TO_KEEP {
            self.clean_splices();
        }
        Ok(())
    }

    pub fn set_predicted_mate_information(&self, record: &mut bam::Record) {
        if !self.output_to_file || !record.is_paired() || record.is_unmapped() {
            return;
        }
        let key = MateKey {
            name: record.qname().to_vec(),
            first_of_pair: record.is_first_in_template(),
            mate_start: record.mpos(),
        };
        if let Some(prediction) = self.mate_changed_reads.get(&key) {
            record.set_mpos(prediction.start);
            if record.aux(b"MC").is_ok() {
                let _ = record.remove_aux(b"MC");
                let _ = record.push_aux(b"MC", Aux::String(&prediction.cigar));
            }
        }
    }

    pub fn flush(&mut self) -> Vec<bam::Record> {
        self.write_reads(0).unwrap_or_default()
    }

    fn fix_split(&self, read: &mut ManagedRead, splice: &Splice) -> Result<()> {
        fix_split_with_options(read, splice, self.options)
    }

    fn encountered_new_contig(&self, first_new_record: &bam::Record) -> bool {
        let Some(first_group) = self.waiting_groups.front() else {
            return false;
        };
        let Some(top_read) = first_group.records.first() else {
            return false;
        };
        !top_read.record.is_unmapped()
            && !first_new_record.is_unmapped()
            && top_read.record.tid() != first_new_record.tid()
    }

    fn write_reads(&mut self, target_queue_size: usize) -> Result<Vec<bam::Record>> {
        let mut ready = Vec::new();
        while self.waiting_reads > target_queue_size {
            let Some(mut group) = self.pop_next_group() else {
                break;
            };
            self.waiting_reads = self.waiting_reads.saturating_sub(group.len());

            if self.output_to_file {
                let mut records = group
                    .records
                    .drain(..)
                    .map(|managed| managed.record)
                    .collect::<Vec<_>>();
                repair_sa_tags(&mut records, Some(&self.contig_names))?;
                ready.extend(records);
            } else if let Some(first) = group.records.first() {
                if !first.record.is_secondary() && first.has_been_overhang_clipped() {
                    self.mate_changed_reads.insert(
                        MateKey {
                            name: first.record.qname().to_vec(),
                            first_of_pair: !first.record.is_first_in_template(),
                            mate_start: first.old_pos,
                        },
                        MatePrediction {
                            start: first.record.pos(),
                            cigar: first.record.cigar().to_string(),
                        },
                    );
                }
            }
        }
        Ok(ready)
    }

    fn pop_next_group(&mut self) -> Option<ReadGroup> {
        let (idx, _) = self
            .waiting_groups
            .iter()
            .enumerate()
            .min_by_key(|(_, group)| group.first_key())?;
        self.waiting_groups.remove(idx)
    }

    fn clean_splices(&mut self) {
        let remove_count = self.splices.len() / 2;
        let to_remove = self
            .splices
            .iter()
            .take(remove_count)
            .cloned()
            .collect::<Vec<_>>();
        for splice in to_remove {
            self.splices.remove(&splice);
        }
    }
}

pub fn transform_record(record: &bam::Record, options: SplitOptions) -> Result<Vec<bam::Record>> {
    transform_record_with_contig_names(record, options, None)
}

pub fn transform_record_with_contig_names(
    record: &bam::Record,
    options: SplitOptions,
    contig_names: Option<&[String]>,
) -> Result<Vec<bam::Record>> {
    let mut base = record.clone();
    if !options.skip_mq_transform && base.mapq() == 255 {
        base.set_mapq(60);
    }

    if base.is_unmapped() || (base.is_secondary() && !options.process_secondary_alignments) {
        if options.mode == SplitMode::Compatibility {
            remove_stale_tags(&mut base, options.mode);
            repair_mc_tag(&mut base)?;
        }
        return Ok(vec![base]);
    }

    let cigar: Vec<Cigar> = base.cigar().iter().copied().collect();
    let plans = split_cigar_at_ref_skips(&cigar)?;
    if plans.is_empty() {
        if options.mode == SplitMode::Compatibility {
            remove_stale_tags(&mut base, options.mode);
            repair_mc_tag(&mut base)?;
        }
        return Ok(vec![base]);
    }

    let original_pos = base.pos();
    let original_flags = base.flags();
    let mut records = Vec::with_capacity(plans.len());
    for (idx, plan) in plans.iter().enumerate() {
        let mut split = base.clone();
        split.set_pos(original_pos + plan.ref_offset);
        split.set_cigar(Some(&plan.cigar));
        split.set_bin(reg2bin(
            split.pos(),
            split.pos() + reference_span(&plan.cigar.0) as i64,
        ));
        remove_stale_tags(&mut split, options.mode);
        if idx > 0 {
            split.set_flags(original_flags | SUPPLEMENTARY_FLAG);
        }
        if options.mode == SplitMode::Compatibility {
            repair_mc_tag(&mut split)?;
        }
        records.push(split);
    }
    if options.mode == SplitMode::Compatibility {
        repair_sa_tags(&mut records, contig_names)?;
    }
    Ok(records)
}

fn tags_to_remove(mode: SplitMode) -> &'static [&'static [u8]] {
    match mode {
        SplitMode::Fast => &TAGS_TO_REMOVE_FAST_MODE,
        SplitMode::Compatibility => &TAGS_TO_REMOVE_COMPATIBILITY_MODE,
    }
}

fn remove_stale_tags(record: &mut bam::Record, mode: SplitMode) {
    for tag in tags_to_remove(mode) {
        let _ = record.remove_aux(tag);
    }
}

fn repair_mc_tag(record: &mut bam::Record) -> Result<()> {
    let Ok(Aux::String(mate_cigar)) = record.aux(b"MC") else {
        return Ok(());
    };
    let mate_cigar = mate_cigar.to_owned();
    let parsed = CigarString::try_from(mate_cigar.as_str())?;
    let plans = split_cigar_at_ref_skips(&parsed.0)?;
    if let Some(first_plan) = plans.first() {
        let repaired = first_plan.cigar.to_string();
        let _ = record.remove_aux(b"MC");
        record.push_aux(b"MC", Aux::String(&repaired))?;
    }
    Ok(())
}

fn repair_sa_tags(records: &mut [bam::Record], contig_names: Option<&[String]>) -> Result<()> {
    if records.len() <= 1 {
        return Ok(());
    }

    let original_sa = records
        .iter()
        .map(|record| match record.aux(b"SA") {
            Ok(Aux::String(value)) => value.to_string(),
            _ => String::new(),
        })
        .collect::<Vec<_>>();
    let entries: Vec<String> = records
        .iter()
        .map(|record| sa_entry(record, contig_names))
        .collect();
    for (idx, record) in records.iter_mut().enumerate() {
        let mut sa = String::new();
        if idx > 0 {
            sa.push_str(&entries[0]);
        }
        sa.push_str(&original_sa[idx]);
        for (entry_idx, entry) in entries.iter().enumerate() {
            if entry_idx == idx || (idx > 0 && entry_idx == 0) {
                continue;
            }
            sa.push_str(entry);
        }
        let _ = record.remove_aux(b"SA");
        if !sa.is_empty() {
            record.push_aux(b"SA", Aux::String(&sa))?;
        }
    }
    Ok(())
}

fn sa_entry(record: &bam::Record, contig_names: Option<&[String]>) -> String {
    let contig = if record.tid() < 0 {
        "*".to_string()
    } else {
        contig_names
            .and_then(|names| names.get(record.tid() as usize))
            .cloned()
            .unwrap_or_else(|| (record.tid() + 1).to_string())
    };
    let pos = if record.pos() < 0 {
        0
    } else {
        record.pos() + 1
    };
    let strand = if record.is_reverse() { "-" } else { "+" };
    let nm = match record.aux(b"NM") {
        Ok(Aux::I8(value)) => value.to_string(),
        Ok(Aux::U8(value)) => value.to_string(),
        Ok(Aux::I16(value)) => value.to_string(),
        Ok(Aux::U16(value)) => value.to_string(),
        Ok(Aux::I32(value)) => value.to_string(),
        Ok(Aux::U32(value)) => value.to_string(),
        Ok(Aux::String(value)) => value.to_string(),
        _ => "*".to_string(),
    };
    format!(
        "{contig},{pos},{strand},{},{},{};",
        record.cigar(),
        record.mapq(),
        nm
    )
}

pub fn transform_record_for_compatibility_manager(
    record: &bam::Record,
    options: SplitOptions,
    manager: &mut OverhangFixingManager,
) -> Result<Vec<bam::Record>> {
    let mut base = record.clone();
    if !options.skip_mq_transform && base.mapq() == 255 {
        base.set_mapq(60);
    }
    manager.set_predicted_mate_information(&mut base);

    if base.is_unmapped() || (base.is_secondary() && !options.process_secondary_alignments) {
        remove_stale_tags(&mut base, options.mode);
        repair_mc_tag(&mut base)?;
        return manager.add_read_group(vec![base]);
    }

    let cigar: Vec<Cigar> = base.cigar().iter().copied().collect();
    let plans = split_cigar_at_ref_skips(&cigar)?;
    if plans.is_empty() {
        remove_stale_tags(&mut base, options.mode);
        repair_mc_tag(&mut base)?;
        return manager.add_read_group(vec![base]);
    }

    let original_pos = base.pos();
    let original_flags = base.flags();
    let mut records = Vec::with_capacity(plans.len());
    for (idx, plan) in plans.iter().enumerate() {
        let mut split = base.clone();
        split.set_pos(original_pos + plan.ref_offset);
        split.set_cigar(Some(&plan.cigar));
        update_record_bin(&mut split);
        remove_stale_tags(&mut split, options.mode);
        if idx > 0 {
            split.set_flags(original_flags | SUPPLEMENTARY_FLAG);
        }
        repair_mc_tag(&mut split)?;
        records.push(split);
    }

    for splice in splice_positions_for_plans(base.tid(), original_pos, &plans) {
        manager.add_splice_position(splice.tid, splice.start, splice.end)?;
    }

    manager.add_read_group(records)
}

pub fn split_cigar_at_ref_skips(cigar: &[Cigar]) -> Result<Vec<SplitPlan>> {
    let mut emitted_before_ref_skip = false;
    let mut first = 0usize;
    let mut section_has_real_op = false;
    let mut ranges: Vec<(usize, usize)> = Vec::new();

    for (idx, element) in cigar.iter().copied().enumerate() {
        if is_real_alignment_op(element) {
            section_has_real_op = true;
        }

        if matches!(element, Cigar::RefSkip(_)) {
            if section_has_real_op {
                ranges.push((first, idx));
                emitted_before_ref_skip = true;
            }
            first = idx + 1;
            section_has_real_op = false;
        }
    }

    if !emitted_before_ref_skip {
        return Ok(Vec::new());
    }
    if first < cigar.len() && section_has_real_op {
        ranges.push((first, cigar.len()));
    }

    ranges
        .into_iter()
        .map(|(start, end)| build_split_plan(cigar, start, end))
        .collect()
}

fn build_split_plan(cigar: &[Cigar], start: usize, end: usize) -> Result<SplitPlan> {
    let mut first = start;
    let mut second = end;

    while first < second && matches!(cigar[first], Cigar::Del(_)) {
        first += 1;
    }
    while second > first && matches!(cigar[second - 1], Cigar::Del(_)) {
        second -= 1;
    }
    if first >= second {
        bail!("cannot split read section with no aligned read bases");
    }

    let ref_offset = cigar[..first]
        .iter()
        .copied()
        .map(reference_len)
        .sum::<u32>() as i64;
    let leading_soft = cigar[..first].iter().copied().map(read_len).sum::<u32>();
    let trailing_soft = cigar[second..].iter().copied().map(read_len).sum::<u32>();
    let leading_hard = leading_hard_clip_len(cigar);
    let trailing_hard = trailing_hard_clip_len(cigar);

    let mut split = Vec::with_capacity(cigar.len() + 4);
    push_cigar(&mut split, Cigar::HardClip(leading_hard));
    push_cigar(&mut split, Cigar::SoftClip(leading_soft));
    for element in cigar[first..second].iter().copied() {
        if !matches!(element, Cigar::HardClip(_)) {
            push_cigar(&mut split, element);
        }
    }
    push_cigar(&mut split, Cigar::SoftClip(trailing_soft));
    push_cigar(&mut split, Cigar::HardClip(trailing_hard));

    Ok(SplitPlan {
        ref_offset,
        cigar: CigarString(split),
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SplicePosition {
    tid: i32,
    start: i64,
    end: i64,
}

fn splice_positions_for_plans(
    tid: i32,
    original_pos: i64,
    plans: &[SplitPlan],
) -> Vec<SplicePosition> {
    plans
        .windows(2)
        .filter_map(|window| {
            let left = &window[0];
            let right = &window[1];
            let left_end = original_pos + left.ref_offset + reference_span(&left.cigar.0) as i64;
            let right_start = original_pos + right.ref_offset;
            if right_start > left_end {
                Some(SplicePosition {
                    tid,
                    start: left_end,
                    end: right_start - 1,
                })
            } else {
                None
            }
        })
        .collect()
}

fn fix_split_with_options(
    read: &mut ManagedRead,
    splice: &Splice,
    options: OverhangOptions,
) -> Result<()> {
    if read.record.is_unmapped() || read.record.tid() != splice.tid {
        return Ok(());
    }
    if !options.process_secondary_alignments && read.record.is_secondary() {
        return Ok(());
    }

    let Some((soft_start, soft_end)) = soft_reference_bounds(&read.record) else {
        return Ok(());
    };
    if soft_end < splice.start || soft_start > splice.end {
        return Ok(());
    }

    let read_bases_len = unclipped_read_bases_len(&read.record);
    let read_bases = read.record.seq().as_bytes();
    if is_left_overhang(soft_start, soft_end, splice.start, splice.end) {
        let overhang = (splice.end - read.record.pos() + 1) as usize;
        let read_start_index = (read.record.pos() - soft_start) as usize;
        let reference_start_index = splice.reference.len().saturating_sub(overhang);
        if overhanging_bases_mismatch(
            &read_bases,
            read_start_index,
            read_bases_len,
            &splice.reference,
            reference_start_index,
            overhang,
            options,
        ) {
            soft_clip_by_read_coordinates(&mut read.record, 0, (splice.end - soft_start) as usize)?;
        }
    } else if is_right_overhang(soft_start, soft_end, splice.start, splice.end) {
        let overhang = (soft_end - splice.start + 1) as usize;
        let span_to_test = (reference_end(&read.record) - splice.start + 1).max(0) as usize;
        if overhanging_bases_mismatch(
            &read_bases,
            read_bases.len().saturating_sub(overhang),
            read_bases_len,
            &splice.reference,
            0,
            span_to_test,
            options,
        ) {
            soft_clip_by_read_coordinates(
                &mut read.record,
                read_bases.len().saturating_sub(overhang),
                read_bases.len().saturating_sub(1),
            )?;
        }
    }
    Ok(())
}

fn is_left_overhang(read_start: i64, read_end: i64, splice_start: i64, splice_end: i64) -> bool {
    read_start <= splice_end && read_start > splice_start && read_end > splice_end
}

fn is_right_overhang(read_start: i64, read_end: i64, splice_start: i64, splice_end: i64) -> bool {
    read_end >= splice_start && read_end < splice_end && read_start < splice_start
}

fn overhanging_bases_mismatch(
    read: &[u8],
    read_start_index: usize,
    read_length: usize,
    reference: &[u8],
    reference_start_index: usize,
    span_to_test: usize,
    options: OverhangOptions,
) -> bool {
    if span_to_test < 1
        || span_to_test > options.max_bases_in_overhang
        || span_to_test > read_length / 2
    {
        return false;
    }
    if read_start_index + span_to_test > read.len()
        || reference_start_index + span_to_test > reference.len()
    {
        return false;
    }

    let mut mismatches = 0usize;
    for idx in 0..span_to_test {
        if !read[read_start_index + idx]
            .eq_ignore_ascii_case(&reference[reference_start_index + idx])
        {
            mismatches += 1;
            if mismatches > options.max_mismatches_in_overhang {
                return true;
            }
        }
    }
    mismatches >= span_to_test.div_ceil(2)
}

fn soft_reference_bounds(record: &bam::Record) -> Option<(i64, i64)> {
    if record.is_unmapped() {
        return None;
    }
    let cigar = record.cigar().iter().copied().collect::<Vec<_>>();
    let start = record.pos() - leading_soft_clip_len(&cigar) as i64;
    let span = cigar
        .iter()
        .copied()
        .map(reference_len_with_soft_clips)
        .sum::<u32>() as i64;
    if span <= 0 {
        None
    } else {
        Some((start, start + span - 1))
    }
}

fn reference_end(record: &bam::Record) -> i64 {
    let cigar = record.cigar().iter().copied().collect::<Vec<_>>();
    record.pos() + reference_span(&cigar) as i64 - 1
}

fn unclipped_read_bases_len(record: &bam::Record) -> usize {
    record
        .cigar()
        .iter()
        .copied()
        .filter(|element| {
            consumes_read_bases(*element)
                && !matches!(element, Cigar::SoftClip(_) | Cigar::HardClip(_))
        })
        .map(|element| element.len() as usize)
        .sum()
}

fn soft_clip_by_read_coordinates(
    record: &mut bam::Record,
    start: usize,
    stop: usize,
) -> Result<()> {
    if record.seq_len() == 0 || (start == 0 && stop + 1 >= record.seq_len()) {
        return Ok(());
    }
    if start > stop || stop >= record.seq_len() || (start > 0 && stop + 1 < record.seq_len()) {
        bail!(
            "invalid read-coordinate soft clip {start}-{stop} for CIGAR {}",
            record.cigar()
        );
    }

    let original: Vec<Cigar> = record.cigar().iter().copied().collect();
    let (mut clipped, ref_shift) = if start == 0 {
        clip_left_by_read_bases(&original, stop + 1)?
    } else {
        (
            clip_right_by_read_bases(&original, record.seq_len() - start)?,
            0,
        )
    };
    let (cleaned, edge_ref_shift) = clean_edge_deletions(clipped);
    clipped = cleaned;
    let new_pos = alignment_start_after_clipping(record.pos(), ref_shift + edge_ref_shift);
    record.set_cigar(Some(&CigarString(clipped)));
    record.set_pos(new_pos);
    update_record_bin(record);
    Ok(())
}

fn clip_left_by_read_bases(cigar: &[Cigar], mut clip_bases: usize) -> Result<(Vec<Cigar>, i64)> {
    let leading_hard = leading_hard_clip_len(cigar);
    let total_read = cigar.iter().copied().map(read_len).sum::<u32>();
    let trailing_clip_start = cigar
        .iter()
        .rposition(|element| !matches!(element, Cigar::HardClip(_) | Cigar::SoftClip(_)))
        .map(|idx| idx + 1)
        .unwrap_or(cigar.len());
    let mut kept = Vec::new();
    let mut soft = 0u32;
    let mut ref_shift = 0i64;
    for element in cigar[..trailing_clip_start].iter().copied() {
        if matches!(element, Cigar::HardClip(_)) {
            continue;
        }
        if matches!(element, Cigar::SoftClip(_)) {
            soft += element.len();
            clip_bases = clip_bases.saturating_sub(element.len() as usize);
            continue;
        }
        let len = read_len(element) as usize;
        if consumes_read_bases(element) && clip_bases > 0 {
            let clipped_here = clip_bases.min(len);
            soft += clipped_here as u32;
            clip_bases -= clipped_here;
            ref_shift += reference_len(with_len(element, clipped_here as u32)) as i64;
            let remaining = len - clipped_here;
            if remaining > 0 {
                push_cigar(&mut kept, with_len(element, remaining as u32));
            }
        } else if clip_bases > 0 && matches!(element, Cigar::Del(_)) {
            ref_shift += element.len() as i64;
        } else {
            push_cigar(&mut kept, element);
        }
    }
    if clip_bases > 0 {
        bail!("left soft clip exceeds read length");
    }
    let mut result = Vec::with_capacity(kept.len() + 3);
    push_cigar(&mut result, Cigar::HardClip(leading_hard));
    push_cigar(&mut result, Cigar::SoftClip(soft.min(total_read)));
    for element in kept {
        push_cigar(&mut result, element);
    }
    for element in cigar[trailing_clip_start..].iter().copied() {
        push_cigar(&mut result, element);
    }
    Ok((result, ref_shift))
}

fn clip_right_by_read_bases(cigar: &[Cigar], clip_bases: usize) -> Result<Vec<Cigar>> {
    let trailing_hard = trailing_hard_clip_len(cigar);
    let total_read = cigar.iter().copied().map(read_len).sum::<u32>();
    let mut clip_bases = clip_bases;
    let leading_clip_count = cigar
        .iter()
        .take_while(|element| matches!(element, Cigar::HardClip(_) | Cigar::SoftClip(_)))
        .count();
    let end = cigar.len().saturating_sub(usize::from(matches!(
        cigar.last(),
        Some(Cigar::HardClip(_))
    )));
    let mut kept_reversed = Vec::new();
    let mut soft = 0u32;

    for element in cigar[leading_clip_count..end].iter().copied().rev() {
        if matches!(element, Cigar::SoftClip(_)) {
            soft += element.len();
            clip_bases = clip_bases.saturating_sub(element.len() as usize);
            continue;
        }
        let len = read_len(element) as usize;
        if consumes_read_bases(element) && clip_bases > 0 {
            let clipped_here = clip_bases.min(len);
            soft += clipped_here as u32;
            clip_bases -= clipped_here;
            let remaining = len - clipped_here;
            if remaining > 0 {
                push_cigar(&mut kept_reversed, with_len(element, remaining as u32));
            }
        } else if clip_bases > 0 && matches!(element, Cigar::Del(_)) {
            continue;
        } else {
            push_cigar(&mut kept_reversed, element);
        }
    }
    if clip_bases > 0 {
        bail!("right soft clip exceeds read length");
    }

    let mut result = Vec::with_capacity(kept_reversed.len() + 3);
    for element in cigar[..leading_clip_count].iter().copied() {
        push_cigar(&mut result, element);
    }
    for element in kept_reversed.into_iter().rev() {
        push_cigar(&mut result, element);
    }
    push_cigar(&mut result, Cigar::SoftClip(soft.min(total_read)));
    push_cigar(&mut result, Cigar::HardClip(trailing_hard));
    Ok(result)
}

fn clean_edge_deletions(mut cigar: Vec<Cigar>) -> (Vec<Cigar>, i64) {
    let mut result = Vec::with_capacity(cigar.len());
    let mut seen_read_or_reference = false;
    let mut leading_ref_shift = 0i64;
    for element in cigar.drain(..) {
        if matches!(element, Cigar::Del(_)) && !seen_read_or_reference {
            leading_ref_shift += element.len() as i64;
            continue;
        }
        if !matches!(element, Cigar::SoftClip(_) | Cigar::HardClip(_)) {
            seen_read_or_reference = true;
        }
        push_cigar(&mut result, element);
    }
    while let Some(last) = result.last().copied() {
        if matches!(last, Cigar::Del(_)) {
            result.pop();
        } else {
            break;
        }
    }
    (result, leading_ref_shift)
}

fn alignment_start_after_clipping(original_pos: i64, ref_shift: i64) -> i64 {
    original_pos + ref_shift
}

fn update_record_bin(record: &mut bam::Record) {
    let cigar = record.cigar().iter().copied().collect::<Vec<_>>();
    record.set_bin(reg2bin(
        record.pos(),
        record.pos() + reference_span(&cigar) as i64,
    ));
}

fn push_cigar(cigar: &mut Vec<Cigar>, element: Cigar) {
    if element.len() == 0 {
        return;
    }

    if let Some(last) = cigar.last_mut() {
        if same_cigar_op(*last, element) {
            *last = with_len(*last, last.len() + element.len());
            return;
        }
    }

    cigar.push(element);
}

fn same_cigar_op(left: Cigar, right: Cigar) -> bool {
    matches!(
        (left, right),
        (Cigar::Match(_), Cigar::Match(_))
            | (Cigar::Ins(_), Cigar::Ins(_))
            | (Cigar::Del(_), Cigar::Del(_))
            | (Cigar::RefSkip(_), Cigar::RefSkip(_))
            | (Cigar::SoftClip(_), Cigar::SoftClip(_))
            | (Cigar::HardClip(_), Cigar::HardClip(_))
            | (Cigar::Pad(_), Cigar::Pad(_))
            | (Cigar::Equal(_), Cigar::Equal(_))
            | (Cigar::Diff(_), Cigar::Diff(_))
    )
}

fn with_len(element: Cigar, len: u32) -> Cigar {
    match element {
        Cigar::Match(_) => Cigar::Match(len),
        Cigar::Ins(_) => Cigar::Ins(len),
        Cigar::Del(_) => Cigar::Del(len),
        Cigar::RefSkip(_) => Cigar::RefSkip(len),
        Cigar::SoftClip(_) => Cigar::SoftClip(len),
        Cigar::HardClip(_) => Cigar::HardClip(len),
        Cigar::Pad(_) => Cigar::Pad(len),
        Cigar::Equal(_) => Cigar::Equal(len),
        Cigar::Diff(_) => Cigar::Diff(len),
    }
}

fn leading_hard_clip_len(cigar: &[Cigar]) -> u32 {
    cigar
        .first()
        .copied()
        .map(|element| match element {
            Cigar::HardClip(len) => len,
            _ => 0,
        })
        .unwrap_or(0)
}

fn leading_soft_clip_len(cigar: &[Cigar]) -> u32 {
    cigar
        .iter()
        .copied()
        .find_map(|element| match element {
            Cigar::SoftClip(len) => Some(len),
            Cigar::HardClip(_) => None,
            _ => Some(0),
        })
        .unwrap_or(0)
}

fn trailing_hard_clip_len(cigar: &[Cigar]) -> u32 {
    cigar
        .last()
        .copied()
        .map(|element| match element {
            Cigar::HardClip(len) => len,
            _ => 0,
        })
        .unwrap_or(0)
}

fn is_real_alignment_op(element: Cigar) -> bool {
    matches!(
        element,
        Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) | Cigar::Ins(_) | Cigar::Del(_)
    )
}

fn read_len(element: Cigar) -> u32 {
    match element {
        Cigar::Match(len)
        | Cigar::Ins(len)
        | Cigar::SoftClip(len)
        | Cigar::Equal(len)
        | Cigar::Diff(len) => len,
        Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => 0,
    }
}

fn consumes_read_bases(element: Cigar) -> bool {
    matches!(
        element,
        Cigar::Match(_) | Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::Equal(_) | Cigar::Diff(_)
    )
}

fn reference_len(element: Cigar) -> u32 {
    match element {
        Cigar::Match(len)
        | Cigar::Del(len)
        | Cigar::RefSkip(len)
        | Cigar::Equal(len)
        | Cigar::Diff(len) => len,
        Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => 0,
    }
}

fn reference_len_with_soft_clips(element: Cigar) -> u32 {
    match element {
        Cigar::SoftClip(len) => len,
        _ => reference_len(element),
    }
}

fn reference_span(cigar: &[Cigar]) -> u32 {
    cigar.iter().copied().map(reference_len).sum::<u32>().max(1)
}

fn reg2bin(beg: i64, end: i64) -> u16 {
    let beg = beg.max(0) as u32;
    let end = (end.max(1) - 1) as u32;

    if (beg >> 14) == (end >> 14) {
        return (4681 + (beg >> 14)) as u16;
    }
    if (beg >> 17) == (end >> 17) {
        return (585 + (beg >> 17)) as u16;
    }
    if (beg >> 20) == (end >> 20) {
        return (73 + (beg >> 20)) as u16;
    }
    if (beg >> 23) == (end >> 23) {
        return (9 + (beg >> 23)) as u16;
    }
    if (beg >> 26) == (end >> 26) {
        return (1 + (beg >> 26)) as u16;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Aux;

    fn cigar_text(cigar: &CigarString) -> String {
        cigar.to_string()
    }

    fn record_with_cigar(cigar: Vec<Cigar>, flags: u16) -> bam::Record {
        let cigar = CigarString(cigar);
        let read_len = cigar.iter().copied().map(read_len).sum::<u32>() as usize;
        let bases = vec![b'A'; read_len];
        let quals = vec![30u8; read_len];
        let mut record = bam::Record::new();
        record.set(b"read1", Some(&cigar), &bases, &quals);
        record.set_tid(0);
        record.set_pos(100);
        record.set_mtid(0);
        record.set_mpos(200);
        record.set_insert_size(150);
        record.set_mapq(255);
        record.set_flags(flags);
        record
    }

    #[test]
    fn splits_simple_n_cigar() {
        let plans =
            split_cigar_at_ref_skips(&[Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)])
                .unwrap();

        assert_eq!(plans.len(), 2);
        assert_eq!(plans[0].ref_offset, 0);
        assert_eq!(cigar_text(&plans[0].cigar), "5M5S");
        assert_eq!(plans[1].ref_offset, 15);
        assert_eq!(cigar_text(&plans[1].cigar), "5S5M");
    }

    #[test]
    fn preserves_soft_and_hard_clips() {
        let plans = split_cigar_at_ref_skips(&[
            Cigar::HardClip(1),
            Cigar::SoftClip(2),
            Cigar::Match(3),
            Cigar::RefSkip(8),
            Cigar::Match(4),
            Cigar::SoftClip(5),
            Cigar::HardClip(6),
        ])
        .unwrap();

        assert_eq!(plans.len(), 2);
        assert_eq!(cigar_text(&plans[0].cigar), "1H2S3M9S6H");
        assert_eq!(plans[1].ref_offset, 11);
        assert_eq!(cigar_text(&plans[1].cigar), "1H5S4M5S6H");
    }

    #[test]
    fn handles_insertions_and_deletions() {
        let plans = split_cigar_at_ref_skips(&[
            Cigar::HardClip(1),
            Cigar::Match(2),
            Cigar::Del(2),
            Cigar::Match(1),
            Cigar::RefSkip(2),
            Cigar::Match(1),
            Cigar::Ins(2),
            Cigar::RefSkip(1),
            Cigar::Match(1),
            Cigar::SoftClip(2),
        ])
        .unwrap();

        assert_eq!(plans.len(), 3);
        assert_eq!(cigar_text(&plans[0].cigar), "1H2M2D1M6S");
        assert_eq!(plans[1].ref_offset, 7);
        assert_eq!(cigar_text(&plans[1].cigar), "1H3S1M2I3S");
        assert_eq!(plans[2].ref_offset, 9);
        assert_eq!(cigar_text(&plans[2].cigar), "1H6S1M2S");
    }

    #[test]
    fn trims_deletions_at_split_edges() {
        let plans = split_cigar_at_ref_skips(&[
            Cigar::Match(4),
            Cigar::Del(3),
            Cigar::RefSkip(5),
            Cigar::Del(2),
            Cigar::Match(4),
        ])
        .unwrap();

        assert_eq!(plans.len(), 2);
        assert_eq!(cigar_text(&plans[0].cigar), "4M4S");
        assert_eq!(plans[1].ref_offset, 14);
        assert_eq!(cigar_text(&plans[1].cigar), "4S4M");
    }

    #[test]
    fn leaves_bogus_leading_n_only_split_unchanged() {
        let plans = split_cigar_at_ref_skips(&[
            Cigar::SoftClip(1),
            Cigar::RefSkip(3),
            Cigar::Match(2),
            Cigar::HardClip(4),
        ])
        .unwrap();

        assert!(plans.is_empty());
    }

    #[test]
    fn transform_splits_paired_read_and_removes_stale_tags() {
        let mut record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x1 | 0x2 | 0x40,
        );
        record.push_aux(b"NM", Aux::I32(1)).unwrap();
        record.push_aux(b"MD", Aux::String("5")).unwrap();
        record.push_aux(b"NH", Aux::I32(1)).unwrap();
        record.push_aux(b"MC", Aux::String("10M")).unwrap();

        let records = transform_record(&record, SplitOptions::default()).unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].mapq(), 60);
        assert_eq!(records[0].pos(), 100);
        assert_eq!(records[1].pos(), 115);
        assert_eq!(records[0].flags() & 0x1, 0x1);
        assert_eq!(records[1].flags() & SUPPLEMENTARY_FLAG, SUPPLEMENTARY_FLAG);
        assert_eq!(records[1].mtid(), 0);
        assert_eq!(records[1].mpos(), 200);
        assert!(records[0].aux(b"NM").is_err());
        assert!(records[0].aux(b"MD").is_err());
        assert!(records[0].aux(b"NH").is_err());
        assert!(records[0].aux(b"MC").is_err());
    }

    #[test]
    fn compatibility_mode_repairs_sa_tags_for_split_family() {
        let mut record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x1 | 0x2 | 0x40,
        );
        record
            .push_aux(b"SA", Aux::String("chr9,7,-,3M,40,2;"))
            .unwrap();
        let records = transform_record(
            &record,
            SplitOptions {
                mode: SplitMode::Compatibility,
                ..SplitOptions::default()
            },
        )
        .unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].flags() & SUPPLEMENTARY_FLAG, 0);
        assert_eq!(records[1].flags() & SUPPLEMENTARY_FLAG, SUPPLEMENTARY_FLAG);
        assert_eq!(
            records[0].aux(b"SA").unwrap(),
            Aux::String("chr9,7,-,3M,40,2;1,116,+,5S5M,60,*;")
        );
        assert_eq!(
            records[1].aux(b"SA").unwrap(),
            Aux::String("1,101,+,5M5S,60,*;chr9,7,-,3M,40,2;")
        );
    }

    #[test]
    fn compatibility_mode_uses_header_contig_names_in_sa_tags() {
        let record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0,
        );
        let contigs = vec!["chr1".to_string()];
        let records = transform_record_with_contig_names(
            &record,
            SplitOptions {
                mode: SplitMode::Compatibility,
                ..SplitOptions::default()
            },
            Some(&contigs),
        )
        .unwrap();

        assert_eq!(
            records[0].aux(b"SA").unwrap(),
            Aux::String("chr1,116,+,5S5M,60,*;")
        );
    }

    #[test]
    fn compatibility_mode_repairs_mate_cigar_to_first_split_segment() {
        let mut record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x1 | 0x2 | 0x40,
        );
        record.push_aux(b"MC", Aux::String("3M7N4M")).unwrap();

        let records = transform_record(
            &record,
            SplitOptions {
                mode: SplitMode::Compatibility,
                ..SplitOptions::default()
            },
        )
        .unwrap();

        assert_eq!(records[0].aux(b"MC").unwrap(), Aux::String("3M4S"));
        assert_eq!(records[1].aux(b"MC").unwrap(), Aux::String("3M4S"));
    }

    #[test]
    fn compatibility_mode_keeps_unsplit_record_repairs_mc_and_removes_stale_tags() {
        let mut record = record_with_cigar(vec![Cigar::Match(5)], 0x1 | 0x2 | 0x40);
        record.push_aux(b"NM", Aux::I32(1)).unwrap();
        record.push_aux(b"MD", Aux::String("5")).unwrap();
        record.push_aux(b"NH", Aux::I32(1)).unwrap();
        record.push_aux(b"MC", Aux::String("3M7N4M")).unwrap();

        let records = transform_record(
            &record,
            SplitOptions {
                mode: SplitMode::Compatibility,
                ..SplitOptions::default()
            },
        )
        .unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].aux(b"MC").unwrap(), Aux::String("3M4S"));
        assert!(records[0].aux(b"NM").is_err());
        assert!(records[0].aux(b"MD").is_err());
        assert!(records[0].aux(b"NH").is_err());
    }

    #[test]
    fn compatibility_mode_skips_secondary_but_repairs_mc() {
        let mut record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x100,
        );
        record.push_aux(b"MC", Aux::String("3M7N4M")).unwrap();

        let records = transform_record(
            &record,
            SplitOptions {
                mode: SplitMode::Compatibility,
                ..SplitOptions::default()
            },
        )
        .unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].flags() & 0x100, 0x100);
        assert_eq!(records[0].cigar().to_string(), "5M10N5M");
        assert_eq!(records[0].aux(b"MC").unwrap(), Aux::String("3M4S"));
    }

    #[test]
    fn overhang_fixing_soft_clips_left_mismatching_overhang() {
        let mut record = record_with_cigar(vec![Cigar::SoftClip(2), Cigar::Match(8)], 0);
        record.set_pos(103);
        let splice = Splice {
            tid: 0,
            start: 100,
            end: 104,
            reference: b"TTTTT".to_vec(),
        };
        let mut read = ManagedRead::new(record);

        fix_split_with_options(&mut read, &splice, OverhangOptions::default()).unwrap();

        assert_eq!(read.record.pos(), 105);
        assert_eq!(read.record.cigar().to_string(), "4S6M");
    }

    #[test]
    fn overhang_fixing_soft_clips_right_mismatching_overhang() {
        let mut record = record_with_cigar(vec![Cigar::Match(8), Cigar::SoftClip(2)], 0);
        record.set_pos(95);
        let splice = Splice {
            tid: 0,
            start: 100,
            end: 106,
            reference: b"TTTTTTT".to_vec(),
        };
        let mut read = ManagedRead::new(record);

        fix_split_with_options(&mut read, &splice, OverhangOptions::default()).unwrap();

        assert_eq!(read.record.pos(), 95);
        assert_eq!(read.record.cigar().to_string(), "5M5S");
    }

    #[test]
    fn overhang_fixing_preserves_opposite_end_soft_clip() {
        let mut record = record_with_cigar(
            vec![Cigar::SoftClip(1), Cigar::Match(72), Cigar::SoftClip(27)],
            0,
        );
        record.set_pos(1_717_961 - 1);

        soft_clip_by_read_coordinates(&mut record, 73, 99).unwrap();

        assert_eq!(record.pos(), 1_717_961 - 1);
        assert_eq!(record.cigar().to_string(), "1S72M27S");
    }

    #[test]
    fn overhang_fixing_left_clip_across_deletion_preserves_read_length() {
        let mut record =
            record_with_cigar(vec![Cigar::Match(30), Cigar::Del(1), Cigar::Match(55)], 0);
        record.set_pos(3_576_525 - 1);

        soft_clip_by_read_coordinates(&mut record, 0, 33).unwrap();

        assert_eq!(record.pos(), 3_576_560 - 1);
        assert_eq!(record.cigar().to_string(), "34S51M");
        assert_eq!(
            record.cigar().iter().copied().map(read_len).sum::<u32>() as usize,
            record.seq_len()
        );
    }

    #[test]
    fn overhang_fixing_left_clip_preserves_trailing_soft_clip_and_insertions() {
        let mut record = record_with_cigar(
            vec![
                Cigar::SoftClip(4),
                Cigar::Match(56),
                Cigar::Ins(2),
                Cigar::Match(21),
                Cigar::SoftClip(17),
            ],
            0,
        );
        record.set_pos(3_577_698 - 1);

        soft_clip_by_read_coordinates(&mut record, 0, 18).unwrap();

        assert_eq!(record.pos(), 3_577_713 - 1);
        assert_eq!(record.cigar().to_string(), "19S41M2I21M17S");
        assert_eq!(
            record.cigar().iter().copied().map(read_len).sum::<u32>() as usize,
            record.seq_len()
        );
    }

    #[test]
    fn overhang_fixing_keeps_matching_overhang() {
        let mut record = record_with_cigar(vec![Cigar::SoftClip(2), Cigar::Match(8)], 0);
        record.set_pos(103);
        let splice = Splice {
            tid: 0,
            start: 100,
            end: 104,
            reference: b"AAAAA".to_vec(),
        };
        let mut read = ManagedRead::new(record);

        fix_split_with_options(&mut read, &splice, OverhangOptions::default()).unwrap();

        assert_eq!(read.record.pos(), 103);
        assert_eq!(read.record.cigar().to_string(), "2S8M");
    }

    #[test]
    fn transform_skips_secondary_by_default() {
        let record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x100,
        );

        let records = transform_record(&record, SplitOptions::default()).unwrap();

        assert_eq!(records.len(), 1);
        assert_eq!(records[0].mapq(), 60);
        assert_eq!(records[0].cigar().to_string(), "5M10N5M");
    }

    #[test]
    fn transform_can_process_secondary_alignments() {
        let record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            0x100,
        );
        let records = transform_record(
            &record,
            SplitOptions {
                process_secondary_alignments: true,
                ..SplitOptions::default()
            },
        )
        .unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].flags() & 0x100, 0x100);
        assert_eq!(records[1].flags() & 0x100, 0x100);
        assert_eq!(records[1].flags() & SUPPLEMENTARY_FLAG, SUPPLEMENTARY_FLAG);
    }

    #[test]
    fn transform_preserves_supplementary_and_unmapped_mate_flags() {
        let record = record_with_cigar(
            vec![Cigar::Match(5), Cigar::RefSkip(10), Cigar::Match(5)],
            SUPPLEMENTARY_FLAG | 0x1 | 0x8,
        );

        let records = transform_record(&record, SplitOptions::default()).unwrap();

        assert_eq!(records.len(), 2);
        assert_eq!(records[0].flags() & SUPPLEMENTARY_FLAG, SUPPLEMENTARY_FLAG);
        assert_eq!(records[1].flags() & SUPPLEMENTARY_FLAG, SUPPLEMENTARY_FLAG);
        assert_eq!(records[0].flags() & 0x8, 0x8);
        assert_eq!(records[1].flags() & 0x8, 0x8);
    }
}

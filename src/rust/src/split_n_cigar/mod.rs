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

include!("transform.rs");
include!("split_logic.rs");
include!("tests_block.rs");

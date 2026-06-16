use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub struct CliError {
    pub message: String,
    pub exit_code: i32,
}

impl CliError {
    fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            exit_code: 2,
        }
    }
}

type Result<T> = std::result::Result<T, CliError>;

#[derive(Clone, Debug)]
struct DictRecord {
    length: u64,
}

#[derive(Clone, Debug)]
struct SequenceDict {
    header_lines: Vec<String>,
    records: Vec<DictRecord>,
    index_by_name: HashMap<String, usize>,
}

impl SequenceDict {
    fn order(&self, contig: &str) -> Option<usize> {
        self.index_by_name.get(contig).copied()
    }

    fn contig_length(&self, contig: &str) -> Option<u64> {
        self.order(contig).map(|idx| self.records[idx].length)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Interval {
    contig: String,
    start: u64,
    end: u64,
}

impl Interval {
    fn len(&self) -> u64 {
        self.end - self.start + 1
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum MergeRule {
    All,
    OverlapOnly,
}

pub fn run(args: Vec<String>) -> Result<()> {
    if args.len() <= 1 || is_help(&args[1]) {
        print_top_help();
        return Ok(());
    }

    match args[1].as_str() {
        "bed-to-interval-list" => run_bed_to_interval_list(&args[2..]),
        "split-intervals" => run_split_intervals(&args[2..]),
        "prepare" => run_prepare(&args[2..]),
        other => Err(CliError::new(format!(
            "unknown subcommand '{other}'. Run rust_interval_tools --help"
        ))),
    }
}

fn run_bed_to_interval_list(args: &[String]) -> Result<()> {
    if args.iter().any(|arg| is_help(arg)) {
        print_bed_help();
        return Ok(());
    }
    let opts = parse_options(
        args,
        &[
            "input-bed",
            "ref-dict",
            "output-interval-list",
            "output-merged-bed",
            "merge-rule",
        ],
    )?;
    let input_bed = required_path(&opts, "input-bed")?;
    let ref_dict = required_path(&opts, "ref-dict")?;
    let output_interval_list = required_path(&opts, "output-interval-list")?;
    let output_merged_bed = optional_path(&opts, "output-merged-bed");
    let merge_rule = parse_merge_rule(opts.get("merge-rule"))?;

    let dict = read_sequence_dict(&ref_dict)?;
    let mut intervals = read_bed_intervals(&input_bed, &dict)?;
    sort_and_merge(&mut intervals, &dict, merge_rule)?;
    write_interval_list(&output_interval_list, &dict, &intervals)?;
    if let Some(path) = output_merged_bed {
        write_bed(&path, &intervals)?;
    }
    Ok(())
}

fn run_split_intervals(args: &[String]) -> Result<()> {
    if args.iter().any(|arg| is_help(arg)) {
        print_split_help();
        return Ok(());
    }
    let opts = parse_options(
        args,
        &[
            "input-interval-list",
            "output-dir",
            "scatter-count",
            "merge-rule",
        ],
    )?;
    let input_interval_list = required_path(&opts, "input-interval-list")?;
    let output_dir = required_path(&opts, "output-dir")?;
    let scatter_count = required_usize(&opts, "scatter-count")?;
    let merge_rule = parse_merge_rule(opts.get("merge-rule"))?;

    let (dict, mut intervals) = read_interval_list(&input_interval_list)?;
    sort_and_merge(&mut intervals, &dict, merge_rule)?;
    write_scattered_interval_lists(&output_dir, &dict, &intervals, scatter_count)?;
    Ok(())
}

fn run_prepare(args: &[String]) -> Result<()> {
    if args.iter().any(|arg| is_help(arg)) {
        print_prepare_help();
        return Ok(());
    }
    let opts = parse_options(
        args,
        &[
            "input-bed",
            "ref-dict",
            "output-merged-bed",
            "output-interval-list",
            "output-dir",
            "scatter-count",
            "merge-rule",
        ],
    )?;
    let input_bed = required_path(&opts, "input-bed")?;
    let ref_dict = required_path(&opts, "ref-dict")?;
    let output_merged_bed = required_path(&opts, "output-merged-bed")?;
    let output_interval_list = required_path(&opts, "output-interval-list")?;
    let output_dir = required_path(&opts, "output-dir")?;
    let scatter_count = required_usize(&opts, "scatter-count")?;
    let merge_rule = parse_merge_rule(opts.get("merge-rule"))?;

    let dict = read_sequence_dict(&ref_dict)?;
    let mut intervals = read_bed_intervals(&input_bed, &dict)?;
    sort_and_merge(&mut intervals, &dict, merge_rule)?;
    write_bed(&output_merged_bed, &intervals)?;
    write_interval_list(&output_interval_list, &dict, &intervals)?;
    write_scattered_interval_lists(&output_dir, &dict, &intervals, scatter_count)?;
    Ok(())
}

fn parse_options(args: &[String], allowed: &[&str]) -> Result<HashMap<String, String>> {
    let mut opts = HashMap::new();
    let mut i = 0;
    while i < args.len() {
        let key = args[i]
            .strip_prefix("--")
            .ok_or_else(|| CliError::new(format!("expected option, found '{}'", args[i])))?;
        if !allowed.contains(&key) {
            return Err(CliError::new(format!("unknown option '--{key}'")));
        }
        let value = args
            .get(i + 1)
            .ok_or_else(|| CliError::new(format!("missing value for '--{key}'")))?;
        if value.starts_with("--") {
            return Err(CliError::new(format!("missing value for '--{key}'")));
        }
        if opts.insert(key.to_string(), value.clone()).is_some() {
            return Err(CliError::new(format!("duplicate option '--{key}'")));
        }
        i += 2;
    }
    Ok(opts)
}

fn required_path(opts: &HashMap<String, String>, key: &str) -> Result<PathBuf> {
    opts.get(key)
        .map(PathBuf::from)
        .ok_or_else(|| CliError::new(format!("missing required option '--{key}'")))
}

fn optional_path(opts: &HashMap<String, String>, key: &str) -> Option<PathBuf> {
    opts.get(key).map(PathBuf::from)
}

fn required_usize(opts: &HashMap<String, String>, key: &str) -> Result<usize> {
    let value = opts
        .get(key)
        .ok_or_else(|| CliError::new(format!("missing required option '--{key}'")))?;
    let parsed = value
        .parse::<usize>()
        .map_err(|_| CliError::new(format!("invalid integer for '--{key}': {value}")))?;
    if parsed == 0 {
        return Err(CliError::new(format!("--{key} must be at least 1")));
    }
    Ok(parsed)
}

fn parse_merge_rule(value: Option<&String>) -> Result<MergeRule> {
    match value.map(String::as_str).unwrap_or("all") {
        "all" => Ok(MergeRule::All),
        "overlap-only" => Ok(MergeRule::OverlapOnly),
        other => Err(CliError::new(format!(
            "invalid --merge-rule '{other}', expected 'all' or 'overlap-only'"
        ))),
    }
}

fn read_sequence_dict(path: &Path) -> Result<SequenceDict> {
    let file = File::open(path)
        .map_err(|err| CliError::new(format!("failed to open {}: {err}", path.display())))?;
    let reader = BufReader::new(file);
    parse_dict_lines(reader.lines(), path)
}

fn parse_dict_lines<I>(lines: I, path: &Path) -> Result<SequenceDict>
where
    I: IntoIterator<Item = std::io::Result<String>>,
{
    let mut header_lines = Vec::new();
    let mut records = Vec::new();
    let mut index_by_name = HashMap::new();

    for line in lines {
        let line =
            line.map_err(|err| CliError::new(format!("failed to read {}: {err}", path.display())))?;
        if !line.starts_with('@') {
            continue;
        }
        if line.starts_with("@SQ\t") {
            let mut name = None;
            let mut length = None;
            for field in line.split('\t').skip(1) {
                if let Some(value) = field.strip_prefix("SN:") {
                    name = Some(value.to_string());
                } else if let Some(value) = field.strip_prefix("LN:") {
                    length = Some(value.parse::<u64>().map_err(|_| {
                        CliError::new(format!("invalid LN value in {}: {value}", path.display()))
                    })?);
                }
            }
            let name = name.ok_or_else(|| {
                CliError::new(format!(
                    "missing SN field in @SQ line in {}",
                    path.display()
                ))
            })?;
            let length = length.ok_or_else(|| {
                CliError::new(format!(
                    "missing LN field in @SQ line in {}",
                    path.display()
                ))
            })?;
            if length == 0 {
                return Err(CliError::new(format!(
                    "zero-length contig '{name}' in {}",
                    path.display()
                )));
            }
            if index_by_name.contains_key(&name) {
                return Err(CliError::new(format!(
                    "duplicate contig '{name}' in {}",
                    path.display()
                )));
            }
            index_by_name.insert(name.clone(), records.len());
            records.push(DictRecord { length });
        }
        header_lines.push(line);
    }

    if records.is_empty() {
        return Err(CliError::new(format!(
            "no @SQ records found in {}",
            path.display()
        )));
    }

    Ok(SequenceDict {
        header_lines,
        records,
        index_by_name,
    })
}

fn read_bed_intervals(path: &Path, dict: &SequenceDict) -> Result<Vec<Interval>> {
    let file = File::open(path)
        .map_err(|err| CliError::new(format!("failed to open {}: {err}", path.display())))?;
    let reader = BufReader::new(file);
    let mut intervals = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line =
            line.map_err(|err| CliError::new(format!("failed to read {}: {err}", path.display())))?;
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track ")
            || trimmed.starts_with("browser ")
        {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 3 {
            return Err(CliError::new(format!(
                "{}:{}: BED row has fewer than 3 columns",
                path.display(),
                line_no + 1
            )));
        }
        let contig = fields[0].to_string();
        let start0 = parse_u64(fields[1], path, line_no + 1, "BED start")?;
        let end0 = parse_u64(fields[2], path, line_no + 1, "BED end")?;
        if end0 <= start0 {
            return Err(CliError::new(format!(
                "{}:{}: BED end must be greater than start",
                path.display(),
                line_no + 1
            )));
        }
        let length = dict.contig_length(&contig).ok_or_else(|| {
            CliError::new(format!(
                "{}:{}: contig '{contig}' is not present in the sequence dictionary",
                path.display(),
                line_no + 1
            ))
        })?;
        if end0 > length {
            return Err(CliError::new(format!(
                "{}:{}: BED interval {contig}:{start0}-{end0} extends past contig length {length}",
                path.display(),
                line_no + 1
            )));
        }
        intervals.push(Interval {
            contig,
            start: start0 + 1,
            end: end0,
        });
    }
    if intervals.is_empty() {
        return Err(CliError::new(format!(
            "no intervals found in {}",
            path.display()
        )));
    }
    Ok(intervals)
}

fn read_interval_list(path: &Path) -> Result<(SequenceDict, Vec<Interval>)> {
    let file = File::open(path)
        .map_err(|err| CliError::new(format!("failed to open {}: {err}", path.display())))?;
    let reader = BufReader::new(file);
    let mut header_lines = Vec::new();
    let mut body_lines = Vec::new();
    for line in reader.lines() {
        let line =
            line.map_err(|err| CliError::new(format!("failed to read {}: {err}", path.display())))?;
        if line.starts_with('@') {
            header_lines.push(Ok(line));
        } else {
            body_lines.push(line);
        }
    }
    let dict = parse_dict_lines(header_lines, path)?;
    let mut intervals = Vec::new();
    for (line_idx, line) in body_lines.iter().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() < 3 {
            return Err(CliError::new(format!(
                "{}: interval row has fewer than 3 columns near body line {}",
                path.display(),
                line_idx + 1
            )));
        }
        let contig = fields[0].to_string();
        let start = parse_u64(fields[1], path, line_idx + 1, "interval start")?;
        let end = parse_u64(fields[2], path, line_idx + 1, "interval end")?;
        validate_interval(path, line_idx + 1, &contig, start, end, &dict)?;
        intervals.push(Interval { contig, start, end });
    }
    if intervals.is_empty() {
        return Err(CliError::new(format!(
            "no intervals found in {}",
            path.display()
        )));
    }
    Ok((dict, intervals))
}

fn parse_u64(value: &str, path: &Path, line_no: usize, label: &str) -> Result<u64> {
    value.parse::<u64>().map_err(|_| {
        CliError::new(format!(
            "{}:{line_no}: invalid {label} value '{value}'",
            path.display()
        ))
    })
}

fn validate_interval(
    path: &Path,
    line_no: usize,
    contig: &str,
    start: u64,
    end: u64,
    dict: &SequenceDict,
) -> Result<()> {
    if start == 0 {
        return Err(CliError::new(format!(
            "{}:{line_no}: interval start must be at least 1",
            path.display()
        )));
    }
    if start > end {
        return Err(CliError::new(format!(
            "{}:{line_no}: interval start is greater than end",
            path.display()
        )));
    }
    let length = dict.contig_length(contig).ok_or_else(|| {
        CliError::new(format!(
            "{}:{line_no}: contig '{contig}' is not present in the sequence dictionary",
            path.display()
        ))
    })?;
    if end > length {
        return Err(CliError::new(format!(
            "{}:{line_no}: interval {contig}:{start}-{end} extends past contig length {length}",
            path.display()
        )));
    }
    Ok(())
}

fn sort_and_merge(
    intervals: &mut Vec<Interval>,
    dict: &SequenceDict,
    merge_rule: MergeRule,
) -> Result<()> {
    for interval in intervals.iter() {
        if dict.order(&interval.contig).is_none() {
            return Err(CliError::new(format!(
                "contig '{}' is not present in the sequence dictionary",
                interval.contig
            )));
        }
    }
    intervals.sort_by(|a, b| {
        dict.order(&a.contig)
            .cmp(&dict.order(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    let mut merged: Vec<Interval> = Vec::with_capacity(intervals.len());
    for interval in intervals.drain(..) {
        if let Some(current) = merged.last_mut() {
            let merge_distance = match merge_rule {
                MergeRule::All => 1,
                MergeRule::OverlapOnly => 0,
            };
            if current.contig == interval.contig
                && interval.start <= current.end.saturating_add(merge_distance)
            {
                current.end = current.end.max(interval.end);
                continue;
            }
        }
        merged.push(interval);
    }
    *intervals = merged;
    Ok(())
}

fn write_interval_list(path: &Path, dict: &SequenceDict, intervals: &[Interval]) -> Result<()> {
    create_parent_dir(path)?;
    let file = File::create(path)
        .map_err(|err| CliError::new(format!("failed to create {}: {err}", path.display())))?;
    let mut writer = BufWriter::new(file);
    for line in &dict.header_lines {
        writeln!(writer, "{line}").map_err(write_error(path))?;
    }
    for interval in intervals {
        writeln!(
            writer,
            "{}\t{}\t{}\t+\t.",
            interval.contig, interval.start, interval.end
        )
        .map_err(write_error(path))?;
    }
    Ok(())
}

fn write_bed(path: &Path, intervals: &[Interval]) -> Result<()> {
    create_parent_dir(path)?;
    let file = File::create(path)
        .map_err(|err| CliError::new(format!("failed to create {}: {err}", path.display())))?;
    let mut writer = BufWriter::new(file);
    for interval in intervals {
        writeln!(
            writer,
            "{}\t{}\t{}",
            interval.contig,
            interval.start - 1,
            interval.end
        )
        .map_err(write_error(path))?;
    }
    Ok(())
}

fn write_scattered_interval_lists(
    output_dir: &Path,
    dict: &SequenceDict,
    intervals: &[Interval],
    scatter_count: usize,
) -> Result<()> {
    if intervals.is_empty() {
        return Err(CliError::new("cannot scatter an empty interval list"));
    }
    fs::create_dir_all(output_dir).map_err(|err| {
        CliError::new(format!(
            "failed to create output directory {}: {err}",
            output_dir.display()
        ))
    })?;
    remove_existing_scatter_outputs(output_dir)?;
    let shards = split_balanced(intervals, scatter_count)?;
    let width = usize::max(4, digits(shards.len().saturating_sub(1)));
    for (idx, shard) in shards.iter().enumerate() {
        let path = output_dir.join(format!("{idx:0width$}-scattered.interval_list"));
        write_interval_list(&path, dict, shard)?;
    }
    Ok(())
}

fn remove_existing_scatter_outputs(output_dir: &Path) -> Result<()> {
    for entry in fs::read_dir(output_dir).map_err(|err| {
        CliError::new(format!(
            "failed to read output directory {}: {err}",
            output_dir.display()
        ))
    })? {
        let entry = entry.map_err(|err| {
            CliError::new(format!(
                "failed to read output directory {}: {err}",
                output_dir.display()
            ))
        })?;
        let path = entry.path();
        if path
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.ends_with("-scattered.interval_list"))
        {
            fs::remove_file(&path).map_err(|err| {
                CliError::new(format!("failed to remove {}: {err}", path.display()))
            })?;
        }
    }
    Ok(())
}

fn split_balanced(intervals: &[Interval], scatter_count: usize) -> Result<Vec<Vec<Interval>>> {
    if scatter_count == 0 {
        return Err(CliError::new("scatter count must be at least 1"));
    }
    if intervals.is_empty() {
        return Err(CliError::new("cannot scatter an empty interval list"));
    }
    let shard_count = usize::min(scatter_count, intervals.len());
    if shard_count == 1 {
        return Ok(vec![intervals.to_vec()]);
    }

    let mut prefix = Vec::with_capacity(intervals.len() + 1);
    prefix.push(0_u128);
    for interval in intervals {
        let next = prefix[prefix.len() - 1] + u128::from(interval.len());
        prefix.push(next);
    }
    let total = prefix[prefix.len() - 1];
    let mut cuts = Vec::with_capacity(shard_count - 1);
    let mut previous_cut = 0_usize;
    for shard_index in 1..shard_count {
        let min_cut = previous_cut + 1;
        let max_cut = intervals.len() - (shard_count - shard_index);
        let desired =
            (total * shard_index as u128 + (shard_count as u128 / 2)) / shard_count as u128;
        let lower = lower_bound(&prefix, desired);
        let mut best = min_cut;
        let mut best_distance = u128::MAX;
        let start = lower.saturating_sub(8).max(min_cut);
        let end = usize::min(lower + 8, max_cut);
        for candidate in start..=end {
            let distance = abs_diff(prefix[candidate], desired);
            if distance < best_distance {
                best = candidate;
                best_distance = distance;
            }
        }
        cuts.push(best);
        previous_cut = best;
    }

    let mut shards = Vec::with_capacity(shard_count);
    let mut start = 0_usize;
    for cut in cuts.into_iter().chain(std::iter::once(intervals.len())) {
        shards.push(intervals[start..cut].to_vec());
        start = cut;
    }
    Ok(shards)
}

fn lower_bound(values: &[u128], target: u128) -> usize {
    let mut left = 0_usize;
    let mut right = values.len();
    while left < right {
        let mid = left + (right - left) / 2;
        if values[mid] < target {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    left
}

fn abs_diff(left: u128, right: u128) -> u128 {
    left.max(right) - left.min(right)
}

fn digits(value: usize) -> usize {
    if value == 0 {
        return 1;
    }
    let mut digits = 0;
    let mut value = value;
    while value > 0 {
        digits += 1;
        value /= 10;
    }
    digits
}

fn create_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).map_err(|err| {
                CliError::new(format!(
                    "failed to create directory {}: {err}",
                    parent.display()
                ))
            })?;
        }
    }
    Ok(())
}

fn write_error(path: &Path) -> impl FnOnce(std::io::Error) -> CliError + '_ {
    move |err| CliError::new(format!("failed to write {}: {err}", path.display()))
}

fn is_help(arg: &str) -> bool {
    arg == "-h" || arg == "--help"
}

fn print_top_help() {
    println!(
        "rust_interval_tools\n\
\n\
Usage:\n\
  rust_interval_tools <subcommand> [options]\n\
\n\
Subcommands:\n\
  bed-to-interval-list  Convert BED to sorted, merged GATK interval_list\n\
  split-intervals       Split an interval_list into balanced GATK shards\n\
  prepare               Convert BED, write merged BED, interval_list, and shards\n\
\n\
Run rust_interval_tools <subcommand> --help for details."
    );
}

fn print_bed_help() {
    println!(
        "Usage:\n\
  rust_interval_tools bed-to-interval-list \\\n\
    --input-bed PATH --ref-dict PATH --output-interval-list PATH \\\n\
    [--output-merged-bed PATH] [--merge-rule all|overlap-only]\n\
\n\
BED coordinates are interpreted as 0-based half-open and written as 1-based inclusive intervals."
    );
}

fn print_split_help() {
    println!(
        "Usage:\n\
  rust_interval_tools split-intervals \\\n\
    --input-interval-list PATH --scatter-count N --output-dir DIR \\\n\
    [--merge-rule all|overlap-only]\n\
\n\
Shard filenames match GATK SplitIntervals: 0000-scattered.interval_list."
    );
}

fn print_prepare_help() {
    println!(
        "Usage:\n\
  rust_interval_tools prepare \\\n\
    --input-bed PATH --ref-dict PATH --output-merged-bed PATH \\\n\
    --output-interval-list PATH --scatter-count N --output-dir DIR \\\n\
    [--merge-rule all|overlap-only]"
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn test_dict() -> SequenceDict {
        let lines = vec![
            Ok("@HD\tVN:1.6\tSO:coordinate".to_string()),
            Ok("@SQ\tSN:chr2\tLN:100".to_string()),
            Ok("@SQ\tSN:chr1\tLN:100".to_string()),
        ];
        parse_dict_lines(lines, Path::new("test.dict")).unwrap()
    }

    #[test]
    fn bed_rows_convert_sort_and_merge_by_dictionary_order() {
        let dict = test_dict();
        let mut intervals = vec![
            Interval {
                contig: "chr1".to_string(),
                start: 1,
                end: 10,
            },
            Interval {
                contig: "chr2".to_string(),
                start: 5,
                end: 8,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 11,
                end: 20,
            },
        ];
        sort_and_merge(&mut intervals, &dict, MergeRule::All).unwrap();
        assert_eq!(
            intervals,
            vec![
                Interval {
                    contig: "chr2".to_string(),
                    start: 5,
                    end: 8,
                },
                Interval {
                    contig: "chr1".to_string(),
                    start: 1,
                    end: 20,
                },
            ]
        );
    }

    #[test]
    fn overlap_only_merge_keeps_bookended_intervals_separate() {
        let dict = test_dict();
        let mut intervals = vec![
            Interval {
                contig: "chr1".to_string(),
                start: 1,
                end: 10,
            },
            Interval {
                contig: "chr1".to_string(),
                start: 11,
                end: 20,
            },
        ];
        sort_and_merge(&mut intervals, &dict, MergeRule::OverlapOnly).unwrap();
        assert_eq!(intervals.len(), 2);
    }

    #[test]
    fn invalid_contig_fails_validation() {
        let dict = test_dict();
        let err =
            validate_interval(Path::new("x.interval_list"), 1, "chr3", 1, 10, &dict).unwrap_err();
        assert!(err.message.contains("not present"));
    }

    #[test]
    fn balanced_split_preserves_total_coverage_without_splitting_intervals() {
        let intervals: Vec<Interval> = (0..10)
            .map(|idx| Interval {
                contig: "chr1".to_string(),
                start: idx * 10 + 1,
                end: idx * 10 + 10,
            })
            .collect();
        let shards = split_balanced(&intervals, 3).unwrap();
        let lengths: Vec<u64> = shards
            .iter()
            .map(|shard| shard.iter().map(Interval::len).sum())
            .collect();
        assert_eq!(lengths.iter().sum::<u64>(), 100);
        assert_eq!(lengths, vec![30, 40, 30]);
    }

    #[test]
    fn write_interval_list_preserves_header_and_coordinates() {
        let dict = test_dict();
        let dir = unique_temp_dir();
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("out.interval_list");
        let intervals = vec![Interval {
            contig: "chr2".to_string(),
            start: 5,
            end: 8,
        }];
        write_interval_list(&path, &dict, &intervals).unwrap();
        let text = fs::read_to_string(path).unwrap();
        assert!(text.contains("@SQ\tSN:chr2\tLN:100\n"));
        assert!(text.ends_with("chr2\t5\t8\t+\t.\n"));
        fs::remove_dir_all(dir).unwrap();
    }

    fn unique_temp_dir() -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("rust_interval_tools_test_{nanos}"))
    }
}

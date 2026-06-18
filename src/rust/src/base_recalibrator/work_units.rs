fn build_work_units(header: &bam::HeaderView, region_bases: u64) -> Result<Vec<WorkUnit>> {
    let mut work_units = Vec::new();
    for tid in 0..header.target_count() {
        let contig = String::from_utf8(header.tid2name(tid).to_vec())
            .with_context(|| format!("BAM header target id {tid} is not valid UTF-8"))?;
        let Some(contig_len) = header.target_len(tid) else {
            continue;
        };
        let mut start0 = 0_u64;
        while start0 < contig_len {
            let end0_exclusive = start0.saturating_add(region_bases).min(contig_len);
            work_units.push(WorkUnit {
                tid,
                contig: contig.clone(),
                start0,
                end0_exclusive,
            });
            start0 = end0_exclusive;
        }
    }
    Ok(work_units)
}

#[allow(clippy::too_many_arguments)]
fn process_work_unit(
    args: &Args,
    worker_index: usize,
    work: &WorkUnit,
    bam_reader: &mut bam::IndexedReader,
    header: &bam::HeaderView,
    rg_info: &ReadGroupInfo,
    reference: &mut ReferenceCache,
    known_sites: &mut KnownSiteCache,
    tables: &mut RecalTables,
) -> Result<()> {
    let start = i64::try_from(work.start0)
        .with_context(|| format!("{} start coordinate exceeds i64", work.contig))?;
    let end = i64::try_from(work.end0_exclusive)
        .with_context(|| format!("{} end coordinate exceeds i64", work.contig))?;
    bam_reader
        .fetch((work.tid as i32, start, end))
        .with_context(|| {
            format!(
                "worker {worker_index} failed to fetch BAM region {}:{}-{}",
                work.contig,
                work.start0 + 1,
                work.end0_exclusive
            )
        })?;

    for result in bam_reader.records() {
        let record = result.with_context(|| {
            format!(
                "worker {worker_index} failed to read BAM region {}:{}-{}",
                work.contig,
                work.start0 + 1,
                work.end0_exclusive
            )
        })?;
        if record.tid() != work.tid as i32 || record.pos() < start || record.pos() >= end {
            continue;
        }
        process_record(
            record,
            header,
            rg_info,
            reference,
            known_sites,
            args,
            tables,
        )?;
    }
    Ok(())
}

fn process_record(
    record: bam::Record,
    header: &bam::HeaderView,
    rg_info: &ReadGroupInfo,
    reference: &mut ReferenceCache,
    known_sites: &mut KnownSiteCache,
    args: &Args,
    tables: &mut RecalTables,
) -> Result<()> {
    if !passes_bqsr_filters(&record) {
        return Ok(());
    }
    let Some(read) = prepare_read(&record, header, rg_info, args)? else {
        return Ok(());
    };
    let reference_bases = reference.bases_for(&read.contig, read.start0, read.end0)?;
    let known = known_sites.intervals_for(&read.contig, read.start0, read.end0)?;
    process_read(&read, reference_bases, known, args, tables)
}

fn passes_bqsr_filters(record: &bam::Record) -> bool {
    !record.is_unmapped()
        && !record.is_secondary()
        && !record.is_duplicate()
        && !record.is_quality_check_failed()
        && record.mapq() != 0
        && record.mapq() != 255
        && record.tid() >= 0
        && record.pos() >= 0
}

fn read_groups_from_header(header: &bam::HeaderView) -> Result<ReadGroupInfo> {
    let text = String::from_utf8_lossy(header.as_bytes());
    let mut identifiers = Vec::new();
    let mut id_to_identifier = HashMap::new();
    let mut identifier_to_index = HashMap::new();
    for line in text.lines() {
        if !line.starts_with("@RG\t") {
            continue;
        }
        let mut id = None;
        let mut platform_unit = None;
        for field in line.split('\t').skip(1) {
            if let Some(value) = field.strip_prefix("ID:") {
                id = Some(value.to_string());
            } else if let Some(value) = field.strip_prefix("PU:") {
                platform_unit = Some(value.to_string());
            }
        }
        let id = id.context("encountered @RG header line without ID")?;
        let identifier = platform_unit.unwrap_or_else(|| id.clone());
        if !identifier_to_index.contains_key(&identifier) {
            identifier_to_index.insert(identifier.clone(), identifiers.len());
            identifiers.push(identifier.clone());
        }
        id_to_identifier.insert(id, identifier);
    }

    if identifiers.is_empty() {
        bail!("input BAM header has no read groups");
    }
    Ok(ReadGroupInfo {
        identifiers,
        id_to_identifier,
    })
}

fn prepare_read(
    record: &bam::Record,
    header: &bam::HeaderView,
    rg_info: &ReadGroupInfo,
    args: &Args,
) -> Result<Option<PreparedRead>> {
    let contig = String::from_utf8_lossy(header.tid2name(record.tid() as u32)).to_string();
    let rg_id = aux_string(record, b"RG")?;
    let rg_identifier = rg_info
        .id_to_identifier
        .get(&rg_id)
        .with_context(|| format!("read references unknown read group '{rg_id}'"))?;
    let rg_index = rg_info
        .identifiers
        .iter()
        .position(|identifier| identifier == rg_identifier)
        .with_context(|| format!("read group '{rg_identifier}' missing from lookup"))?;

    let bases = record.seq().as_bytes();
    let quals = if args.use_original_qualities {
        original_qualities(record)?.unwrap_or_else(|| record.qual().to_vec())
    } else {
        record.qual().to_vec()
    };
    if quals.len() != bases.len() {
        bail!(
            "read {} has {} bases but {} qualities",
            String::from_utf8_lossy(record.qname()),
            bases.len(),
            quals.len()
        );
    }

    let mut clipped_bases = Vec::with_capacity(bases.len());
    let mut clipped_quals = Vec::with_capacity(quals.len());
    let mut cigar = Vec::new();
    let mut read_pos = 0_usize;
    let mut ref_len = 0_u64;

    for op in &record.cigar() {
        match *op {
            Cigar::Match(len) => {
                copy_read_segment(
                    &bases,
                    &quals,
                    read_pos,
                    len,
                    &mut clipped_bases,
                    &mut clipped_quals,
                );
                read_pos += len as usize;
                ref_len += u64::from(len);
                cigar.push(SimpleCigar::Match(len));
            }
            Cigar::Equal(len) => {
                copy_read_segment(
                    &bases,
                    &quals,
                    read_pos,
                    len,
                    &mut clipped_bases,
                    &mut clipped_quals,
                );
                read_pos += len as usize;
                ref_len += u64::from(len);
                cigar.push(SimpleCigar::Equal(len));
            }
            Cigar::Diff(len) => {
                copy_read_segment(
                    &bases,
                    &quals,
                    read_pos,
                    len,
                    &mut clipped_bases,
                    &mut clipped_quals,
                );
                read_pos += len as usize;
                ref_len += u64::from(len);
                cigar.push(SimpleCigar::Diff(len));
            }
            Cigar::Ins(len) => {
                copy_read_segment(
                    &bases,
                    &quals,
                    read_pos,
                    len,
                    &mut clipped_bases,
                    &mut clipped_quals,
                );
                read_pos += len as usize;
                cigar.push(SimpleCigar::Ins(len));
            }
            Cigar::Del(len) => {
                ref_len += u64::from(len);
                cigar.push(SimpleCigar::Del(len));
            }
            Cigar::RefSkip(len) => {
                ref_len += u64::from(len);
                cigar.push(SimpleCigar::RefSkip(len));
            }
            Cigar::SoftClip(len) => {
                read_pos += len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    if clipped_bases.is_empty() || ref_len == 0 {
        return Ok(None);
    }
    let start0 = record.pos() as u64;
    let end0 = start0 + ref_len - 1;

    Ok(Some(PreparedRead {
        contig,
        start0,
        end0,
        bases: clipped_bases,
        quals: clipped_quals,
        cigar,
        rg_index,
        is_reverse: record.is_reverse(),
        is_second_in_pair: record.is_paired() && record.is_last_in_template(),
    }))
}

fn copy_read_segment(
    bases: &[u8],
    quals: &[u8],
    read_pos: usize,
    len: u32,
    clipped_bases: &mut Vec<u8>,
    clipped_quals: &mut Vec<u8>,
) {
    let end = read_pos + len as usize;
    clipped_bases.extend_from_slice(&bases[read_pos..end]);
    clipped_quals.extend_from_slice(&quals[read_pos..end]);
}

fn aux_string(record: &bam::Record, tag: &[u8]) -> Result<String> {
    match record.aux(tag) {
        Ok(Aux::String(value)) => Ok(value.to_string()),
        Ok(Aux::Char(value)) => Ok((value as char).to_string()),
        Ok(other) => bail!(
            "aux tag {} has unsupported value type {:?}",
            String::from_utf8_lossy(tag),
            other
        ),
        Err(_) => bail!(
            "read {} is missing required aux tag {}",
            String::from_utf8_lossy(record.qname()),
            String::from_utf8_lossy(tag)
        ),
    }
}

fn original_qualities(record: &bam::Record) -> Result<Option<Vec<u8>>> {
    match record.aux(b"OQ") {
        Ok(Aux::String(value)) => Ok(Some(value.bytes().map(|b| b.saturating_sub(33)).collect())),
        Ok(_) => bail!(
            "read {} has non-string OQ tag",
            String::from_utf8_lossy(record.qname())
        ),
        Err(_) => Ok(None),
    }
}

impl ReferenceCache {
    fn new(reference: &PathBuf, chunk_size: u64) -> Result<Self> {
        let reader = faidx::Reader::from_path(reference)
            .with_context(|| format!("failed to open FASTA index for {}", reference.display()))?;
        Ok(Self {
            reader,
            contig: None,
            start0: 0,
            end0: 0,
            bases: Vec::new(),
            chunk_size,
        })
    }

    fn bases_for(&mut self, contig: &str, start0: u64, end0: u64) -> Result<&[u8]> {
        let cache_hit =
            self.contig.as_deref() == Some(contig) && start0 >= self.start0 && end0 <= self.end0;
        if !cache_hit {
            self.load_chunk(contig, start0, end0)?;
        }
        let offset = usize::try_from(start0 - self.start0)
            .with_context(|| format!("reference offset overflow for {contig}:{start0}-{end0}"))?;
        let length = usize::try_from(end0 - start0 + 1)
            .with_context(|| format!("reference length overflow for {contig}:{start0}-{end0}"))?;
        Ok(&self.bases[offset..offset + length])
    }

    fn load_chunk(&mut self, contig: &str, start0: u64, end0: u64) -> Result<()> {
        self.contig = Some(contig.to_string());
        self.start0 = start0;
        let contig_len = self.reader.fetch_seq_len(contig);
        if contig_len == 0 {
            bail!("contig '{contig}' is not present in the reference FASTA");
        }
        self.end0 = end0
            .max(start0.saturating_add(self.chunk_size).saturating_sub(1))
            .min(contig_len - 1);
        self.bases = self
            .reader
            .fetch_seq_string(contig, self.start0 as usize, self.end0 as usize)
            .with_context(|| {
                format!(
                    "failed to fetch reference {}:{}-{}",
                    contig,
                    self.start0 + 1,
                    self.end0 + 1
                )
            })?
            .into_bytes();
        Ok(())
    }
}

impl KnownSiteCache {
    fn new(paths: &[PathBuf], chunk_size: u64) -> Result<Self> {
        let mut readers = Vec::new();
        for path in paths {
            let reader = bcf::IndexedReader::from_path(path).with_context(|| {
                format!("failed to open indexed known-sites VCF {}", path.display())
            })?;
            readers.push(KnownSiteReader {
                path: path.clone(),
                reader,
                unfetchable_contigs: HashSet::new(),
            });
        }
        Ok(Self {
            readers,
            contig: None,
            start0: 0,
            end0: 0,
            intervals: Vec::new(),
            chunk_size,
        })
    }

    fn intervals_for(&mut self, contig: &str, start0: u64, end0: u64) -> Result<&[KnownInterval]> {
        let cache_hit =
            self.contig.as_deref() == Some(contig) && start0 >= self.start0 && end0 <= self.end0;
        if !cache_hit {
            self.load_chunk(contig, start0, end0)?;
        }
        Ok(&self.intervals)
    }

    fn load_chunk(&mut self, contig: &str, start0: u64, end0: u64) -> Result<()> {
        self.contig = Some(contig.to_string());
        self.start0 = start0;
        self.end0 = end0.max(start0.saturating_add(self.chunk_size).saturating_sub(1));
        self.intervals.clear();

        for known_reader in &mut self.readers {
            if known_reader.unfetchable_contigs.contains(contig) {
                continue;
            }
            let rid = match known_reader.reader.header().name2rid(contig.as_bytes()) {
                Ok(rid) => rid,
                Err(_) => continue,
            };
            if let Err(err) = known_reader
                .reader
                .fetch(rid, self.start0, Some(self.end0 + 1))
            {
                eprintln!(
                    "warning: skipping known-sites VCF {} on unfetchable contig {}: {}",
                    known_reader.path.display(),
                    contig,
                    err
                );
                known_reader.unfetchable_contigs.insert(contig.to_string());
                continue;
            }
            for record_result in known_reader.reader.records() {
                let record = record_result.with_context(|| {
                    format!(
                        "failed to read known-sites VCF {}",
                        known_reader.path.display()
                    )
                })?;
                let rec_start = record.pos().max(0) as u64;
                let rec_len = record.rlen().max(1) as u64;
                let rec_end = rec_start + rec_len - 1;
                if rec_end >= self.start0 && rec_start <= self.end0 {
                    self.intervals.push(KnownInterval {
                        start0: rec_start,
                        end0: rec_end,
                    });
                }
            }
        }
        self.intervals
            .sort_by(|a, b| a.start0.cmp(&b.start0).then(a.end0.cmp(&b.end0)));
        Ok(())
    }
}


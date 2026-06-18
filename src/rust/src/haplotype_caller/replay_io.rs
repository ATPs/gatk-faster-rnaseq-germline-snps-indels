fn sort_replay_rows(output: &mut ReplayWorkerOutput, dict: &SequenceDict) -> Result<()> {
    for row in &output.active_regions {
        if dict.order(&row.contig).is_none() {
            bail!(
                "contig '{}' is not present in the sequence dictionary",
                row.contig
            );
        }
    }
    for row in &output.active_loci {
        if dict.order(&row.contig).is_none() {
            bail!(
                "contig '{}' is not present in the sequence dictionary",
                row.contig
            );
        }
    }
    for row in &output.events {
        if dict.order(&row.chrom).is_none() {
            bail!(
                "contig '{}' is not present in the sequence dictionary",
                row.chrom
            );
        }
    }
    output.active_regions.sort_by(|a, b| {
        dict.order(&a.contig)
            .cmp(&dict.order(&b.contig))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    output.active_loci.sort_by(|a, b| {
        dict.order(&a.contig)
            .cmp(&dict.order(&b.contig))
            .then(a.pos.cmp(&b.pos))
            .then(a.region.cmp(&b.region))
    });
    output.read_observations.sort_by(|a, b| {
        a.region
            .cmp(&b.region)
            .then(a.pos.cmp(&b.pos))
            .then(a.read.cmp(&b.read))
            .then(a.kind.cmp(&b.kind))
            .then(a.qpos.cmp(&b.qpos))
    });
    output.events.sort_by(|a, b| {
        dict.order(&a.chrom)
            .cmp(&dict.order(&b.chrom))
            .then(a.pos.cmp(&b.pos))
            .then(a.event_type.cmp(&b.event_type))
            .then(a.alleles.cmp(&b.alleles))
    });
    output.haplotypes.sort_by(|a, b| {
        a.region
            .cmp(&b.region)
            .then(a.stage.cmp(&b.stage))
            .then(a.haplotype.cmp(&b.haplotype))
    });
    output.pairhmms.sort_by(|a, b| {
        a.region
            .cmp(&b.region)
            .then(a.read_index.cmp(&b.read_index))
            .then(a.haplotype.cmp(&b.haplotype))
    });
    Ok(())
}

fn write_replay_tables(
    prefix: &Path,
    output: &ReplayWorkerOutput,
    genotype_rows: &[ReplayGenotypeRow],
) -> Result<()> {
    write_replay_active_regions(&replay_prefixed_path(prefix, "active_regions.tsv"), output)?;
    write_replay_active_loci(&replay_prefixed_path(prefix, "active_loci.tsv"), output)?;
    write_replay_read_observations(
        &replay_prefixed_path(prefix, "read_observations.tsv"),
        output,
    )?;
    write_replay_events(&replay_prefixed_path(prefix, "events.tsv"), output)?;
    write_replay_genotypes(
        &replay_prefixed_path(prefix, "genotypes.tsv"),
        genotype_rows,
    )?;
    write_replay_haplotypes(&replay_prefixed_path(prefix, "haplotypes.tsv"), output)?;
    write_replay_pairhmms(&replay_prefixed_path(prefix, "pairhmm.tsv"), output)?;
    write_empty_allele_likelihoods(prefix)?;
    Ok(())
}

fn write_replay_active_regions(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "contig\tstart\tend\tregion\tactive_start\tactive_end\tpadded_start\tpadded_end\tobserved_loci\tactive_loci\tcandidate_events\tmax_alt_fraction\tmean_alt_fraction"
    )?;
    for row in &output.active_regions {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}",
            row.contig,
            row.start,
            row.end,
            row.region,
            row.active_start,
            row.active_end,
            row.padded_start,
            row.padded_end,
            row.observed_loci,
            row.active_loci,
            row.candidate_events,
            row.max_alt_fraction,
            row.mean_alt_fraction
        )?;
    }
    Ok(())
}

fn write_replay_active_loci(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "contig\tpos\tregion\tref_base\tdepth\tsnp_alt_count\tsnp_best_alt\tindel_alt_count\tindel_best_alt\talt_fraction\tactive_probability_proxy"
    )?;
    for row in &output.active_loci {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}",
            row.contig,
            row.pos,
            row.region,
            row.ref_base as char,
            row.depth,
            row.snp_alt_count,
            row.snp_best_alt,
            row.indel_alt_count,
            row.indel_best_alt,
            row.alt_fraction,
            row.active_probability_proxy
        )?;
    }
    Ok(())
}

fn write_replay_read_observations(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "region\tread\tkind\tpos\tqpos\tallele\tadjusted_quality\tmapq\tstrand"
    )?;
    for row in &output.read_observations {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.region,
            row.read,
            row.kind,
            row.pos,
            row.qpos,
            row.allele,
            row.adjusted_quality,
            row.mapq,
            row.strand
        )?;
    }
    Ok(())
}

fn write_replay_events(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "region\tevent\tchrom\tpos\ttype\talleles\traw\tdepth\tref_count\talt_count\tqual\tgt"
    )?;
    for row in &output.events {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.region,
            row.event,
            row.chrom,
            row.pos,
            row.event_type,
            row.alleles,
            row.raw,
            row.depth,
            row.ref_count,
            row.alt_count,
            row.qual,
            row.gt
        )?;
    }
    Ok(())
}

fn write_replay_genotypes(path: &Path, rows: &[ReplayGenotypeRow]) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "chrom\tpos\tref\talt\tqual\tfilter\tgt\tgq\tdp\tad_ref\tad_alt\tfs\tqd\tpl\tdb"
    )?;
    for row in rows {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.2}\t{}\t{}",
            row.chrom,
            row.pos,
            row.ref_allele,
            row.alt,
            row.qual,
            row.filter,
            row.gt,
            row.gq,
            row.dp,
            row.ad_ref,
            row.ad_alt,
            row.fs,
            row.qd,
            row.pl,
            row.db
        )?;
    }
    Ok(())
}
fn write_replay_haplotypes(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "region\tstage\thaplotype\tspan_start\tspan_end\tkmer\tlength\tcigar\tis_ref\tbases"
    )?;
    for row in &output.haplotypes {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.region,
            row.stage,
            row.haplotype,
            row.span_start,
            row.span_end,
            row.kmer,
            row.length,
            row.cigar,
            row.is_ref,
            row.bases
        )?;
    }
    Ok(())
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct LocalHaplotype {
    bases: Vec<u8>,
    is_ref: bool,
    cigar: String,
    event_indices: Vec<usize>,
}

fn push_cigar(cigar: &mut Vec<(u32, char)>, len: u32, op: char) {
    if len == 0 {
        return;
    }
    if let Some(last) = cigar.last_mut() {
        if last.1 == op {
            last.0 += len;
            return;
        }
    }
    cigar.push((len, op));
}

fn format_cigar(cigar: &[(u32, char)]) -> String {
    let mut s = String::with_capacity(cigar.len() * 4);
    for &(len, op) in cigar {
        s.push_str(&len.to_string());
        s.push(op);
    }
    s
}

fn build_local_haplotypes(
    _contig: &str,
    region_start: u64,
    region_end: u64,
    ref_bases: &[u8],
    candidate_events: &[VariantCall],
    max_haplotypes: usize,
) -> Vec<LocalHaplotype> {
    let ref_hap = LocalHaplotype {
        bases: ref_bases.to_vec(),
        is_ref: true,
        cigar: format!("{}M", ref_bases.len()),
        event_indices: Vec::new(),
    };

    let mut haplotypes = vec![ref_hap];

    let mut valid_events = Vec::new();
    for event in candidate_events {
        let event_end = event.pos + event.ref_allele.len() as u64 - 1;
        if event.pos >= region_start && event_end <= region_end {
            valid_events.push(event);
        }
    }

    if valid_events.len() > 7 {
        valid_events.truncate(7);
    }

    let n_events = valid_events.len();
    for mask in 1..(1 << n_events) {
        if haplotypes.len() >= max_haplotypes {
            break;
        }

        let mut overlap = false;
        let mut last_end = 0;
        let mut selected_events = Vec::new();
        let mut event_indices = Vec::new();
        for i in 0..n_events {
            if (mask & (1 << i)) != 0 {
                let ev = valid_events[i];
                if ev.pos <= last_end {
                    overlap = true;
                    break;
                }
                last_end = ev.pos + ev.ref_allele.len() as u64 - 1;
                selected_events.push(ev);
                event_indices.push(i);
            }
        }

        if overlap {
            continue;
        }

        let mut bases = Vec::with_capacity(ref_bases.len());
        let mut cigar_ops = Vec::new();
        let mut ref_offset = 0;
        let mut current_pos = region_start;

        for ev in &selected_events {
            let dist = ev.pos.saturating_sub(current_pos) as usize;
            if dist > 0 {
                bases.extend_from_slice(&ref_bases[ref_offset..ref_offset + dist]);
                push_cigar(&mut cigar_ops, dist as u32, 'M');
                ref_offset += dist;
                current_pos += dist as u64;
            }

            bases.extend_from_slice(&ev.alt_allele);
            let match_len = (ev.alt_allele.len() as u32).min(ev.ref_allele.len() as u32);
            push_cigar(&mut cigar_ops, match_len, 'M');

            if ev.alt_allele.len() > ev.ref_allele.len() {
                push_cigar(
                    &mut cigar_ops,
                    (ev.alt_allele.len() - ev.ref_allele.len()) as u32,
                    'I',
                );
            } else if ev.ref_allele.len() > ev.alt_allele.len() {
                push_cigar(
                    &mut cigar_ops,
                    (ev.ref_allele.len() - ev.alt_allele.len()) as u32,
                    'D',
                );
            }

            ref_offset += ev.ref_allele.len();
            current_pos += ev.ref_allele.len() as u64;
        }

        let rem = ref_bases.len().saturating_sub(ref_offset);
        if rem > 0 {
            bases.extend_from_slice(&ref_bases[ref_offset..]);
            push_cigar(&mut cigar_ops, rem as u32, 'M');
        }

        haplotypes.push(LocalHaplotype {
            bases,
            is_ref: false,
            cigar: format_cigar(&cigar_ops),
            event_indices,
        });
    }

    haplotypes
}

fn write_replay_pairhmms(path: &Path, output: &ReplayWorkerOutput) -> Result<()> {
    create_parent_dir(path)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(
        writer,
        "region\tread\thaplotype\tread_index\tcigar\tmapq\tloc\tunclipped_loc\tlength\tscore"
    )?;
    for row in &output.pairhmms {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.region,
            row.read,
            row.haplotype,
            row.read_index,
            row.cigar,
            row.mapq,
            row.loc,
            row.unclipped_loc,
            row.length,
            row.score
        )?;
    }
    Ok(())
}

fn write_empty_allele_likelihoods(prefix: &Path) -> Result<()> {
    let mut allele_likelihoods = BufWriter::new(File::create(replay_prefixed_path(
        prefix,
        "allele_likelihoods.tsv",
    ))?);
    writeln!(
        allele_likelihoods,
        "region\tevent\tmatrix\tread\tread_index\tallele\tscore"
    )?;
    Ok(())
}

fn replay_prefixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    PathBuf::from(format!("{}.{}", prefix.display(), suffix))
}

fn open_vcf_writer(path: &Path) -> Result<Box<dyn Write>> {
    if is_gzip_path(path) {
        let writer = bgzf::Writer::from_path(path)
            .with_context(|| format!("creating bgzipped VCF {}", path.display()))?;
        Ok(Box::new(writer))
    } else {
        let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
        Ok(Box::new(BufWriter::new(file)))
    }
}

fn write_variant_record(writer: &mut dyn Write, variant: &VariantCall) -> Result<()> {
    let ac = variant.alt_allele_count();
    let af = f64::from(ac) / 2.0;
    let qd = if variant.depth == 0 {
        0.0
    } else {
        fix_too_high_qd(f64::from(variant.qual) / f64::from(variant.depth))
    };
    let gt = variant.genotype();
    let pl = genotype_likelihoods(variant);
    let ref_allele = allele_string(&variant.ref_allele)?;
    let alt_allele = allele_string(&variant.alt_allele)?;
    let id = variant.id.as_deref().unwrap_or(".");
    let db = if variant.db { ";DB" } else { "" };
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\tPASS\tAC={};AF={:.3};AN=2;DP={}{};FS={:.3};QD={:.2}\tGT:AD:DP:GQ:PL\t{}:{},{}:{}:{}:{}",
        variant.contig,
        variant.pos,
        id,
        ref_allele,
        alt_allele,
        variant.qual,
        ac,
        af,
        variant.depth,
        db,
        variant.fs,
        qd,
        gt,
        variant.ref_count,
        variant.alt_count,
        variant.depth,
        variant.gq(),
        pl
    )?;
    Ok(())
}

fn genotype_likelihoods(variant: &VariantCall) -> String {
    format!("{},{},{}", variant.pl[0], variant.pl[1], variant.pl[2])
}

fn fix_too_high_qd(qd: f64) -> f64 {
    if qd < 35.0 {
        qd
    } else {
        30.0
    }
}

fn allele_string(allele: &[u8]) -> Result<&str> {
    std::str::from_utf8(allele).context("allele contains non-UTF-8 bases")
}

fn sample_name_from_bam(path: &Path) -> Result<String> {
    let reader =
        bam::Reader::from_path(path).with_context(|| format!("opening {}", path.display()))?;
    let header = String::from_utf8_lossy(reader.header().as_bytes());
    for line in header.lines() {
        if !line.starts_with("@RG\t") {
            continue;
        }
        for field in line.split('\t').skip(1) {
            if let Some(sample) = field.strip_prefix("SM:") {
                if !sample.is_empty() {
                    return Ok(sample.to_string());
                }
            }
        }
    }
    Ok("SAMPLE".to_string())
}

fn is_gzip_path(path: &Path) -> bool {
    path.extension().is_some_and(|extension| extension == "gz")
}

fn write_tabix_index(path: &Path, threads: usize) -> Result<()> {
    let c_path = CString::new(path.as_os_str().as_bytes())?;
    let thread_count = threads.clamp(1, i32::MAX as usize) as i32;
    let result = unsafe {
        htslib::tbx_index_build3(
            c_path.as_ptr(),
            ptr::null(),
            0,
            thread_count,
            &htslib::tbx_conf_vcf,
        )
    };
    if result != 0 {
        bail!("failed to create tabix index for {}", path.display());
    }
    Ok(())
}

fn create_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)
                .with_context(|| format!("creating directory {}", parent.display()))?;
        }
    }
    Ok(())
}

fn normalize_base(base: u8) -> u8 {
    base.to_ascii_uppercase()
}

fn is_acgt(base: u8) -> bool {
    matches!(base, b'A' | b'C' | b'G' | b'T')
}

fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn base_from_index(index: usize) -> u8 {
    match index {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => unreachable!("invalid base index"),
    }
}

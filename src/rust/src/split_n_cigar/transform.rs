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


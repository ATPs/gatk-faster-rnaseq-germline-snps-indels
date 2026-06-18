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


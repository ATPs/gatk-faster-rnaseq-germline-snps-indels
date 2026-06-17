use rust_htslib::bam::record::{Cigar, CigarString};
use std::cmp::max;

pub struct SWParameters {
    pub match_value: i32,
    pub mismatch_penalty: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
}

impl Default for SWParameters {
    fn default() -> Self {
        Self {
            match_value: 1,
            mismatch_penalty: -4,
            gap_open_penalty: -6,
            gap_extend_penalty: -1,
        }
    }
}

#[derive(PartialEq)]
pub enum SWOverhangStrategy {
    SoftClip,
    Ignore,
    Indel,
    LeadingIndel,
}

pub struct SWAlignmentResult {
    pub cigar: CigarString,
    pub alignment_offset: i32,
}

#[derive(Clone, Copy, PartialEq)]
enum State {
    Match,
    Insertion,
    Deletion,
    Clip,
}

pub fn align(
    reference: &[u8],
    alternate: &[u8],
    parameters: &SWParameters,
    overhang_strategy: SWOverhangStrategy,
) -> SWAlignmentResult {
    if reference.is_empty() || alternate.is_empty() {
        panic!("Non-null, non-empty sequences are required for Smith-Waterman");
    }

    // Try substring matching for softclip/ignore
    if overhang_strategy == SWOverhangStrategy::SoftClip || overhang_strategy == SWOverhangStrategy::Ignore {
        if let Some(pos) = find_subsequence(reference, alternate) {
            return SWAlignmentResult {
                cigar: CigarString(vec![Cigar::Match(alternate.len() as u32)]),
                alignment_offset: pos as i32,
            };
        }
    }

    let n = reference.len() + 1;
    let m = alternate.len() + 1;

    let mut sw = vec![vec![0; m]; n];
    let mut btrack = vec![vec![0; m]; n];

    calculate_matrix(reference, alternate, &mut sw, &mut btrack, &overhang_strategy, parameters);
    calculate_cigar(&sw, &btrack, &overhang_strategy)
}

fn find_subsequence(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() { return Some(0); }
    if haystack.len() < needle.len() { return None; }
    // find last index to match java's lastIndexOf, though usually first is fine.
    // Let's match lastIndexOf just in case.
    for i in (0..=haystack.len() - needle.len()).rev() {
        if &haystack[i..i+needle.len()] == needle {
            return Some(i);
        }
    }
    None
}

fn calculate_matrix(
    reference: &[u8],
    alternate: &[u8],
    sw: &mut [Vec<i32>],
    btrack: &mut [Vec<i32>],
    overhang_strategy: &SWOverhangStrategy,
    parameters: &SWParameters,
) {
    let ncol = alternate.len() + 1;
    let nrow = reference.len() + 1;

    let matrix_min_cutoff: i32 = -100_000_000;
    let low_init_value = i32::MIN / 2;

    let mut best_gap_v = vec![low_init_value; ncol + 1];
    let mut gap_size_v = vec![0; ncol + 1];
    let mut best_gap_h = vec![low_init_value; nrow + 1];
    let mut gap_size_h = vec![0; nrow + 1];

    if *overhang_strategy == SWOverhangStrategy::Indel || *overhang_strategy == SWOverhangStrategy::LeadingIndel {
        sw[0][1] = parameters.gap_open_penalty;
        let mut cur_val = parameters.gap_open_penalty;
        for i in 2..ncol {
            cur_val += parameters.gap_extend_penalty;
            sw[0][i] = cur_val;
        }

        sw[1][0] = parameters.gap_open_penalty;
        cur_val = parameters.gap_open_penalty;
        for i in 2..nrow {
            cur_val += parameters.gap_extend_penalty;
            sw[i][0] = cur_val;
        }
    }

    let w_open = parameters.gap_open_penalty;
    let w_extend = parameters.gap_extend_penalty;
    let w_match = parameters.match_value;
    let w_mismatch = parameters.mismatch_penalty;

    for i in 1..nrow {
        let a_base = reference[i - 1];
        for j in 1..ncol {
            let b_base = alternate[j - 1];
            
            let step_diag = sw[i - 1][j - 1] + if a_base == b_base { w_match } else { w_mismatch };

            let mut prev_gap = sw[i - 1][j] + w_open;
            best_gap_v[j] = best_gap_v[j].saturating_add(w_extend);
            if prev_gap > best_gap_v[j] {
                best_gap_v[j] = prev_gap;
                gap_size_v[j] = 1;
            } else {
                gap_size_v[j] += 1;
            }
            let step_down = best_gap_v[j];
            let kd = gap_size_v[j];

            prev_gap = sw[i][j - 1] + w_open;
            best_gap_h[i] = best_gap_h[i].saturating_add(w_extend);
            if prev_gap > best_gap_h[i] {
                best_gap_h[i] = prev_gap;
                gap_size_h[i] = 1;
            } else {
                gap_size_h[i] += 1;
            }
            let step_right = best_gap_h[i];
            let ki = gap_size_h[i];

            let diag_highest_or_equal = step_diag >= step_down && step_diag >= step_right;

            if diag_highest_or_equal {
                sw[i][j] = max(matrix_min_cutoff, step_diag);
                btrack[i][j] = 0;
            } else if step_right >= step_down {
                sw[i][j] = max(matrix_min_cutoff, step_right);
                btrack[i][j] = -ki;
            } else {
                sw[i][j] = max(matrix_min_cutoff, step_down);
                btrack[i][j] = kd;
            }
        }
    }
}

fn calculate_cigar(
    sw: &[Vec<i32>],
    btrack: &[Vec<i32>],
    overhang_strategy: &SWOverhangStrategy,
) -> SWAlignmentResult {
    let ref_length = sw.len() - 1;
    let alt_length = sw[0].len() - 1;

    let mut p1 = 0;
    let mut p2 = 0;
    let mut maxscore = i32::MIN;
    let mut segment_length = 0;

    if *overhang_strategy == SWOverhangStrategy::Indel {
        p1 = ref_length;
        p2 = alt_length;
    } else {
        p2 = alt_length;
        for i in 1..=ref_length {
            let cur_score = sw[i][alt_length];
            if cur_score >= maxscore {
                p1 = i;
                maxscore = cur_score;
            }
        }
        
        if *overhang_strategy != SWOverhangStrategy::LeadingIndel {
            for j in 1..=alt_length {
                let cur_score = sw[ref_length][j];
                if cur_score > maxscore || (cur_score == maxscore && (ref_length as i32 - j as i32).abs() < (p1 as i32 - p2 as i32).abs()) {
                    p1 = ref_length;
                    p2 = j;
                    maxscore = cur_score;
                    segment_length = alt_length - j;
                }
            }
        }
    }

    let mut lce = Vec::new();
    if segment_length > 0 && *overhang_strategy == SWOverhangStrategy::SoftClip {
        lce.push(make_element(State::Clip, segment_length));
        segment_length = 0;
    }

    let mut state = State::Match;
    
    while p1 > 0 && p2 > 0 {
        let btr = btrack[p1][p2];
        let new_state;
        let mut step_length = 1;
        
        if btr > 0 {
            new_state = State::Deletion;
            step_length = btr as usize;
        } else if btr < 0 {
            new_state = State::Insertion;
            step_length = (-btr) as usize;
        } else {
            new_state = State::Match;
        }

        match new_state {
            State::Match => { p1 -= 1; p2 -= 1; },
            State::Insertion => { p2 -= step_length; },
            State::Deletion => { p1 -= step_length; },
            _ => unreachable!(),
        }

        if new_state == state {
            segment_length += step_length;
        } else {
            if segment_length > 0 {
                lce.push(make_element(state, segment_length));
            }
            segment_length = step_length;
            state = new_state;
        }
    }

    let alignment_offset;

    if *overhang_strategy == SWOverhangStrategy::SoftClip {
        lce.push(make_element(state, segment_length));
        if p2 > 0 {
            lce.push(make_element(State::Clip, p2));
        }
        alignment_offset = p1 as i32;
    } else if *overhang_strategy == SWOverhangStrategy::Ignore {
        lce.push(make_element(state, segment_length + p2));
        alignment_offset = p1 as i32 - p2 as i32;
    } else {
        lce.push(make_element(state, segment_length));
        if p1 > 0 {
            lce.push(make_element(State::Deletion, p1));
        } else if p2 > 0 {
            lce.push(make_element(State::Insertion, p2));
        }
        alignment_offset = 0;
    }

    lce.reverse();
    SWAlignmentResult {
        cigar: CigarString(lce),
        alignment_offset,
    }
}

fn make_element(state: State, length: usize) -> Cigar {
    let len = length as u32;
    match state {
        State::Match => Cigar::Match(len),
        State::Insertion => Cigar::Ins(len),
        State::Deletion => Cigar::Del(len),
        State::Clip => Cigar::SoftClip(len),
    }
}

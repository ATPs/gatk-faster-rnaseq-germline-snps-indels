pub const TRISTATE_CORRECTION: f64 = 3.0;

#[inline]
pub fn qual_to_error_prob(q: u8) -> f64 {
    10.0f64.powf(-(q as f64) / 10.0)
}

#[inline]
pub fn qual_to_prob(q: u8) -> f64 {
    1.0 - qual_to_error_prob(q)
}

pub struct PairHmmTransition {
    pub m2m: f64,
    pub m2i: f64,
    pub m2d: f64,
    pub i2m: f64,
    pub i2i: f64,
    pub d2m: f64,
    pub d2d: f64,
}

pub fn qual_to_trans_probs(ins_qual: u8, del_qual: u8, gcp: u8) -> PairHmmTransition {
    let m2i = qual_to_error_prob(ins_qual);
    let m2d = qual_to_error_prob(del_qual);
    let m2m = 1.0 - m2i - m2d;
    let i2i = qual_to_error_prob(gcp);
    let i2m = 1.0 - i2i;
    let d2d = qual_to_error_prob(gcp);
    let d2m = 1.0 - d2d;

    PairHmmTransition {
        m2m: m2m.max(0.0), // Prevent negative from approximation limits
        m2i,
        m2d,
        i2m,
        i2i,
        d2m,
        d2d,
    }
}

pub fn compute_read_likelihood_given_haplotype(
    hap_bases: &[u8],
    read_bases: &[u8],
    read_quals: &[u8],
    read_ins_quals: &[u8],
    read_del_quals: &[u8],
    gcp: u8,
) -> f64 {
    let read_len = read_bases.len();
    let hap_len = hap_bases.len();

    let mut transitions = Vec::with_capacity(read_len);
    for i in 0..read_len {
        transitions.push(qual_to_trans_probs(read_ins_quals[i], read_del_quals[i], gcp));
    }

    let initial_value = 1.0 / (hap_len as f64);

    // M, I, D matrices. Indexed by [read_pos][hap_pos]
    let mut match_matrix = vec![vec![0.0; hap_len + 1]; read_len + 1];
    let mut insertion_matrix = vec![vec![0.0; hap_len + 1]; read_len + 1];
    let mut deletion_matrix = vec![vec![0.0; hap_len + 1]; read_len + 1];

    for j in 0..=hap_len {
        deletion_matrix[0][j] = initial_value;
    }

    for i in 1..=read_len {
        let t = &transitions[i - 1];
        let rb = read_bases[i - 1];
        let rq = read_quals[i - 1];
        let p_match = qual_to_prob(rq);
        let p_mismatch = qual_to_error_prob(rq) / TRISTATE_CORRECTION;

        for j in 1..=hap_len {
            let hb = hap_bases[j - 1];
            let prior = if rb == hb || rb == b'N' || hb == b'N' {
                p_match
            } else {
                p_mismatch
            };

            match_matrix[i][j] = prior * (
                match_matrix[i - 1][j - 1] * t.m2m +
                insertion_matrix[i - 1][j - 1] * t.i2m +
                deletion_matrix[i - 1][j - 1] * t.d2m
            );
            
            insertion_matrix[i][j] = match_matrix[i - 1][j] * t.m2i + insertion_matrix[i - 1][j] * t.i2i;
            deletion_matrix[i][j] = match_matrix[i][j - 1] * t.m2d + deletion_matrix[i][j - 1] * t.d2d;
        }
    }

    let mut final_sum = 0.0;
    for j in 1..=hap_len {
        final_sum += match_matrix[read_len][j] + insertion_matrix[read_len][j];
    }

    if final_sum > 0.0 {
        final_sum.log10()
    } else {
        -1000.0 // log10(0) approximation
    }
}

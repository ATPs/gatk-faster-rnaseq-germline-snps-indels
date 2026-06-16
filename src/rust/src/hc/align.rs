#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SmithWatermanScoring {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
}

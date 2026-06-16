#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReferenceConfidenceMode {
    None,
    Gvcf,
    BpResolution,
}

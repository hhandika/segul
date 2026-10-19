//! Pairwise genetic distances between aligned sequences.
//!
//! Supported models:
//! 1. Uncorrected p-distance (DNA and amino acid)
//! 2. Jukes-Cantor (JC69; DNA and amino acid)
//! 3. Kimura 2-parameter (K2P; DNA)
//! 4. Tamura-Nei (TN93; DNA)
//! 5. Poisson correction (amino acid)
//!
//! Missing data are handled using pairwise deletion.
//! A site is skipped for a pair if either sequence has a gap,
//! missing data, or an ambiguous character at that site.
use crate::helper::types::DataType;

/// Encoded value for gaps, missing data, and ambiguous characters.
pub const INVALID: u8 = u8::MAX;

// DNA encoding. Transitions are A <-> G (purines) and C <-> T (pyrimidines).
const A: u8 = 0;
const C: u8 = 1;
const G: u8 = 2;
const T: u8 = 3;

/// The 20 standard amino acid residues.
const AA_RESIDUES: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

const DNA_TABLE: [u8; 256] = build_dna_table();
const AA_TABLE: [u8; 256] = build_aa_table();

const fn build_dna_table() -> [u8; 256] {
    let mut table = [INVALID; 256];
    let bases = [(b'A', A), (b'C', C), (b'G', G), (b'T', T), (b'U', T)];
    let mut i = 0;
    while i < bases.len() {
        let (base, code) = bases[i];
        table[base as usize] = code;
        table[base.to_ascii_lowercase() as usize] = code;
        i += 1;
    }
    table
}

const fn build_aa_table() -> [u8; 256] {
    let mut table = [INVALID; 256];
    let mut i = 0;
    while i < AA_RESIDUES.len() {
        let residue = AA_RESIDUES[i];
        table[residue as usize] = i as u8;
        table[residue.to_ascii_lowercase() as usize] = i as u8;
        i += 1;
    }
    table
}

/// Encode a sequence for distance computation.
/// Unambiguous characters are encoded as indices,
/// everything else (gaps, missing data, ambiguous characters) as [`INVALID`].
/// Matching is case-insensitive. For DNA, `U` is treated as `T`.
pub fn encode(seq: &str, datatype: &DataType) -> Vec<u8> {
    let table = match datatype {
        DataType::Aa => &AA_TABLE,
        _ => &DNA_TABLE,
    };
    seq.bytes().map(|b| table[b as usize]).collect()
}

/// Distance models.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistanceModel {
    /// Uncorrected p-distance: the proportion of differing sites.
    Uncorrected,
    /// Jukes-Cantor (1969) correction.
    JukesCantor,
    /// Kimura (1980) 2-parameter correction. DNA only.
    Kimura2P,
    /// Tamura-Nei (1993) correction. DNA only.
    /// Base frequencies are estimated from the whole alignment.
    TamuraNei,
    /// Poisson correction. Amino acid only.
    Poisson,
}

impl DistanceModel {
    /// Check whether the model can be used for the data type.
    pub fn is_valid_for(&self, datatype: &DataType) -> bool {
        match datatype {
            DataType::Dna => !matches!(self, Self::Poisson),
            DataType::Aa => matches!(self, Self::Uncorrected | Self::JukesCantor | Self::Poisson),
            DataType::Ignore => false,
        }
    }

    /// Compute the distance from pairwise site counts.
    /// Returns `None` if no sites are compared
    /// or if the distance is undefined (saturated).
    pub fn distance(
        &self,
        counts: &PairCounts,
        datatype: &DataType,
        freqs: &BaseFreqs,
    ) -> Option<f64> {
        if counts.sites == 0 {
            return None;
        }
        let sites = counts.sites as f64;
        let p = counts.diffs as f64 / sites;
        let d = match self {
            Self::Uncorrected => Some(p),
            Self::JukesCantor => {
                let b = match datatype {
                    DataType::Aa => 19.0 / 20.0,
                    _ => 3.0 / 4.0,
                };
                neg_ln(1.0 - p / b).map(|v| b * v)
            }
            Self::Kimura2P => {
                let ts = (counts.ts_purine + counts.ts_pyrimidine) as f64 / sites;
                let tv = counts.transversions as f64 / sites;
                Some(0.5 * neg_ln(1.0 - 2.0 * ts - tv)? + 0.25 * neg_ln(1.0 - 2.0 * tv)?)
            }
            Self::TamuraNei => {
                let p1 = counts.ts_purine as f64 / sites;
                let p2 = counts.ts_pyrimidine as f64 / sites;
                let q = counts.transversions as f64 / sites;
                freqs.tamura_nei(p1, p2, q)
            }
            Self::Poisson => neg_ln(1.0 - p),
        };
        // Avoid writing negative zero
        d.map(|d| if d <= 0.0 { 0.0 } else { d })
    }
}

impl std::fmt::Display for DistanceModel {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::Uncorrected => write!(f, "p-distance"),
            Self::JukesCantor => write!(f, "JC69"),
            Self::Kimura2P => write!(f, "K2P"),
            Self::TamuraNei => write!(f, "TN93"),
            Self::Poisson => write!(f, "Poisson"),
        }
    }
}

impl std::str::FromStr for DistanceModel {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "p" => Ok(Self::Uncorrected),
            "jc69" => Ok(Self::JukesCantor),
            "k2p" => Ok(Self::Kimura2P),
            "tn93" => Ok(Self::TamuraNei),
            "poisson" => Ok(Self::Poisson),
            _ => Err(format!("{} is not a valid distance model", s)),
        }
    }
}

/// Negative natural log. `None` if the argument is not positive.
fn neg_ln(x: f64) -> Option<f64> {
    if x > 0.0 { Some(-x.ln()) } else { None }
}

/// Site counts for a pair of encoded sequences.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct PairCounts {
    /// Number of sites compared (both characters valid).
    pub sites: usize,
    /// Number of differing sites.
    pub diffs: usize,
    /// A <-> G transitions. DNA only.
    pub ts_purine: usize,
    /// C <-> T transitions. DNA only.
    pub ts_pyrimidine: usize,
    /// Transversions. DNA only.
    pub transversions: usize,
}

impl PairCounts {
    /// Count compared and differing sites between two encoded sequences.
    /// Sites with an [`INVALID`] character in either sequence are skipped.
    pub fn new(seq_1: &[u8], seq_2: &[u8], datatype: &DataType) -> Self {
        let mut counts = Self::default();
        let is_dna = !matches!(datatype, DataType::Aa);
        seq_1
            .iter()
            .zip(seq_2)
            .filter(|(a, b)| **a != INVALID && **b != INVALID)
            .for_each(|(&a, &b)| {
                counts.sites += 1;
                if a == b {
                    return;
                }
                counts.diffs += 1;
                if is_dna {
                    match (a.min(b), a.max(b)) {
                        (A, G) => counts.ts_purine += 1,
                        (C, T) => counts.ts_pyrimidine += 1,
                        _ => counts.transversions += 1,
                    }
                }
            });
        counts
    }
}

/// DNA base frequencies (A, C, G, T) used by the TN93 model.
#[derive(Debug, Clone, PartialEq)]
pub struct BaseFreqs {
    pub a: f64,
    pub c: f64,
    pub g: f64,
    pub t: f64,
}

impl Default for BaseFreqs {
    /// Equal base frequencies.
    fn default() -> Self {
        Self {
            a: 0.25,
            c: 0.25,
            g: 0.25,
            t: 0.25,
        }
    }
}

impl BaseFreqs {
    /// Estimate base frequencies from encoded DNA sequences.
    /// Falls back to equal frequencies if there are no valid bases.
    pub fn from_encoded(seqs: &[Vec<u8>]) -> Self {
        let mut counts = [0usize; 4];
        seqs.iter()
            .flatten()
            .filter(|b| **b <= T)
            .for_each(|b| counts[*b as usize] += 1);
        let total: usize = counts.iter().sum();
        if total == 0 {
            return Self::default();
        }
        let total = total as f64;
        Self {
            a: counts[A as usize] as f64 / total,
            c: counts[C as usize] as f64 / total,
            g: counts[G as usize] as f64 / total,
            t: counts[T as usize] as f64 / total,
        }
    }

    /// Tamura-Nei (1993) distance given the proportion of
    /// purine transitions (`p1`), pyrimidine transitions (`p2`),
    /// and transversions (`q`).
    fn tamura_nei(&self, p1: f64, p2: f64, q: f64) -> Option<f64> {
        let ag = self.a * self.g;
        let ct = self.c * self.t;
        if ag <= 0.0 || ct <= 0.0 {
            return None;
        }
        let r = self.a + self.g;
        let y = self.c + self.t;
        let term_r = 2.0 * ag / r * neg_ln(1.0 - r / (2.0 * ag) * p1 - q / (2.0 * r))?;
        let term_y = 2.0 * ct / y * neg_ln(1.0 - y / (2.0 * ct) * p2 - q / (2.0 * y))?;
        let term_q = 2.0 * (r * y - ag * y / r - ct * r / y) * neg_ln(1.0 - q / (2.0 * r * y))?;
        Some(term_r + term_y + term_q)
    }
}

/// Summary statistics of pairwise distances.
/// Fields are `None` if there are no values.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct DistanceSummary {
    pub min: Option<f64>,
    pub max: Option<f64>,
    pub mean: Option<f64>,
    pub median: Option<f64>,
    /// Sample standard deviation. `None` if there are fewer than two values.
    pub sd: Option<f64>,
}

impl DistanceSummary {
    pub fn new(values: &[f64]) -> Self {
        if values.is_empty() {
            return Self::default();
        }
        let mut sorted = values.to_vec();
        sorted.sort_by(|a, b| a.total_cmp(b));
        let n = sorted.len();
        let mean = sorted.iter().sum::<f64>() / n as f64;
        let median = if n.is_multiple_of(2) {
            (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
        } else {
            sorted[n / 2]
        };
        let sd = if n > 1 {
            let ss: f64 = sorted.iter().map(|v| (v - mean).powi(2)).sum();
            Some((ss / (n - 1) as f64).sqrt())
        } else {
            None
        };
        Self {
            min: sorted.first().copied(),
            max: sorted.last().copied(),
            mean: Some(mean),
            median: Some(median),
            sd,
        }
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

    fn counts(seq_1: &str, seq_2: &str, datatype: &DataType) -> PairCounts {
        PairCounts::new(&encode(seq_1, datatype), &encode(seq_2, datatype), datatype)
    }

    fn distance(model: DistanceModel, seq_1: &str, seq_2: &str) -> Option<f64> {
        let counts = counts(seq_1, seq_2, &DataType::Dna);
        model.distance(&counts, &DataType::Dna, &BaseFreqs::default())
    }

    // 10 sites: 1 A<->G transition, 1 C<->T transition, 1 transversion
    const SEQ_1: &str = "AACGTACGTA";
    const SEQ_2: &str = "GATGTACGTC";

    #[test]
    fn test_encode_dna() {
        assert_eq!(
            encode("ACGTUacgtu", &DataType::Dna),
            vec![0, 1, 2, 3, 3, 0, 1, 2, 3, 3]
        );
        assert!(
            encode("-?NRYn.", &DataType::Dna)
                .iter()
                .all(|b| *b == INVALID)
        );
    }

    #[test]
    fn test_encode_aa() {
        assert_eq!(encode("ARv", &DataType::Aa), vec![0, 1, 19]);
        assert!(
            encode("-?XBZJ*", &DataType::Aa)
                .iter()
                .all(|b| *b == INVALID)
        );
    }

    #[test]
    fn test_pair_counts() {
        let res = counts(SEQ_1, SEQ_2, &DataType::Dna);
        let expected = PairCounts {
            sites: 10,
            diffs: 3,
            ts_purine: 1,
            ts_pyrimidine: 1,
            transversions: 1,
        };
        assert_eq!(res, expected);
    }

    #[test]
    fn test_pair_counts_pairwise_deletion() {
        // Only sites 0 and 4 are compared
        let res = counts("A-C?GT", "GANAgR", &DataType::Dna);
        assert_eq!(res.sites, 2);
        assert_eq!(res.diffs, 1);
        assert_eq!(res.ts_purine, 1);
    }

    #[test]
    fn test_uncorrected() {
        let d = distance(DistanceModel::Uncorrected, SEQ_1, SEQ_2).unwrap();
        assert_approx_eq!(d, 0.3);
    }

    #[test]
    fn test_jukes_cantor() {
        let d = distance(DistanceModel::JukesCantor, SEQ_1, SEQ_2).unwrap();
        assert_approx_eq!(d, -0.75 * (1.0f64 - 0.4).ln());
    }

    #[test]
    fn test_kimura_2p() {
        let d = distance(DistanceModel::Kimura2P, SEQ_1, SEQ_2).unwrap();
        let expected = -0.5 * (1.0f64 - 0.4 - 0.1).ln() - 0.25 * (1.0f64 - 0.2).ln();
        assert_approx_eq!(d, expected);
    }

    #[test]
    fn test_tamura_nei_equals_k2p_with_equal_freqs() {
        let tn93 = distance(DistanceModel::TamuraNei, SEQ_1, SEQ_2).unwrap();
        let k2p = distance(DistanceModel::Kimura2P, SEQ_1, SEQ_2).unwrap();
        assert_approx_eq!(tn93, k2p);
    }

    #[test]
    fn test_tamura_nei_unequal_freqs() {
        let seqs = vec![encode(SEQ_1, &DataType::Dna), encode(SEQ_2, &DataType::Dna)];
        let freqs = BaseFreqs::from_encoded(&seqs);
        assert_approx_eq!(freqs.a, 0.3);
        assert_approx_eq!(freqs.c, 0.2);
        assert_approx_eq!(freqs.g, 0.25);
        assert_approx_eq!(freqs.t, 0.25);
        let counts = PairCounts::new(&seqs[0], &seqs[1], &DataType::Dna);
        let d = DistanceModel::TamuraNei
            .distance(&counts, &DataType::Dna, &freqs)
            .unwrap();
        let (r, y, ag, ct) = (0.55, 0.45, 0.075, 0.05);
        let expected = -2.0 * ag / r * (1.0f64 - r / (2.0 * ag) * 0.1 - 0.1 / (2.0 * r)).ln()
            - 2.0 * ct / y * (1.0f64 - y / (2.0 * ct) * 0.1 - 0.1 / (2.0 * y)).ln()
            - 2.0 * (r * y - ag * y / r - ct * r / y) * (1.0f64 - 0.1 / (2.0 * r * y)).ln();
        assert_approx_eq!(d, expected);
    }

    #[test]
    fn test_identical_sequences() {
        for model in [
            DistanceModel::Uncorrected,
            DistanceModel::JukesCantor,
            DistanceModel::Kimura2P,
            DistanceModel::TamuraNei,
        ] {
            let d = distance(model, "ACGT", "acgt").unwrap();
            assert_eq!(d.to_string(), "0");
        }
    }

    #[test]
    fn test_saturated_distance() {
        // p = 0.75
        assert_eq!(distance(DistanceModel::JukesCantor, "ACGT", "CATT"), None);
        assert!(distance(DistanceModel::Uncorrected, "ACGT", "CATT").is_some());
    }

    #[test]
    fn test_no_compared_sites() {
        assert_eq!(distance(DistanceModel::Uncorrected, "AC--", "--GT"), None);
    }

    #[test]
    fn test_amino_acid_models() {
        let counts = counts("ARNDCQEGHI", "ARNDCQEGHL", &DataType::Aa);
        assert_eq!(counts.transversions, 0);
        let freqs = BaseFreqs::default();
        let p = DistanceModel::Uncorrected
            .distance(&counts, &DataType::Aa, &freqs)
            .unwrap();
        assert_approx_eq!(p, 0.1);
        let jc = DistanceModel::JukesCantor
            .distance(&counts, &DataType::Aa, &freqs)
            .unwrap();
        assert_approx_eq!(jc, -0.95 * (1.0f64 - 0.1 / 0.95).ln());
        let poisson = DistanceModel::Poisson
            .distance(&counts, &DataType::Aa, &freqs)
            .unwrap();
        assert_approx_eq!(poisson, -(0.9f64).ln());
    }

    #[test]
    fn test_model_validity() {
        assert!(DistanceModel::Kimura2P.is_valid_for(&DataType::Dna));
        assert!(!DistanceModel::Kimura2P.is_valid_for(&DataType::Aa));
        assert!(!DistanceModel::TamuraNei.is_valid_for(&DataType::Aa));
        assert!(DistanceModel::Poisson.is_valid_for(&DataType::Aa));
        assert!(!DistanceModel::Poisson.is_valid_for(&DataType::Dna));
        assert!(DistanceModel::JukesCantor.is_valid_for(&DataType::Aa));
        assert!(!DistanceModel::Uncorrected.is_valid_for(&DataType::Ignore));
    }

    #[test]
    fn test_model_from_str() {
        assert_eq!("k2p".parse::<DistanceModel>(), Ok(DistanceModel::Kimura2P));
        assert!("logdet".parse::<DistanceModel>().is_err());
    }

    #[test]
    fn test_distance_summary() {
        let summary = DistanceSummary::new(&[0.4, 0.1, 0.3, 0.2]);
        assert_eq!(summary.min, Some(0.1));
        assert_eq!(summary.max, Some(0.4));
        assert_approx_eq!(summary.mean.unwrap(), 0.25);
        assert_approx_eq!(summary.median.unwrap(), 0.25);
        assert_approx_eq!(summary.sd.unwrap(), (0.05f64 / 3.0).sqrt());
        assert_eq!(DistanceSummary::new(&[0.5]).sd, None);
        assert_eq!(DistanceSummary::new(&[]), DistanceSummary::default());
    }
}

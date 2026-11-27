//! Technology- and library-aware pass, warn, and fail checks of read quality.
//!
//! Thresholds are shared by every viewer of `read-report.json`, so a sample
//! gets the same status in the terminal and in the GUI.
use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

use super::cohort::PairedReadSummary;
use super::report::SampleReadReport;

/// Sequencing regime or library profile used to evaluate quality checks.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum QcProfile {
    /// Infer profile per sample from header platform and read lengths.
    #[default]
    Auto,
    /// Short-read shotgun / WGS libraries.
    ShortWgs,
    /// Short-read amplicon, target capture (UCE/exon), or RNA-seq.
    Amplicon,
    /// Oxford Nanopore long reads.
    LongOnt,
    /// PacBio HiFi / CCS long reads.
    LongHifi,
}

impl QcProfile {
    pub fn label(self) -> &'static str {
        match self {
            Self::Auto => "Auto",
            Self::ShortWgs => "Short (WGS)",
            Self::Amplicon => "Amplicon / Capture",
            Self::LongOnt => "Long (ONT)",
            Self::LongHifi => "Long (HiFi)",
        }
    }
}

impl fmt::Display for QcProfile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::Auto => "auto",
            Self::ShortWgs => "short-wgs",
            Self::Amplicon => "amplicon",
            Self::LongOnt => "long-ont",
            Self::LongHifi => "long-hifi",
        };
        write!(f, "{name}")
    }
}

impl FromStr for QcProfile {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "short-wgs" | "short" | "wgs" => Ok(Self::ShortWgs),
            "amplicon" | "capture" | "rna" => Ok(Self::Amplicon),
            "long-ont" | "ont" | "nanopore" => Ok(Self::LongOnt),
            "long-hifi" | "hifi" | "pacbio" => Ok(Self::LongHifi),
            _ => Err(format!("{} is not a valid QC profile", s)),
        }
    }
}

/// Outcome of a quality check, from best to worst.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum QcStatus {
    Pass,
    Warn,
    Fail,
}

impl QcStatus {
    pub fn label(&self) -> &'static str {
        match self {
            QcStatus::Pass => "Pass",
            QcStatus::Warn => "Warn",
            QcStatus::Fail => "Fail",
        }
    }
}

/// Thresholds of a quality check and what it measures.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QcRule {
    pub description: &'static str,
    /// Values past this warn.
    pub warn: f64,
    /// Values past this fail.
    pub fail: f64,
    pub higher_is_better: bool,
    /// Suffix of values, such as "%".
    pub unit: &'static str,
}

impl QcRule {
    pub fn status_of(&self, value: f64) -> QcStatus {
        if self.higher_is_better {
            if value < self.fail {
                QcStatus::Fail
            } else if value < self.warn {
                QcStatus::Warn
            } else {
                QcStatus::Pass
            }
        } else if value >= self.fail {
            QcStatus::Fail
        } else if value >= self.warn {
            QcStatus::Warn
        } else {
            QcStatus::Pass
        }
    }

    /// The range of values for each status, from pass to fail.
    pub fn ranges(&self) -> [(QcStatus, String); 3] {
        let warn = self.show(self.warn);
        let fail = self.show(self.fail);
        if self.higher_is_better {
            [
                (QcStatus::Pass, format!("{warn} or higher")),
                (QcStatus::Warn, format!("{fail} to below {warn}")),
                (QcStatus::Fail, format!("Below {fail}")),
            ]
        } else {
            [
                (QcStatus::Pass, format!("Below {warn}")),
                (QcStatus::Warn, format!("{warn} to below {fail}")),
                (QcStatus::Fail, format!("{fail} or higher")),
            ]
        }
    }

    /// `value` with the unit, without decimals for whole numbers.
    fn show(&self, value: f64) -> String {
        format!("{value}{}", self.unit)
    }
}

pub const Q30_RULE: QcRule = QcRule {
    description: "Share of bases with a quality score of 30 or more, \
        a 1 in 1,000 chance of a wrong base call.",
    warn: 80.0,
    fail: 70.0,
    higher_is_better: true,
    unit: "%",
};

pub const ONT_Q15_RULE: QcRule = QcRule {
    description: "Share of long-read bases with a quality score of 15 or more.",
    warn: 80.0,
    fail: 60.0,
    higher_is_better: true,
    unit: "%",
};

pub const HIFI_Q20_RULE: QcRule = QcRule {
    description: "Share of HiFi bases with a quality score of 20 or more.",
    warn: 90.0,
    fail: 80.0,
    higher_is_better: true,
    unit: "%",
};

pub const MEAN_QUALITY_RULE: QcRule = QcRule {
    description: "Mean quality score of every base.",
    warn: 28.0,
    fail: 20.0,
    higher_is_better: true,
    unit: "",
};

pub const ONT_MEAN_QUALITY_RULE: QcRule = QcRule {
    description: "Mean quality score of every base for Oxford Nanopore reads.",
    warn: 18.0,
    fail: 12.0,
    higher_is_better: true,
    unit: "",
};

pub const N_CONTENT_RULE: QcRule = QcRule {
    description: "Share of bases called as N, where the sequencer could not tell the base.",
    warn: 5.0,
    fail: 20.0,
    higher_is_better: false,
    unit: "%",
};

pub const LONG_N_CONTENT_RULE: QcRule = QcRule {
    description: "Share of long-read bases called as N.",
    warn: 1.0,
    fail: 5.0,
    higher_is_better: false,
    unit: "%",
};

pub const DUPLICATION_RULE: QcRule = QcRule {
    description: "Share of reads that repeat another read. High duplication can mean \
        low library complexity or PCR over-amplification.",
    warn: 20.0,
    fail: 50.0,
    higher_is_better: false,
    unit: "%",
};

pub const AMPLICON_DUPLICATION_RULE: QcRule = QcRule {
    description: "Share of reads that repeat another read in an amplicon or target-capture library.",
    warn: 85.0,
    fail: 98.0,
    higher_is_better: false,
    unit: "%",
};

pub const ADAPTERS_RULE: QcRule = QcRule {
    description: "Share of reads with a known adapter sequence. \
        High values suggest the reads need trimming.",
    warn: 5.0,
    fail: 10.0,
    higher_is_better: false,
    unit: "%",
};

pub const POLY_G_RULE: QcRule = QcRule {
    description: "Share of reads ending in a run of G, which two-color sequencers \
        call when the signal is lost.",
    warn: 5.0,
    fail: 10.0,
    higher_is_better: false,
    unit: "%",
};

pub const OVERREPRESENTED_RULE: QcRule = QcRule {
    description: "Share of reads from the most common read start. Overrepresented \
        sequences can be adapters, primers, or contamination.",
    warn: 0.1,
    fail: 1.0,
    higher_is_better: false,
    unit: "%",
};

pub const AMPLICON_OVERREPRESENTED_RULE: QcRule = QcRule {
    description: "Share of reads from the most common read start in an amplicon or capture library.",
    warn: 10.0,
    fail: 25.0,
    higher_is_better: false,
    unit: "%",
};

pub const LOW_COMPLEXITY_RULE: QcRule = QcRule {
    description: "Share of reads with low dinucleotide complexity, such as repeats or homopolymers.",
    warn: 10.0,
    fail: 25.0,
    higher_is_better: false,
    unit: "%",
};

pub const CHASTITY_FAILED_RULE: QcRule = QcRule {
    description: "Share of reads flagged as failing the Illumina chastity filter (:Y:).",
    warn: 1.0,
    fail: 5.0,
    higher_is_better: false,
    unit: "%",
};

pub const PAIR_SYNC_RULE: QcRule = QcRule {
    description: "Absolute read-count difference between R1 and R2 mates.",
    warn: 1.0,
    fail: 1.0,
    higher_is_better: false,
    unit: "",
};

pub const PAIR_Q30_DROP_RULE: QcRule = QcRule {
    description: "Drop in Q30 base share from R1 to R2. Large drops indicate read-2 phasing or chemistry issues.",
    warn: 10.0,
    fail: 20.0,
    higher_is_better: false,
    unit: "%",
};

/// One quality check of a sample or pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QcCheck {
    pub name: &'static str,
    /// The measured value, in the unit of the rule.
    pub value: f64,
    pub status: QcStatus,
    pub rule: QcRule,
}

impl QcCheck {
    pub fn new(name: &'static str, value: f64, rule: QcRule) -> Self {
        Self {
            name,
            value,
            status: rule.status_of(value),
            rule,
        }
    }

    /// The value as shown with the check, such as "92.4%".
    pub fn value_text(&self) -> String {
        format!("{:.1}{}", self.value, self.rule.unit)
    }
}

/// Resolve [`QcProfile::Auto`] into a concrete profile for `sample`.
pub fn effective_profile(sample: &SampleReadReport, requested: QcProfile) -> QcProfile {
    if requested != QcProfile::Auto {
        return requested;
    }
    let m = &sample.metrics;
    match m.platform.as_deref() {
        Some("Nanopore") => QcProfile::LongOnt,
        Some("PacBio") => {
            if m.q20_pct >= 75.0 || m.mean_qscore >= 22.0 {
                QcProfile::LongHifi
            } else {
                QcProfile::LongOnt
            }
        }
        Some("Illumina") => QcProfile::ShortWgs,
        _ => {
            if m.n50_read_length >= 1_000 {
                if m.q20_pct >= 75.0 && m.mean_qscore >= 22.0 {
                    QcProfile::LongHifi
                } else {
                    QcProfile::LongOnt
                }
            } else {
                QcProfile::ShortWgs
            }
        }
    }
}

/// Quality checks of `sample` using [`QcProfile::Auto`].
pub fn qc_checks(sample: &SampleReadReport) -> Vec<QcCheck> {
    qc_checks_with_profile(sample, QcProfile::Auto)
}

/// Quality checks of `sample` under `profile`. Checks that need Complete mode
/// statistics are left out when those statistics are missing.
pub fn qc_checks_with_profile(sample: &SampleReadReport, profile: QcProfile) -> Vec<QcCheck> {
    let resolved = effective_profile(sample, profile);
    let m = &sample.metrics;
    let mut checks = match resolved {
        QcProfile::LongOnt => vec![
            QcCheck::new("Q15 bases", m.q15_pct, ONT_Q15_RULE),
            QcCheck::new("Mean quality", m.mean_qscore, ONT_MEAN_QUALITY_RULE),
            QcCheck::new("N content", 100.0 * m.n_content, LONG_N_CONTENT_RULE),
        ],
        QcProfile::LongHifi => vec![
            QcCheck::new("Q20 bases", m.q20_pct, HIFI_Q20_RULE),
            QcCheck::new("Mean quality", m.mean_qscore, MEAN_QUALITY_RULE),
            QcCheck::new("N content", 100.0 * m.n_content, LONG_N_CONTENT_RULE),
        ],
        QcProfile::ShortWgs | QcProfile::Amplicon | QcProfile::Auto => vec![
            QcCheck::new("Q30 bases", m.q30_pct, Q30_RULE),
            QcCheck::new("Mean quality", m.mean_qscore, MEAN_QUALITY_RULE),
            QcCheck::new("N content", 100.0 * m.n_content, N_CONTENT_RULE),
        ],
    };

    if matches!(resolved, QcProfile::ShortWgs | QcProfile::Amplicon) {
        let dup_rule = if resolved == QcProfile::Amplicon {
            AMPLICON_DUPLICATION_RULE
        } else {
            DUPLICATION_RULE
        };
        let poly_g_value = if sample.homopolymers.is_some() && m.poly_g_reads_pct.is_some() {
            Some(m.poly_g_tail_reads_pct)
        } else {
            m.poly_g_reads_pct
        };
        for (name, value, rule) in [
            ("Duplication", m.duplicate_pct_est, dup_rule),
            ("Adapters", m.adapter_reads_pct, ADAPTERS_RULE),
            ("Poly-G", poly_g_value, POLY_G_RULE),
        ] {
            if let Some(value) = value {
                checks.push(QcCheck::new(name, value, rule));
            }
        }
        if m.overrepresented_count.is_some() {
            let top = sample.overrepresented.first().map_or(0.0, |o| o.percent);
            let rule = if resolved == QcProfile::Amplicon {
                AMPLICON_OVERREPRESENTED_RULE
            } else {
                OVERREPRESENTED_RULE
            };
            checks.push(QcCheck::new("Overrepresented", top, rule));
        }
    }

    // Only emit complexity check on reports that computed it.
    if sample.homopolymers.is_some() {
        checks.push(QcCheck::new(
            "Low complexity",
            m.low_complexity_reads_pct,
            LOW_COMPLEXITY_RULE,
        ));
    }
    if m.chastity_failed_reads > 0 {
        checks.push(QcCheck::new(
            "Chastity failed",
            m.chastity_failed_pct,
            CHASTITY_FAILED_RULE,
        ));
    }
    checks
}

/// Quality checks for a paired-end R1/R2 sample pair.
pub fn pair_qc_checks(pair: &PairedReadSummary) -> Vec<QcCheck> {
    vec![
        QcCheck::new(
            "Mate count diff",
            pair.read_count_diff.unsigned_abs() as f64,
            PAIR_SYNC_RULE,
        ),
        QcCheck::new(
            "R2 Q30 drop",
            pair.q30_drop_pct.max(0.0),
            PAIR_Q30_DROP_RULE,
        ),
    ]
}

/// The worst status among `checks`, or pass when there are none.
pub fn overall_status(checks: &[QcCheck]) -> QcStatus {
    checks
        .iter()
        .map(|c| c.status)
        .max()
        .unwrap_or(QcStatus::Pass)
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::stats::fastq::report::ReadReport;

    #[test]
    fn test_status_higher_is_better() {
        assert_eq!(Q30_RULE.status_of(80.0), QcStatus::Pass);
        assert_eq!(Q30_RULE.status_of(79.9), QcStatus::Warn);
        assert_eq!(Q30_RULE.status_of(70.0), QcStatus::Warn);
        assert_eq!(Q30_RULE.status_of(69.9), QcStatus::Fail);
    }

    #[test]
    fn test_status_lower_is_better() {
        assert_eq!(ADAPTERS_RULE.status_of(4.9), QcStatus::Pass);
        assert_eq!(ADAPTERS_RULE.status_of(5.0), QcStatus::Warn);
        assert_eq!(ADAPTERS_RULE.status_of(10.0), QcStatus::Fail);
    }

    #[test]
    fn test_ranges() {
        let ranges = Q30_RULE.ranges();
        assert_eq!(ranges[0], (QcStatus::Pass, String::from("80% or higher")));
        assert_eq!(ranges[1].1, "70% to below 80%");
        assert_eq!(ranges[2].1, "Below 70%");
        let ranges = OVERREPRESENTED_RULE.ranges();
        assert_eq!(ranges[0].1, "Below 0.1%");
        assert_eq!(ranges[2].1, "1% or higher");
        assert_eq!(MEAN_QUALITY_RULE.ranges()[0].1, "28 or higher");
    }

    #[test]
    fn test_checks() {
        let report = ReadReport::from_path(std::path::Path::new(
            "../../tests/files/plot/read-report.json",
        ))
        .expect("Failed reading fixture");
        let sample = &report.samples[2];
        let checks = qc_checks(sample);
        let names: Vec<&str> = checks.iter().map(|c| c.name).collect();
        assert_eq!(names[..3], ["Q30 bases", "Mean quality", "N content"]);
        assert!(names.contains(&"Duplication"));
        let mut quick = sample.clone();
        quick.metrics.duplicate_pct_est = None;
        quick.metrics.adapter_reads_pct = None;
        quick.metrics.poly_g_reads_pct = None;
        quick.metrics.overrepresented_count = None;
        assert_eq!(qc_checks(&quick).len(), 3);
    }

    #[test]
    fn test_profile_resolution_and_ont_thresholds() {
        let mut ont = SampleReadReport::default();
        ont.metrics.platform = Some("Nanopore".to_string());
        ont.metrics.n50_read_length = 8_000;
        ont.metrics.q15_pct = 85.0;
        ont.metrics.mean_qscore = 19.0;
        assert_eq!(effective_profile(&ont, QcProfile::Auto), QcProfile::LongOnt);
        let checks = qc_checks_with_profile(&ont, QcProfile::Auto);
        assert_eq!(checks[0].name, "Q15 bases");
        assert_eq!(overall_status(&checks), QcStatus::Pass);

        let mut amp = SampleReadReport::default();
        amp.metrics.q30_pct = 90.0;
        amp.metrics.mean_qscore = 35.0;
        amp.metrics.duplicate_pct_est = Some(75.0);
        assert_eq!(
            overall_status(&qc_checks_with_profile(&amp, QcProfile::ShortWgs)),
            QcStatus::Fail
        );
        assert_eq!(
            overall_status(&qc_checks_with_profile(&amp, QcProfile::Amplicon)),
            QcStatus::Pass
        );
    }

    #[test]
    fn test_pair_checks() {
        let mut pair = PairedReadSummary {
            in_sync: true,
            read_count_diff: 0,
            q30_drop_pct: 4.0,
            ..Default::default()
        };
        assert_eq!(overall_status(&pair_qc_checks(&pair)), QcStatus::Pass);
        pair.in_sync = false;
        pair.read_count_diff = 3;
        assert_eq!(overall_status(&pair_qc_checks(&pair)), QcStatus::Fail);
    }

    #[test]
    fn test_overall_status() {
        assert_eq!(overall_status(&[]), QcStatus::Pass);
        let checks = [
            QcCheck::new("Q30 bases", 75.0, Q30_RULE),
            QcCheck::new("Adapters", 1.0, ADAPTERS_RULE),
        ];
        assert_eq!(overall_status(&checks), QcStatus::Warn);
        assert_eq!(checks[0].value_text(), "75.0%");
    }
}

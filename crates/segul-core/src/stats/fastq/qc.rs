//! Pass, warn, and fail checks of read quality, following FastQC
//! conventions.
//!
//! The thresholds are shared by every viewer of `read-report.json`, so a
//! sample gets the same status in the terminal and in the GUI.

use super::report::SampleReadReport;

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

    /// [value] with the unit, without decimals for whole numbers.
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

pub const MEAN_QUALITY_RULE: QcRule = QcRule {
    description: "Mean quality score of every base, following the FastQC quality zones.",
    warn: 28.0,
    fail: 20.0,
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

pub const DUPLICATION_RULE: QcRule = QcRule {
    description: "Share of reads that repeat another read. High duplication can mean \
        low library complexity or PCR over-amplification.",
    warn: 20.0,
    fail: 50.0,
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

/// One quality check of a sample.
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

/// Quality checks of [sample]. Checks that need Complete mode statistics
/// are left out when those statistics are missing.
pub fn qc_checks(sample: &SampleReadReport) -> Vec<QcCheck> {
    let m = &sample.metrics;
    let mut checks = vec![
        QcCheck::new("Q30 bases", m.q30_pct, Q30_RULE),
        QcCheck::new("Mean quality", m.mean_qscore, MEAN_QUALITY_RULE),
        QcCheck::new("N content", 100.0 * m.n_content, N_CONTENT_RULE),
    ];
    let optional = [
        ("Duplication", m.duplicate_pct_est, DUPLICATION_RULE),
        ("Adapters", m.adapter_reads_pct, ADAPTERS_RULE),
        ("Poly-G", m.poly_g_reads_pct, POLY_G_RULE),
    ];
    for (name, value, rule) in optional {
        if let Some(value) = value {
            checks.push(QcCheck::new(name, value, rule));
        }
    }
    if m.overrepresented_count.is_some() {
        let top = sample.overrepresented.first().map_or(0.0, |o| o.percent);
        checks.push(QcCheck::new("Overrepresented", top, OVERREPRESENTED_RULE));
    }
    checks
}

/// The worst status among [checks], or pass when there are none.
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
        let report =
            ReadReport::from_path(std::path::Path::new("tests/files/plot/read-report.json"))
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

//! Summary reports as overview values, charts, and tables for the viewer.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::stats::align_report::{AlignReport, SiteTypeSummary};
use crate::stats::distance_report::{DistanceReport, PairStats};
use crate::stats::fastq::report::ReadReport;

/// Most bars in the pairwise distance histogram.
const DISTANCE_BINS: usize = 20;

/// A report ready to show: overview values and charts, one table per
/// kind of item, such as loci and taxa, and any square matrices.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct SummaryData {
    /// What the report summarizes, such as "Alignment summary".
    pub(super) kind: &'static str,
    /// Short counts shown beside the file name.
    pub(super) headline: String,
    pub(super) overview: Vec<(String, String)>,
    pub(super) charts: Vec<Chart>,
    pub(super) tables: Vec<SummaryTable>,
    /// Shown after the tables.
    pub(super) matrices: Vec<SummaryMatrix>,
}

/// Values for every pair of items, such as distances between samples.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct SummaryMatrix {
    /// Name of the matrix, such as "Matrix".
    pub(super) name: &'static str,
    /// Plural name of the items, such as "samples".
    pub(super) item: &'static str,
    pub(super) names: Vec<String>,
    /// `values[i][j]` is the value of items `i` and `j`.
    pub(super) values: Vec<Vec<Option<f64>>>,
    pub(super) format: Format,
}

/// Labeled values drawn as bars on the overview.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Chart {
    pub(super) title: String,
    pub(super) labels: Vec<String>,
    pub(super) values: Vec<f64>,
    pub(super) format: Format,
}

/// Items such as loci, with one value per item in each column.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct SummaryTable {
    /// Plural name of the items, such as "Loci".
    pub(super) name: &'static str,
    /// Heading of the name column, such as "Locus".
    pub(super) item: &'static str,
    pub(super) names: Vec<String>,
    pub(super) columns: Vec<Column>,
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct Column {
    pub(super) label: &'static str,
    /// One value per item. Missing values are statistics a run did not
    /// compute.
    pub(super) values: Vec<Option<f64>>,
    pub(super) format: Format,
}

/// How values are written.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Format {
    /// Whole numbers with thousands separators.
    Count,
    /// One decimal.
    Decimal,
    /// Four decimals, for genetic distances.
    Distance,
}

impl Format {
    pub(super) fn show(self, value: f64) -> String {
        match self {
            Format::Count => thousands(value.round() as i128),
            Format::Decimal => format!("{:.1}", value),
            Format::Distance => format!("{:.4}", value),
        }
    }
}

/// [value] with commas between groups of three digits.
pub(super) fn thousands(value: i128) -> String {
    let digits = value.unsigned_abs().to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, ch) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    if value < 0 {
        out.insert(0, '-');
    }
    out
}

impl SummaryData {
    /// Read an alignment, read, or distance report, telling them apart by
    /// content.
    pub(super) fn from_path(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("Failed reading {}", path.display()))?;
        if let Ok(report) = AlignReport::from_json(&text) {
            return Ok(Self::from_alignment(&report));
        }
        if let Ok(report) = ReadReport::from_json(&text) {
            return Ok(Self::from_reads(&report));
        }
        if let Ok(report) = DistanceReport::from_json(&text) {
            return Ok(Self::from_distance(&report));
        }
        bail!(
            "{} is not a segul summary report. Use alignment-report.json from \
            `segul align summary`, read-report.json from \
            `segul read summary --mode complete`, or distance-report.json from \
            `segul align distance`.",
            path.display()
        )
    }

    pub(super) fn from_alignment(report: &AlignReport) -> Self {
        let s = &report.summary;
        let is_dna = report.datatype == "dna";
        let pct = |p: f64| format!("{:.1}%", 100.0 * p);
        let count = |n: usize| thousands(n as i128);
        let range = |min: usize, mean: f64, max: usize| {
            format!("{} to {}, mean {:.1}", count(min), count(max), mean)
        };
        let site_type = |t: &SiteTypeSummary, loci: &str| {
            format!(
                "{} ({}), {} {}",
                count(t.sites),
                pct(t.prop_sites),
                count(t.loci),
                loci
            )
        };
        let mut overview = vec![
            (
                String::from("Data type"),
                String::from(if is_dna { "DNA" } else { "Amino acid" }),
            ),
            (String::from("Taxa"), count(s.taxa)),
            (String::from("Loci"), count(s.loci)),
            (String::from("Sites"), count(s.sites)),
            (String::from("Characters"), count(s.chars)),
            (
                String::from("Missing data"),
                format!("{} ({})", count(s.missing), pct(s.prop_missing)),
            ),
        ];
        if let (Some(gc), Some(at)) = (s.gc, s.at) {
            overview.push((String::from("GC content"), pct(gc)));
            overview.push((String::from("AT content"), pct(at)));
        }
        overview.extend([
            (
                String::from("Locus length"),
                range(s.min_sites, s.mean_sites, s.max_sites),
            ),
            (
                String::from("Taxa per locus"),
                range(s.min_taxa, s.mean_taxa, s.max_taxa),
            ),
            (
                String::from("Conserved sites"),
                site_type(&s.conserved, "loci fully conserved"),
            ),
            (
                String::from("Variable sites"),
                site_type(&s.variable, "loci with any"),
            ),
            (
                String::from("Informative sites"),
                site_type(&s.informative, "loci with any"),
            ),
        ]);
        let charts = vec![
            Chart {
                title: String::from("Loci with at least this share of taxa"),
                labels: report
                    .completeness
                    .iter()
                    .map(|c| format!("{}%", c.percent))
                    .collect(),
                values: report.completeness.iter().map(|c| c.loci as f64).collect(),
                format: Format::Count,
            },
            Chart {
                title: String::from("Characters"),
                labels: report.chars.keys().map(|c| character_label(c)).collect(),
                values: report.chars.values().map(|&c| c as f64).collect(),
                format: Format::Count,
            },
        ];
        let loci = &report.loci;
        let taxa = &report.taxa;
        let counts = |values: Vec<usize>| values.into_iter().map(|v| Some(v as f64)).collect();
        let percents =
            |values: Vec<Option<f64>>| values.into_iter().map(|v| v.map(|p| 100.0 * p)).collect();
        let mut loci_columns = vec![
            Column {
                label: "Taxa",
                values: counts(loci.iter().map(|l| l.ntax).collect()),
                format: Format::Count,
            },
            Column {
                label: "Sites",
                values: counts(loci.iter().map(|l| l.sites).collect()),
                format: Format::Count,
            },
            Column {
                label: "Conserved",
                values: counts(loci.iter().map(|l| l.conserved).collect()),
                format: Format::Count,
            },
            Column {
                label: "Variable",
                values: counts(loci.iter().map(|l| l.variable).collect()),
                format: Format::Count,
            },
            Column {
                label: "Informative",
                values: counts(loci.iter().map(|l| l.informative).collect()),
                format: Format::Count,
            },
            Column {
                label: "Missing %",
                values: percents(loci.iter().map(|l| Some(l.prop_missing)).collect()),
                format: Format::Decimal,
            },
        ];
        let mut taxa_columns = vec![
            Column {
                label: "Loci",
                values: counts(taxa.iter().map(|t| t.loci).collect()),
                format: Format::Count,
            },
            Column {
                label: "Characters",
                values: counts(taxa.iter().map(|t| t.chars).collect()),
                format: Format::Count,
            },
            Column {
                label: "Missing",
                values: counts(taxa.iter().map(|t| t.missing).collect()),
                format: Format::Count,
            },
            Column {
                label: "Missing %",
                values: percents(taxa.iter().map(|t| Some(t.prop_missing)).collect()),
                format: Format::Decimal,
            },
        ];
        if is_dna {
            loci_columns.push(Column {
                label: "GC %",
                values: percents(loci.iter().map(|l| l.gc).collect()),
                format: Format::Decimal,
            });
            taxa_columns.push(Column {
                label: "GC %",
                values: percents(taxa.iter().map(|t| t.gc).collect()),
                format: Format::Decimal,
            });
        }
        Self {
            kind: "Alignment summary",
            headline: format!(
                "{} taxa · {} loci · {}",
                count(s.taxa),
                count(s.loci),
                if is_dna { "DNA" } else { "Amino acid" }
            ),
            overview,
            charts,
            tables: vec![
                SummaryTable {
                    name: "Loci",
                    item: "Locus",
                    names: loci.iter().map(|l| l.name.clone()).collect(),
                    columns: loci_columns,
                },
                SummaryTable {
                    name: "Taxa",
                    item: "Taxon",
                    names: taxa.iter().map(|t| t.name.clone()).collect(),
                    columns: taxa_columns,
                },
            ],
            matrices: Vec::new(),
        }
    }

    pub(super) fn from_reads(report: &ReadReport) -> Self {
        let samples = &report.samples;
        let reads: u64 = samples.iter().map(|s| s.metrics.read_count).sum();
        let bases: u64 = samples.iter().map(|s| s.metrics.base_count).sum();
        let q30: u64 = samples.iter().map(|s| s.metrics.q30_bases).sum();
        let ratio = |part: u64, whole: u64| {
            if whole == 0 {
                0.0
            } else {
                part as f64 / whole as f64
            }
        };
        let overview = vec![
            (String::from("Samples"), thousands(samples.len() as i128)),
            (String::from("Reads"), thousands(reads as i128)),
            (String::from("Bases"), thousands(bases as i128)),
            (
                String::from("Mean read length"),
                format!("{:.1}", ratio(bases, reads)),
            ),
            (
                String::from("Q30 bases"),
                format!("{:.1}%", 100.0 * ratio(q30, bases)),
            ),
        ];
        // Reads of every sample by mean quality, summed.
        let bins = samples
            .iter()
            .map(|s| s.read_quality_distribution.len())
            .max()
            .unwrap_or(0);
        let mut quality = vec![0.0; bins];
        for sample in samples {
            for (q, &n) in sample.read_quality_distribution.iter().enumerate() {
                quality[q] += n as f64;
            }
        }
        // Leave out the empty low scores, which are most of the range.
        let first = quality.iter().position(|&n| n > 0.0).unwrap_or(0);
        let charts = vec![Chart {
            title: String::from("Reads by mean quality score, all samples"),
            labels: (first..bins).map(|q| q.to_string()).collect(),
            values: quality[first..].to_vec(),
            format: Format::Count,
        }];
        let metric = |label: &'static str, format: Format, value: &dyn Fn(usize) -> Option<f64>| {
            let values: Vec<Option<f64>> = (0..samples.len()).map(value).collect();
            (values.iter().any(Option::is_some)).then_some(Column {
                label,
                values,
                format,
            })
        };
        let m = |i: usize| &samples[i].metrics;
        let columns: Vec<Column> = [
            metric("Reads", Format::Count, &|i| Some(m(i).read_count as f64)),
            metric("Bases", Format::Count, &|i| Some(m(i).base_count as f64)),
            metric("Mean length", Format::Decimal, &|i| {
                Some(m(i).mean_read_length)
            }),
            metric("N50", Format::Count, &|i| Some(m(i).n50_read_length as f64)),
            metric("GC %", Format::Decimal, &|i| Some(100.0 * m(i).gc_content)),
            metric("N %", Format::Decimal, &|i| Some(100.0 * m(i).n_content)),
            metric("Q20 %", Format::Decimal, &|i| Some(m(i).q20_pct)),
            metric("Q30 %", Format::Decimal, &|i| Some(m(i).q30_pct)),
            metric("Mean Q", Format::Decimal, &|i| Some(m(i).mean_qscore)),
            metric("Low-Q reads", Format::Count, &|i| {
                Some(m(i).low_quality_reads as f64)
            }),
            metric("Duplication %", Format::Decimal, &|i| {
                m(i).duplicate_pct_est
            }),
            metric("Adapters %", Format::Decimal, &|i| m(i).adapter_reads_pct),
            metric("Poly-G %", Format::Decimal, &|i| m(i).poly_g_reads_pct),
        ]
        .into_iter()
        .flatten()
        .collect();
        Self {
            kind: "Read summary",
            headline: format!(
                "{} samples · {} reads",
                thousands(samples.len() as i128),
                thousands(reads as i128)
            ),
            overview,
            charts,
            tables: vec![SummaryTable {
                name: "Samples",
                item: "Sample",
                names: samples.iter().map(|s| s.sample.clone()).collect(),
                columns,
            }],
            matrices: Vec::new(),
        }
    }

    pub(super) fn from_distance(report: &DistanceReport) -> Self {
        let s = &report.summary;
        let is_dna = report.datatype == "dna";
        let count = |n: usize| thousands(n as i128);
        let distance = |d: Option<f64>| d.map_or_else(|| String::from("-"), |d| format!("{d:.4}"));
        let mut overview = vec![
            (
                String::from("Data type"),
                String::from(if is_dna { "DNA" } else { "Amino acid" }),
            ),
            (String::from("Model"), report.model.clone()),
            (String::from("Loci"), count(s.loci)),
            (String::from("Samples"), count(s.samples)),
            (String::from("Pairs"), count(s.pairs)),
            (String::from("Undefined pairs"), count(s.undefined_pairs)),
            (
                String::from("Distance"),
                format!(
                    "{} to {}, mean {} (SD {})",
                    distance(s.min),
                    distance(s.max),
                    distance(s.mean),
                    distance(s.sd)
                ),
            ),
            (
                String::from("Median locus mean"),
                distance(s.median_locus_mean),
            ),
            (
                String::from("Compared sites"),
                s.mean_compared_sites
                    .map_or_else(|| String::from("-"), |n| format!("{n:.1} per pair")),
            ),
        ];
        if let Some(groups) = s.groups {
            overview.push((String::from("Groups"), count(groups)));
            overview.push((String::from("Unassigned samples"), count(s.unassigned)));
        }
        let (labels, values) = distance_bars(report);
        let charts = vec![Chart {
            title: String::from("Pairwise distances, all loci"),
            labels,
            values,
            format: Format::Count,
        }];
        let counts = |values: Vec<usize>| values.into_iter().map(|v| Some(v as f64)).collect();
        let column = |label: &'static str, values: Vec<Option<f64>>| Column {
            label,
            values,
            format: Format::Distance,
        };
        let loci = &report.loci;
        let loci_columns = vec![
            Column {
                label: "Samples",
                values: counts(loci.iter().map(|l| l.ntax).collect()),
                format: Format::Count,
            },
            Column {
                label: "Sites",
                values: counts(loci.iter().map(|l| l.nchar).collect()),
                format: Format::Count,
            },
            Column {
                label: "Pairs",
                values: counts(loci.iter().map(|l| l.pairs).collect()),
                format: Format::Count,
            },
            Column {
                label: "Undefined",
                values: counts(loci.iter().map(|l| l.undefined_pairs).collect()),
                format: Format::Count,
            },
            column("Mean", loci.iter().map(|l| l.mean).collect()),
            column("Median", loci.iter().map(|l| l.median).collect()),
            column("SD", loci.iter().map(|l| l.sd).collect()),
            column("Min", loci.iter().map(|l| l.min).collect()),
            column("Max", loci.iter().map(|l| l.max).collect()),
            Column {
                label: "Compared sites",
                values: loci.iter().map(|l| l.mean_compared_sites).collect(),
                format: Format::Decimal,
            },
        ];
        let samples = &report.samples;
        let mut sample_columns = vec![
            Column {
                label: "Loci",
                values: counts(samples.iter().map(|s| s.loci).collect()),
                format: Format::Count,
            },
            Column {
                label: "Pairs",
                values: counts(samples.iter().map(|s| s.stats.pairs).collect()),
                format: Format::Count,
            },
            column("Mean", samples.iter().map(|s| s.stats.mean).collect()),
            column("SD", samples.iter().map(|s| s.stats.sd).collect()),
            column("Min", samples.iter().map(|s| s.stats.min).collect()),
            column("Max", samples.iter().map(|s| s.stats.max).collect()),
        ];
        if report.matrix.is_some() {
            sample_columns.push(column(
                "Nearest",
                samples.iter().map(|s| s.nearest_distance).collect(),
            ));
        }
        if report.matrix.is_some() && report.groups.is_some() {
            sample_columns.push(column(
                "Max within",
                samples.iter().map(|s| s.max_within).collect(),
            ));
            sample_columns.push(column(
                "Min between",
                samples.iter().map(|s| s.min_between).collect(),
            ));
        }
        let mut tables = vec![
            SummaryTable {
                name: "Loci",
                item: "Locus",
                names: loci.iter().map(|l| l.name.clone()).collect(),
                columns: loci_columns,
            },
            SummaryTable {
                name: "Samples",
                item: "Sample",
                names: samples.iter().map(|s| s.name.clone()).collect(),
                columns: sample_columns,
            },
        ];
        let mut matrices = Vec::new();
        if let Some(matrix) = &report.matrix {
            matrices.push(SummaryMatrix {
                name: "Matrix",
                item: "samples",
                names: matrix.names.clone(),
                values: matrix.mean.clone(),
                format: Format::Distance,
            });
        }
        if let Some(groups) = &report.groups {
            let names: Vec<String> = groups.groups.iter().map(|g| g.name.clone()).collect();
            let mut means = vec![vec![None; names.len()]; names.len()];
            let position = |name: &str| names.iter().position(|n| n == name);
            for pair in &groups.pairs {
                if let (Some(i), Some(j)) = (position(&pair.group_1), position(&pair.group_2)) {
                    means[i][j] = pair.stats.mean;
                    means[j][i] = pair.stats.mean;
                }
            }
            // Lowest mean distance to any other group.
            let nearest: Vec<Option<f64>> = (0..names.len())
                .map(|i| {
                    (0..names.len())
                        .filter(|&j| j != i)
                        .filter_map(|j| means[i][j])
                        .reduce(f64::min)
                })
                .collect();
            let within = |value: fn(&PairStats) -> Option<f64>| {
                groups.groups.iter().map(|g| value(&g.within)).collect()
            };
            tables.push(SummaryTable {
                name: "Groups",
                item: "Group",
                names: names.clone(),
                columns: vec![
                    Column {
                        label: "Samples",
                        values: counts(groups.groups.iter().map(|g| g.samples).collect()),
                        format: Format::Count,
                    },
                    Column {
                        label: "Pairs",
                        values: counts(groups.groups.iter().map(|g| g.within.pairs).collect()),
                        format: Format::Count,
                    },
                    column("Mean within", within(|w| w.mean)),
                    column("SD within", within(|w| w.sd)),
                    column("Max within", within(|w| w.max)),
                    column("Nearest group", nearest),
                ],
            });
            matrices.push(SummaryMatrix {
                name: "Group matrix",
                item: "groups",
                names,
                values: means,
                format: Format::Distance,
            });
        }
        Self {
            kind: "Distance summary",
            headline: format!(
                "{} samples · {} loci · {}",
                count(s.samples),
                count(s.loci),
                report.model
            ),
            overview,
            charts,
            tables,
            matrices,
        }
    }
}

/// Labels and counts of the pairwise distance histogram in at most
/// [`DISTANCE_BINS`] bars of a round width.
fn distance_bars(report: &DistanceReport) -> (Vec<String>, Vec<f64>) {
    let histogram = &report.histogram;
    let counts = &histogram.counts;
    if counts.is_empty() || histogram.bin_width <= 0.0 {
        return (Vec::new(), Vec::new());
    }
    let width = histogram.bin_width;
    let step = round_step(counts.len() as f64 * width / DISTANCE_BINS as f64).max(width);
    let per_bar = (step / width).round().max(1.0) as usize;
    let decimals = (-step.log10().floor()).max(0.0) as usize;
    counts
        .chunks(per_bar)
        .enumerate()
        .map(|(i, bins)| {
            let from = i as f64 * per_bar as f64 * width;
            let to = from + per_bar as f64 * width;
            (
                format!("{from:.decimals$}-{to:.decimals$}"),
                bins.iter().sum::<usize>() as f64,
            )
        })
        .unzip()
}

/// A width of 1, 2, or 5 times a power of ten, at least [raw].
pub(super) fn round_step(raw: f64) -> f64 {
    if raw <= 0.0 {
        return 1.0;
    }
    let scale = 10_f64.powf(raw.log10().floor());
    [1.0, 2.0, 5.0, 10.0]
        .into_iter()
        .map(|s| s * scale)
        .find(|&s| s >= raw)
        .unwrap_or(10.0 * scale)
}

impl SummaryData {}

/// Chart label of an alignment character, naming gaps (`-`) and missing
/// data (`?`).
fn character_label(character: &str) -> String {
    match character {
        "-" => String::from("- (gap)"),
        "?" => String::from("? (missing)"),
        _ => character.to_string(),
    }
}

#[cfg(test)]
mod test {
    use super::*;

    const ALIGN: &str = "crates/segul-core/tests/files/plot/alignment-report.json";
    const READS: &str = "crates/segul-core/tests/files/plot/read-report.json";
    const DISTANCE: &str = "crates/segul-core/tests/files/plot/distance-report.json";

    #[test]
    fn test_character_label() {
        assert_eq!(character_label("-"), "- (gap)");
        assert_eq!(character_label("?"), "? (missing)");
        assert_eq!(character_label("A"), "A");
    }

    #[test]
    fn test_alignment_report() {
        let data = SummaryData::from_path(Path::new(ALIGN)).expect("Failed reading fixture");
        assert_eq!(data.kind, "Alignment summary");
        assert_eq!(data.headline, "3 taxa · 4 loci · DNA");
        let loci = &data.tables[0];
        assert_eq!(loci.names, vec!["gene_1", "gene_2", "gene_3", "gene_4"]);
        assert_eq!(loci.columns[1].label, "Sites");
        assert_eq!(loci.columns.last().unwrap().label, "GC %");
        assert_eq!(data.tables[1].names.len(), 3);
        assert!(data.overview.iter().any(|(k, _)| k == "GC content"));
    }

    #[test]
    fn test_read_report() {
        let data = SummaryData::from_path(Path::new(READS)).expect("Failed reading fixture");
        assert_eq!(data.kind, "Read summary");
        assert_eq!(data.tables.len(), 1);
        let samples = &data.tables[0];
        assert!(!samples.names.is_empty());
        assert_eq!(samples.columns[0].label, "Reads");
        assert!(samples.columns.iter().any(|c| c.label == "Q30 %"));
        assert!(!data.charts[0].values.is_empty());
    }

    #[test]
    fn test_distance_report() {
        let data = SummaryData::from_path(Path::new(DISTANCE)).expect("Failed reading fixture");
        assert_eq!(data.kind, "Distance summary");
        assert_eq!(data.headline, "6 samples · 3 loci · K2P");
        let names: Vec<&str> = data.tables.iter().map(|t| t.name).collect();
        assert_eq!(names, ["Loci", "Samples", "Groups"]);
        let samples = &data.tables[1];
        assert_eq!(samples.names[0], "Mus_musculus_A");
        assert_eq!(samples.columns.last().unwrap().label, "Min between");
        let groups = &data.tables[2];
        assert_eq!(groups.names.len(), 3);
        assert_eq!(data.matrices.len(), 2);
        let matrix = &data.matrices[0];
        assert_eq!(matrix.names.len(), 6);
        assert_eq!(matrix.values[0][0], Some(0.0));
        let means = &data.matrices[1];
        assert_eq!(means.values[0][1], means.values[1][0]);
        assert!(means.values[0][0].is_some());
        let chart = &data.charts[0];
        assert!(chart.values.len() <= DISTANCE_BINS);
        assert_eq!(chart.values.iter().sum::<f64>(), 31.0);
        assert_eq!(chart.labels[0], "0.00-0.05");
        assert!(data.overview.iter().any(|(k, v)| k == "Groups" && v == "3"));
    }

    #[test]
    fn test_other_json() {
        let dir = tempfile::TempDir::with_prefix("view").unwrap();
        let path = dir.path().join("other.json");
        fs::write(&path, "{\"a\": 1}").unwrap();
        let err = SummaryData::from_path(&path).unwrap_err().to_string();
        assert!(err.contains("not a segul summary report"));
    }

    #[test]
    fn test_thousands() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1_234_567), "1,234,567");
        assert_eq!(thousands(-1_000), "-1,000");
        assert_eq!(Format::Decimal.show(2.345), "2.3");
        assert_eq!(Format::Distance.show(0.123456), "0.1235");
        assert_eq!(round_step(0.04), 0.05);
    }
}

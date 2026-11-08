//! MAF summary plots.
//!
//! Species charts read `<prefix>_per_species.csv` and file charts read
//! `<prefix>_per_file.csv`, both written by `segul maf summary`.
use std::path::{Path, PathBuf};

use anyhow::{Result, bail};
use indexmap::IndexSet;

use super::PlotPreset;
use super::batch::PlotWriter;
use super::preset::{bar_spec, file_label, plot_presets};
use super::spec::{ChartKind, HeatmapData, PlotSpec, PlotStyle, Series};
use crate::parser::summary_csv::SummaryCsv;

// Rows that pool every input file in the per-species summary.
const ALL_FILES: &str = "all";

plot_presets! {
    /// Charts available for MAF summaries.
    MafPlot {
        SpeciesCoverage => ("species-coverage", "Reference coverage per species (per_species.csv)"),
        SpeciesGaps => ("species-gaps", "Gap proportion per species (per_species.csv)"),
        SpeciesBlocks => ("species-blocks", "Proportion of blocks with each species (per_species.csv)"),
        CoverageHeatmap => ("coverage-heatmap", "Reference coverage by species and file (per_species.csv)"),
        FileCoverage => ("file-coverage", "Reference genome coverage per file (per_file.csv)"),
        FileGaps => ("file-gaps", "Gap proportion per file (per_file.csv)"),
        FileBlocks => ("file-blocks", "Number of blocks per file (per_file.csv)"),
    }
}

/// Summary file read by a MAF chart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MafSource {
    PerFile,
    PerSpecies,
}

impl MafPlot {
    pub fn source(&self) -> MafSource {
        match self {
            Self::FileCoverage | Self::FileGaps | Self::FileBlocks => MafSource::PerFile,
            _ => MafSource::PerSpecies,
        }
    }
}

/// Write every MAF summary chart. Returns the written files.
pub fn plot_maf_summary(per_file: &Path, per_species: &Path, writer: &PlotWriter) -> Vec<PathBuf> {
    MafPlot::all()
        .iter()
        .filter_map(|preset| {
            let input = match preset.source() {
                MafSource::PerFile => per_file,
                MafSource::PerSpecies => per_species,
            };
            let spec = MafPlotter::new(input, writer.style()).build(preset);
            writer.write(preset.name(), spec)
        })
        .collect()
}

/// Builds MAF summary charts.
pub struct MafPlotter<'a> {
    input: &'a Path,
    style: &'a PlotStyle,
}

impl<'a> MafPlotter<'a> {
    pub fn new(input: &'a Path, style: &'a PlotStyle) -> Self {
        Self { input, style }
    }

    pub fn build(&self, preset: &MafPlot) -> Result<PlotSpec> {
        let csv = SummaryCsv::from_path(self.input)?;
        if csv.is_empty() {
            bail!("No rows in {}", self.input.display());
        }
        match preset {
            MafPlot::SpeciesCoverage => {
                self.species_bar(&csv, "ref_coverage", "Reference coverage")
            }
            MafPlot::SpeciesGaps => self.species_bar(&csv, "gap_proportion", "Gap proportion"),
            MafPlot::SpeciesBlocks => {
                self.species_bar(&csv, "block_proportion", "Proportion of blocks")
            }
            MafPlot::CoverageHeatmap => self.coverage_heatmap(&csv),
            MafPlot::FileCoverage => {
                self.file_bar(&csv, "ref_genome_coverage", "Reference genome coverage")
            }
            MafPlot::FileGaps => self.file_bar(&csv, "gap_proportion", "Gap proportion"),
            MafPlot::FileBlocks => self.file_bar(&csv, "blocks", "Number of blocks"),
        }
    }

    /// Species bar chart. Uses the pooled rows when there are several files.
    fn species_bar(&self, csv: &SummaryCsv, column: &str, label: &str) -> Result<PlotSpec> {
        let paths = csv.column_str("path")?;
        let pooled = paths.iter().any(|p| p == ALL_FILES);
        let keep: Vec<bool> = paths.iter().map(|p| !pooled || p == ALL_FILES).collect();
        let species = select(&csv.column_str("name")?, &keep);
        let values = select(&csv.column_f64(column)?, &keep);
        Ok(bar_spec(
            &format!("{} per species", label),
            ("Species", label),
            species,
            vec![Series::categorical(label, values)],
            self.style,
        ))
    }

    fn coverage_heatmap(&self, csv: &SummaryCsv) -> Result<PlotSpec> {
        let paths = csv.column_str("path")?;
        let species = csv.column_str("name")?;
        let coverage = csv.column_f64("ref_coverage")?;
        let files: IndexSet<&String> = paths.iter().filter(|p| *p != ALL_FILES).collect();
        let names: IndexSet<&String> = species.iter().collect();
        let mut cells = vec![f64::NAN; names.len() * files.len()];
        for ((path, name), value) in paths.iter().zip(&species).zip(&coverage) {
            if let (Some(col), Some(row)) = (files.get_index_of(path), names.get_index_of(name)) {
                cells[row * files.len() + col] = *value;
            }
        }
        let mut spec = PlotSpec::new(
            "Reference coverage by species",
            ChartKind::Heatmap,
            "File",
            "Species",
        );
        spec.heatmap = Some(HeatmapData::new(
            names.iter().map(|n| n.to_string()).collect(),
            files.iter().map(|f| file_label(f)).collect(),
            cells,
        )?);
        spec.style = self.style.clone();
        Ok(spec)
    }

    fn file_bar(&self, csv: &SummaryCsv, column: &str, label: &str) -> Result<PlotSpec> {
        let files = csv
            .column_str("path")?
            .iter()
            .map(|p| file_label(p))
            .collect();
        Ok(bar_spec(
            &format!("{} per file", label),
            ("File", label),
            files,
            vec![Series::categorical(label, csv.column_f64(column)?)],
            self.style,
        ))
    }
}

fn select<T: Clone>(values: &[T], keep: &[bool]) -> Vec<T> {
    values
        .iter()
        .zip(keep)
        .filter(|(_, k)| **k)
        .map(|(v, _)| v.clone())
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::plot::fonts::PlotFonts;
    use crate::plot::spec::PlotFormat;

    const SPECIES: &str = "tests/files/plot/maf_summary_per_species.csv";
    const FILES: &str = "tests/files/plot/maf_summary_per_file.csv";

    fn build(input: &str, preset: MafPlot) -> PlotSpec {
        let style = PlotStyle::default();
        let spec = MafPlotter::new(Path::new(input), &style)
            .build(&preset)
            .unwrap();
        spec.validate().unwrap();
        spec
    }

    #[test]
    fn test_all_presets_build() {
        for preset in MafPlot::all() {
            let input = match preset.source() {
                MafSource::PerFile => FILES,
                MafSource::PerSpecies => SPECIES,
            };
            build(input, *preset);
        }
    }

    #[test]
    fn test_plot_maf_summary() {
        let dir = tempfile::TempDir::with_prefix("plot").unwrap();
        let style = PlotStyle::default();
        let fonts = PlotFonts::new();
        let writer = PlotWriter::new(dir.path(), None, &PlotFormat::Svg, &style, &fonts);
        let written = plot_maf_summary(Path::new(FILES), Path::new(SPECIES), &writer);
        assert_eq!(written.len(), MafPlot::all().len());
        assert!(dir.path().join("file-coverage.svg").is_file());
    }

    #[test]
    fn test_species_coverage() {
        let spec = build(SPECIES, MafPlot::SpeciesCoverage);
        assert_eq!(spec.categories.len(), 5);
        assert_eq!(spec.series[0].y[4], 0.7544);
    }

    #[test]
    fn test_species_pooled_rows() {
        let csv = SummaryCsv::from_reader(
            "path,name,ref_coverage\na.maf,hg,1\nb.maf,hg,0.5\nall,hg,0.75\n".as_bytes(),
        )
        .unwrap();
        let style = PlotStyle::default();
        let plotter = MafPlotter::new(Path::new(SPECIES), &style);
        let spec = plotter
            .species_bar(&csv, "ref_coverage", "Reference coverage")
            .unwrap();
        assert_eq!(spec.series[0].y, vec![0.75]);
        let heatmap = plotter.coverage_heatmap(&csv).unwrap().heatmap.unwrap();
        assert_eq!(heatmap.cols, vec!["a", "b"]);
        assert_eq!(heatmap.get(0, 1), Some(0.5));
    }

    #[test]
    fn test_file_coverage() {
        let spec = build(FILES, MafPlot::FileCoverage);
        assert_eq!(spec.categories, vec!["simple"]);
        assert_eq!(spec.series[0].y, vec![3.5952e-7]);
    }
}

//! Write a batch of charts into one directory.
//!
//! Used by the summary commands to write every chart that applies to their
//! output. A chart that cannot be built (for example, GC content for amino
//! acid data) is skipped with a warning, so one missing column does not stop
//! the rest.
use std::path::{Path, PathBuf};

use anyhow::Result;

use super::fonts::PlotFonts;
use super::render::PlotRenderer;
use super::spec::{PlotFormat, PlotSpec, PlotStyle};

/// Renders charts to `<output_dir>/<prefix>_<name>.<ext>`.
pub struct PlotWriter<'a> {
    output_dir: &'a Path,
    prefix: Option<&'a str>,
    format: &'a PlotFormat,
    style: &'a PlotStyle,
    fonts: &'a PlotFonts,
}

impl<'a> PlotWriter<'a> {
    pub fn new(
        output_dir: &'a Path,
        prefix: Option<&'a str>,
        format: &'a PlotFormat,
        style: &'a PlotStyle,
        fonts: &'a PlotFonts,
    ) -> Self {
        Self {
            output_dir,
            prefix,
            format,
            style,
            fonts,
        }
    }

    pub fn style(&self) -> &'a PlotStyle {
        self.style
    }

    /// Render one chart. Returns `None`, with a warning, when the chart
    /// cannot be built or written.
    pub fn write(&self, name: &str, spec: Result<PlotSpec>) -> Option<PathBuf> {
        let output = self.output_path(name);
        let result =
            spec.and_then(|spec| PlotRenderer::new(&spec, self.fonts).render(&output, self.format));
        match result {
            Ok(()) => Some(output),
            Err(error) => {
                log::warn!("Skipped chart {}: {:#}", name, error);
                None
            }
        }
    }

    fn output_path(&self, name: &str) -> PathBuf {
        let file_name = match self.prefix {
            Some(prefix) => format!("{}_{}.{}", prefix, name, self.format.extension()),
            None => format!("{}.{}", name, self.format.extension()),
        };
        self.output_dir.join(file_name)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::plot::spec::{ChartKind, Series};

    #[test]
    fn test_write_and_skip() {
        let dir = tempfile::TempDir::with_prefix("plot").unwrap();
        let output_dir = dir.path().join("plots");
        let style = PlotStyle::default();
        let fonts = PlotFonts::new();
        let writer = PlotWriter::new(&output_dir, Some("run"), &PlotFormat::Svg, &style, &fonts);
        let mut spec = PlotSpec::new("t", ChartKind::Scatter, "x", "y");
        spec.series.push(Series::xy("s", vec![1.0], vec![1.0]));
        let path = writer.write("scatter", Ok(spec)).unwrap();
        assert_eq!(path, output_dir.join("run_scatter.svg"));
        assert!(path.is_file());
        assert!(
            writer
                .write("broken", Err(anyhow::anyhow!("missing column")))
                .is_none()
        );
        assert!(!output_dir.join("run_broken.svg").exists());
    }
}

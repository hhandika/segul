//! Shared handling for the `--plot` flag of the summary commands.
use std::path::{Path, PathBuf};

use colored::Colorize;
use segul_core::plot::batch::PlotWriter;
use segul_core::plot::fonts::PlotFonts;
use segul_core::plot::spec::{PlotFormat, PlotStyle};

use crate::cli::args::PlotFlags;
use crate::helper::utils;

const PLOT_DIR: &str = "plots";

pub(in crate::cli) trait PlotCli {
    fn plot_style(&self, flags: &PlotFlags) -> PlotStyle {
        PlotStyle {
            palette: flags.palette,
            colormap: flags.colormap,
            cvd: flags.cvd,
            ..Default::default()
        }
    }

    /// Fonts for PNG text: the `--font` file, otherwise the system fonts.
    /// SVG uses the viewer's fonts, so none are loaded.
    fn plot_fonts(&self, flags: &PlotFlags) -> Option<PlotFonts> {
        if flags.plot_format == PlotFormat::Svg {
            return Some(PlotFonts::new());
        }
        let fonts = match &flags.font {
            Some(path) => {
                let mut fonts = PlotFonts::new();
                if let Err(error) = fonts.add_font_file(path) {
                    log::error!("{:#}", error);
                    return None;
                }
                fonts
            }
            None => PlotFonts::system(),
        };
        if fonts.is_empty() {
            log::error!("No fonts found for PNG charts. Use --font <FILE> or --plot-format svg");
            return None;
        }
        Some(fonts)
    }

    /// Write the charts produced by `plot` into `<output>/plots`.
    /// Does nothing unless `--plot` is set.
    fn write_plots<F>(&self, flags: &PlotFlags, output: &Path, prefix: Option<&str>, plot: F)
    where
        F: FnOnce(&PlotWriter) -> Vec<PathBuf>,
    {
        if !flags.plot {
            return;
        }
        let Some(fonts) = self.plot_fonts(flags) else {
            log::warn!("Skipped writing charts. The summary files are complete.");
            return;
        };
        let plot_dir = output.join(PLOT_DIR);
        let style = self.plot_style(flags);
        let writer = PlotWriter::new(&plot_dir, prefix, &flags.plot_format, &style, &fonts);
        let spin = utils::set_spinner();
        spin.set_message("Writing charts...");
        let written = plot(&writer);
        spin.finish_with_message("Finished writing charts!\n");
        log::info!("{}", "Plots".yellow());
        log::info!("{:18}: {}", "Plot dir", plot_dir.display());
        log::info!("{:18}: {}", "Format", flags.plot_format);
        log::info!("{:18}: {}", "Palette", flags.palette);
        log::info!("{:18}: {}", "Colormap", flags.colormap);
        if let Some(cvd) = &flags.cvd {
            log::info!("{:18}: {}", "Simulated CVD", cvd);
        }
        if let Some(family) = fonts.family() {
            log::info!("{:18}: {}", "Font", family);
        }
        log::info!("{:18}: {}\n", "Charts", utils::fmt_num(&written.len()));
    }
}

//! Fonts for PNG and PDF output.
//!
//! SVG output names the generic `sans-serif` family and needs no fonts.
//! PNG and PDF draw text with the fonts held here: font files supplied by
//! the caller (for example, a font bundled with an app), the system fonts,
//! or both. Fonts added by the caller take precedence over system fonts.
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result, bail};
use usvg::fontdb::{self, Database};

// Common sans-serif families, tried in order when only system fonts are loaded.
const SYSTEM_SANS_SERIF: [&str; 7] = [
    "DejaVu Sans",
    "Noto Sans",
    "Liberation Sans",
    "Helvetica",
    "Arial",
    "Segoe UI",
    "Roboto",
];

/// Fonts used to draw text in PNG and PDF plots.
///
/// # Example
/// ```rust
/// use segul_core::plot::fonts::PlotFonts;
///
/// // In an app, pass the bytes of a bundled font instead.
/// let fonts = PlotFonts::system();
/// let mut custom = PlotFonts::new();
/// assert!(custom.is_empty());
/// assert!(custom.add_font_data(Vec::new()).is_err());
/// ```
#[derive(Debug, Clone, Default)]
pub struct PlotFonts {
    db: Arc<Database>,
    // Family used for the SVG's generic `sans-serif`.
    family: Option<String>,
}

impl PlotFonts {
    /// An empty font set. Add fonts before rendering PNG or PDF.
    pub fn new() -> Self {
        Self::default()
    }

    /// Fonts installed on the system.
    pub fn system() -> Self {
        let mut fonts = Self::new();
        fonts.add_system_fonts();
        fonts
    }

    /// Fonts from TrueType or OpenType file contents.
    /// The first font sets the family used for plot text.
    pub fn from_data(fonts: Vec<Vec<u8>>) -> Result<Self> {
        let mut plot_fonts = Self::new();
        for data in fonts {
            plot_fonts.add_font_data(data)?;
        }
        Ok(plot_fonts)
    }

    /// Add a font from TrueType or OpenType file contents.
    /// The first font added this way sets the family used for plot text.
    pub fn add_font_data(&mut self, data: Vec<u8>) -> Result<()> {
        let db = Arc::make_mut(&mut self.db);
        let ids = db.load_font_source(fontdb::Source::Binary(Arc::new(data)));
        let family = ids
            .first()
            .and_then(|id| db.face(*id))
            .and_then(|face| face.families.first())
            .map(|(name, _)| name.clone());
        match family {
            Some(family) => {
                if self.family.is_none() {
                    db.set_sans_serif_family(family.as_str());
                    self.family = Some(family);
                }
                Ok(())
            }
            None => bail!("Failed reading font data. Use a TrueType or OpenType font"),
        }
    }

    /// Add a font from a TrueType or OpenType file.
    pub fn add_font_file(&mut self, path: &Path) -> Result<()> {
        let data =
            std::fs::read(path).with_context(|| format!("Failed reading {}", path.display()))?;
        self.add_font_data(data)
            .with_context(|| format!("Failed loading font {}", path.display()))
    }

    /// Add the fonts installed on the system.
    pub fn add_system_fonts(&mut self) {
        let db = Arc::make_mut(&mut self.db);
        db.load_system_fonts();
        if self.family.is_none() {
            let family = SYSTEM_SANS_SERIF
                .iter()
                .find(|name| has_family(db, name))
                .map(|name| name.to_string())
                .or_else(|| {
                    db.faces()
                        .next()
                        .and_then(|face| face.families.first())
                        .map(|(name, _)| name.clone())
                });
            if let Some(family) = &family {
                db.set_sans_serif_family(family.as_str());
            }
            self.family = family;
        }
    }

    pub fn is_empty(&self) -> bool {
        self.db.is_empty()
    }

    /// Family used for plot text, if any font is loaded.
    pub fn family(&self) -> Option<&str> {
        self.family.as_deref()
    }

    pub(super) fn svg_options(&self) -> Result<usvg::Options<'static>> {
        let Some(family) = &self.family else {
            bail!("No fonts available for plot text. Add a font file or the system fonts");
        };
        Ok(usvg::Options {
            font_family: family.clone(),
            fontdb: self.db.clone(),
            ..Default::default()
        })
    }
}

fn has_family(db: &Database, name: &str) -> bool {
    db.faces()
        .any(|face| face.families.iter().any(|(family, _)| family == name))
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_empty_fonts() {
        let fonts = PlotFonts::new();
        assert!(fonts.is_empty());
        assert!(fonts.family().is_none());
        assert!(fonts.svg_options().is_err());
    }

    #[test]
    fn test_invalid_font_data() {
        assert!(PlotFonts::from_data(vec![b"not a font".to_vec()]).is_err());
    }

    #[test]
    fn test_system_fonts() {
        let fonts = PlotFonts::system();
        if !fonts.is_empty() {
            assert!(fonts.family().is_some());
            assert!(fonts.svg_options().is_ok());
        }
    }
}

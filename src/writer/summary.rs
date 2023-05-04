use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use colored::Colorize;

use crate::helper::alphabet;
use crate::helper::types::{DataType, TaxonRecords};
use crate::helper::utils;
use crate::stats::sequence::{CharMatrix, CharSummary, Completeness, SiteSummary, Sites, Taxa};
use crate::writer::FileWriter;

trait Alphabet {
    fn match_alphabet(&self, datatype: &DataType) -> &str {
        match datatype {
            DataType::Dna => alphabet::DNA_STR_UPPERCASE,
            DataType::Aa => alphabet::AA_STR_UPPERCASE,
            _ => unreachable!("Invalid data types! Use dna or aa only"),
        }
    }
}

impl Alphabet for CsvWriter<'_> {}
impl FileWriter for CsvWriter<'_> {}

pub struct CsvWriter<'a> {
    output: &'a Path,
    prefix: &'a Option<String>,
    datatype: &'a DataType,
}

impl<'a> CsvWriter<'a> {
    pub fn new(output: &'a Path, prefix: &'a Option<String>, datatype: &'a DataType) -> Self {
        Self {
            output,
            prefix,
            datatype,
        }
    }

    pub fn write_per_locus_summary(&self, file: &Path, summary: &Taxa) -> Result<()> {
        let default_prefix = file
            .file_stem()
            .and_then(OsStr::to_str)
            .expect("Failed to parse input file stem");
        let output = self.create_output_fname(default_prefix);
        let mut writer = self.create_output_file(&output)?;
        write!(
            writer,
            "taxon,\
        total_chars, \
        missing_data,\
        proportion_missing_data\
    "
        )?;
        if DataType::Dna == *self.datatype {
            write!(
                writer,
                ",gc_content\
                ,at_content\
                ,nucleotides"
            )?;
        }
        let alphabet = self.match_alphabet(self.datatype);
        self.write_alphabet_header(&mut writer, alphabet)?;
        summary.records.iter().for_each(|(taxon, chars)| {
            write!(
                writer,
                "{},{},{}",
                taxon, chars.total_chars, chars.missing_data
            )
            .expect("Failed taxon summary stats");
            write!(
                writer,
                ",{}",
                chars.missing_data as f64 / chars.total_chars as f64
            )
            .expect("Failed to write taxon summary stats");
            if DataType::Dna == *self.datatype {
                write!(
                    writer,
                    ",{},{},{}",
                    chars.gc_count as f64 / chars.total_chars as f64,
                    chars.at_count as f64 / chars.total_chars as f64,
                    chars.nucleotides as f64 / chars.total_chars as f64
                )
                .expect("Failed to write taxon summary stats");
            }
            alphabet.chars().for_each(|ch| {
                write!(writer, ",{}", chars.chars.get(&ch).unwrap_or(&0))
                    .expect("Failed getting character summary stats");
            });
            writeln!(writer).expect("Failed writing per locus summary stats");
        });
        Ok(())
    }

    pub fn write_taxon_summary(
        &self,
        taxon_summary: &BTreeMap<String, TaxonRecords>,
    ) -> Result<()> {
        let default_prefix = "taxon_summary";
        let output = self.create_output_fname(default_prefix);
        let mut writer = self.create_output_file(&output)?;
        write!(
            writer,
            "taxon,\
        locus_counts, \
        total_chars, \
        missing_data,\
        proportion_missing_data\
    "
        )?;
        if DataType::Dna == *self.datatype {
            write!(
                writer,
                ",gc_content\
            ,at_content\
            ,nucleotides"
            )?;
        }
        let alphabet = self.match_alphabet(self.datatype);
        self.write_alphabet_header(&mut writer, alphabet)?;
        taxon_summary.iter().for_each(|(taxon, counts)| {
            write!(
                writer,
                "{},{},{},{}",
                taxon, counts.locus_counts, counts.total_chars, counts.missing_data
            )
            .expect("Failed taxon summary stats");
            write!(
                writer,
                ",{}",
                counts.missing_data as f64 / counts.total_chars as f64
            )
            .expect("Failed to write taxon summary stats");
            if DataType::Dna == *self.datatype {
                write!(
                    writer,
                    ",{},{},{}",
                    counts.gc_count as f64 / counts.total_chars as f64,
                    counts.at_count as f64 / counts.total_chars as f64,
                    counts.nucleotides as f64 / counts.total_chars as f64
                )
                .expect("Failed to write taxon summary stats");
            }
            alphabet.chars().for_each(|ch| {
                write!(writer, ",{}", counts.char_counts.get(&ch).unwrap_or(&0))
                    .expect("Failed taxon summary stats");
            });
            writeln!(writer).expect("Failed taxon summary stats");
        });
        writer.flush()?;
        Ok(())
    }

    pub fn write_locus_summary(&self, stats: &[(Sites, CharMatrix, Taxa)]) -> Result<()> {
        let default_prefix = "locus_summary";
        let output = self.create_output_fname(default_prefix);
        let mut writer = self.create_output_file(&output)?;
        let alphabet = self.match_alphabet(self.datatype);
        self.write_locus_header(&mut writer, alphabet)?;
        stats.iter().for_each(|(site, chars, _)| {
            self.write_locus_content(&mut writer, site, chars, alphabet)
                .unwrap();
        });
        writer.flush()?;
        log::info!("{}", "Output".yellow());
        log::info!("{:18}: {}", "Output dir", self.output.display());

        Ok(())
    }

    fn create_output_fname(&self, default_prefix: &str) -> PathBuf {
        match self.prefix {
            Some(fname) => {
                let out_name = format!("{}_{}", fname, default_prefix);
                self.output.join(out_name).with_extension("csv")
            }
            None => self.output.join(default_prefix).with_extension("csv"),
        }
    }

    fn write_locus_header<W: Write>(&self, writer: &mut W, alphabet: &str) -> Result<()> {
        write!(
            writer,
            "path,\
            locus,\
            ntax,\
            chars_count,\
            site_count,\
            conserved_sites,\
            proportion_cons_sites,\
            variable_sites,\
            proportion_var_sites,\
            parsimony_informative_sites,\
            proportion_pars_inf_sites,\
            missing_data,\
            proportion_missing_data\
        "
        )?;
        if DataType::Dna == *self.datatype {
            write!(
                writer,
                ",gc_content\
                ,at_content\
                ,nucleotides"
            )?;
        }
        self.write_alphabet_header(writer, alphabet)?;
        Ok(())
    }

    fn write_alphabet_header<W: Write>(&self, writer: &mut W, alphabet: &str) -> Result<()> {
        alphabet
            .chars()
            .for_each(|ch| write!(writer, ",{}", ch).unwrap());
        writeln!(writer)?;
        Ok(())
    }

    fn write_locus_content<W: Write>(
        &self,
        writer: &mut W,
        site: &Sites,
        chars: &CharMatrix,
        alphabet: &str,
    ) -> Result<()> {
        write!(
            writer,
            "{},{},{},{},",
            site.path.display(),
            site.path
                .file_stem()
                .and_then(OsStr::to_str)
                .with_context(|| format!(
                    "Failed getting locus name for {}",
                    site.path.display()
                ))?,
            chars.ntax,
            chars.chars.total_chars
        )?;

        // Site stats
        write!(
            writer,
            "{},{},{},{},{},{},{},",
            site.counts,
            site.conserved,
            site.prop_cons,
            site.variable,
            site.prop_var,
            site.pars_inf,
            site.prop_pinf
        )?;

        // Missing data
        write!(writer, "{},", chars.chars.missing_data)?;
        write!(writer, "{}", chars.chars.prop_missing_data)?;

        // We move comma position for accurate printing
        if DataType::Dna == *self.datatype {
            // GC content
            write!(
                writer,
                ",{}",
                chars.chars.gc_count as f64 / chars.chars.total_chars as f64
            )?;
            // AT content
            write!(
                writer,
                ",{}",
                chars.chars.at_count as f64 / chars.chars.total_chars as f64
            )?;
            write!(writer, ",{}", chars.chars.nucleotides)?;
        }

        // Characters
        alphabet.chars().for_each(|ch| {
            write!(writer, ",{}", chars.chars.chars.get(&ch).unwrap_or(&0)).unwrap();
        });
        writeln!(writer)?;
        Ok(())
    }
}

impl Alphabet for SummaryWriter<'_> {}

pub struct SummaryWriter<'s> {
    site: &'s SiteSummary,
    chars: &'s CharSummary,
    complete: &'s Completeness,
    datatype: &'s DataType,
}

impl<'s> SummaryWriter<'s> {
    pub fn new(
        site: &'s SiteSummary,
        chars: &'s CharSummary,
        complete: &'s Completeness,
        datatype: &'s DataType,
    ) -> Self {
        Self {
            site,
            chars,
            complete,
            datatype,
        }
    }

    pub fn print_summary(&self) -> Result<()> {
        log::info!("\n{}", "General Summary".yellow());
        self.write_gen_sum();
        log::info!("{}", "Alignment Summary".yellow());
        self.write_aln_sum();
        log::info!("{}", "Taxon Summary".yellow());
        self.write_tax_sum();

        log::info!("{}", "Character Count".yellow());
        self.write_char_count();

        log::info!("{}", "Data Matrix Completeness".yellow());
        self.write_matrix_comp();

        log::info!("{}", "Conserved Sequences".yellow());
        self.write_cons_seq();

        log::info!("{}", "Variable Sequences".yellow());
        self.write_var_seq();

        log::info!("{}", "Parsimony Informative".yellow());
        self.write_pars_inf();
        Ok(())
    }

    fn write_gen_sum(&self) {
        log::info!(
            "{:18}: {}",
            "Total taxa",
            utils::fmt_num(&self.complete.total_tax)
        );
        log::info!(
            "{:18}: {}",
            "Total loci",
            utils::fmt_num(&self.site.total_loci)
        );
        log::info!(
            "{:18}: {}",
            "Total sites",
            utils::fmt_num(&self.site.total_sites)
        );
        log::info!(
            "{:18}: {}",
            "Missing data",
            utils::fmt_num(&self.chars.missing_data)
        );
        log::info!(
            "{:18}: {:.2}%",
            "%Missing data",
            &self.chars.prop_missing_data * 100.0
        );

        match self.datatype {
            DataType::Dna => self.write_dna_sum(),
            DataType::Aa => log::info!(
                "{:18}: {}\n",
                "Characters",
                utils::fmt_num(&self.chars.total_chars)
            ),
            _ => panic!("Please specify datatype"),
        }
    }

    fn write_dna_sum(&self) {
        log::info!("{:18}: {:.2}", "GC content", self.chars.gc_content);
        log::info!("{:18}: {:.2}", "AT content", self.chars.at_content);
        log::info!(
            "{:18}: {}",
            "Characters",
            utils::fmt_num(&self.chars.total_chars)
        );
        log::info!(
            "{:18}: {}\n",
            "Nucleotides",
            utils::fmt_num(&self.chars.total_nucleotides)
        );
    }

    fn write_aln_sum(&self) {
        log::info!(
            "{:18}: {} bp",
            "Min length",
            utils::fmt_num(&self.site.min_sites)
        );
        log::info!(
            "{:18}: {} bp",
            "Max length",
            utils::fmt_num(&self.site.max_sites)
        );
        log::info!("{:18}: {:.2} bp\n", "Mean length", &self.site.mean_sites);
    }

    fn write_tax_sum(&self) {
        log::info!("{:18}: {}", "Min taxa", utils::fmt_num(&self.chars.min_tax));
        log::info!("{:18}: {}", "Max taxa", utils::fmt_num(&self.chars.max_tax));
        log::info!("{:18}: {:.2}\n", "Mean taxa", self.chars.mean_tax);
    }

    fn write_char_count(&self) {
        let alphabet = self.match_alphabet(self.datatype);
        alphabet.chars().for_each(|ch| {
            if let Some(count) = self.chars.chars.get(&ch) {
                log::info!("{:18}: {}", ch, utils::fmt_num(count));
            }
        });
        log::info!("");
    }

    fn write_matrix_comp(&self) {
        self.complete
            .completeness
            .iter()
            .for_each(|(percent, ntax)| {
                let percent_str = format!("{}% taxa", percent);
                log::info!("{:18}: {}", percent_str, utils::fmt_num(ntax))
            });
        log::info!("");
    }

    fn write_cons_seq(&self) {
        log::info!(
            "{:18}: {}",
            "Con. loci",
            utils::fmt_num(&self.site.cons_loci)
        );
        log::info!(
            "{:18}: {:.2}%",
            "%Con. loci",
            self.site.prop_cons_loci * 100.0
        );
        log::info!(
            "{:18}: {}",
            "Con. sites",
            utils::fmt_num(&self.site.total_cons_site)
        );
        log::info!("{:18}: {:.2}%", "%Con. sites", &self.site.prop_cons_site);
        log::info!(
            "{:18}: {}",
            "Min con. sites",
            utils::fmt_num(&self.site.min_cons_site)
        );
        log::info!(
            "{:18}: {}",
            "Max con. sites",
            utils::fmt_num(&self.site.max_cons_site)
        );
        log::info!(
            "{:18}: {:.2}\n",
            "Mean con. sites",
            &self.site.mean_cons_site
        );
    }

    fn write_var_seq(&self) {
        log::info!(
            "{:18}: {}",
            "Var. loci",
            utils::fmt_num(&self.site.var_loci)
        );
        log::info!(
            "{:18}: {:.2}%",
            "%Var. loci",
            self.site.prop_var_loci * 100.0
        );
        log::info!(
            "{:18}: {}",
            "Var. sites",
            utils::fmt_num(&self.site.total_var_site)
        );
        log::info!("{:18}: {:.2}%", "%Var. sites", &self.site.prop_var_site);
        log::info!(
            "{:18}: {}",
            "Min var. sites",
            utils::fmt_num(&self.site.min_var_site)
        );
        log::info!(
            "{:18}: {}",
            "Max var. sites",
            utils::fmt_num(&self.site.max_var_site)
        );
        log::info!(
            "{:18}: {:.2}\n",
            "Mean var. sites",
            &self.site.mean_var_site
        );
    }

    fn write_pars_inf(&self) {
        log::info!(
            "{:18}: {}",
            "Inf. loci",
            utils::fmt_num(&self.site.inf_loci)
        );
        log::info!(
            "{:18}: {:.2}%",
            "%Inf. loci",
            self.site.prop_inf_loci * 100.0
        );
        log::info!(
            "{:18}: {}",
            "Inf. sites",
            utils::fmt_num(&self.site.total_inf_site)
        );
        log::info!("{:18}: {:.2}%", "%Inf. sites", &self.site.prop_inf_site);
        log::info!(
            "{:18}: {}",
            "Min inf. sites",
            utils::fmt_num(&self.site.min_inf_site)
        );
        log::info!(
            "{:18}: {}",
            "Max inf. sites",
            utils::fmt_num(&self.site.max_inf_site)
        );
        log::info!(
            "{:18}: {:.2}\n",
            "Mean inf. sites",
            &self.site.mean_inf_site
        );
    }
}

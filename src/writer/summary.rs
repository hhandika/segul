//! Summary writer for alignment files
use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Ok, Result};
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
    prefix: Option<&'a str>,
    datatype: &'a DataType,
}

impl<'a> CsvWriter<'a> {
    pub fn new(output: &'a Path, prefix: Option<&'a str>, datatype: &'a DataType) -> Self {
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
                ",GC_content\
                ,AT_content\
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

impl FileWriter for SummaryWriter<'_> {}

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

    pub fn write(&self, output_path: &Path, prefix: Option<&str>) -> Result<()> {
        let default_prefix = "alignment_summary";
        let output = match prefix {
            Some(fname) => {
                let out_name = format!("{}_{}", fname, default_prefix);
                output_path.join(out_name).with_extension("txt")
            }
            None => output_path.join(default_prefix).with_extension("txt"),
        };
        let mut writer = self.create_output_file(&output)?;
        self.write_to_file(&mut writer)?;
        writer.flush()?;
        Ok(())
    }

    fn write_to_file<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.write_gen_sum(writer)?;
        self.write_aln_sum(writer)?;
        self.write_tax_sum(writer)?;
        self.write_char_count(writer)?;
        self.write_matrix_comp(writer)?;
        self.write_cons_seq(writer)?;
        self.write_var_seq(writer)?;
        self.write_pars_inf(writer)?;
        Ok(())
    }

    fn write_gen_sum<W: Write>(&self, writer: &mut W) -> Result<()> {
        let summary = "General Summary";
        writeln!(writer, "{}", summary)?;
        let total_taxa = format!(
            "{:18}: {}",
            "Total taxa",
            utils::fmt_num(&self.complete.total_tax)
        );
        writeln!(writer, "{}", total_taxa)?;

        let total_loci = format!(
            "{:18}: {}",
            "Total loci",
            utils::fmt_num(&self.site.total_loci)
        );
        writeln!(writer, "{}", total_loci)?;

        let total_sites = format!(
            "{:18}: {}",
            "Total sites",
            utils::fmt_num(&self.site.total_sites)
        );

        writeln!(writer, "{}", total_sites)?;

        let total_chars = format!(
            "{:18}: {}",
            "Total characters",
            utils::fmt_num(&self.chars.total_chars)
        );
        writeln!(writer, "{}", total_chars)?;

        let missing_data = format!(
            "{:18}: {}",
            "Missing data",
            utils::fmt_num(&self.chars.missing_data)
        );
        writeln!(writer, "{}", missing_data)?;
        let prop_missing_data = format!(
            "{:18}: {:.2}%",
            "%Missing data",
            self.chars.prop_missing_data * 100.0
        );
        writeln!(writer, "{}", prop_missing_data)?;

        match self.datatype {
            DataType::Dna => self.write_dna_sum(writer)?,
            DataType::Aa => self.write_aa_sum(writer)?,
            _ => unreachable!("Invalid data types! Use dna or aa only"),
        };

        Ok(())
    }

    fn write_aa_sum<W: Write>(&self, writer: &mut W) -> Result<()> {
        let characters = format!(
            "{:18}: {}",
            "Characters",
            utils::fmt_num(&self.chars.total_chars)
        );
        writeln!(writer, "{}", characters)?;

        Ok(())
    }

    fn write_dna_sum<W: Write>(&self, writer: &mut W) -> Result<()> {
        let gc_content = format!("{:18}: {:.2}", "GC content", self.chars.gc_content);
        writeln!(writer, "{}", gc_content)?;

        let at_content = format!("{:18}: {:.2}", "AT content", self.chars.at_content);

        writeln!(writer, "{}", at_content)?;

        let characters = format!(
            "{:18}: {}",
            "Characters",
            utils::fmt_num(&self.chars.total_chars)
        );
        writeln!(writer, "{}", characters)?;

        let nucleotides = format!(
            "{:18}: {}",
            "Nucleotides",
            utils::fmt_num(&self.chars.total_nucleotides)
        );
        writeln!(writer, "{}", nucleotides)?;

        Ok(())
    }

    fn write_aln_sum<W: Write>(&self, writer: &mut W) -> Result<()> {
        let aln_summary = "Alignment Summary";
        writeln!(writer, "\n{}", aln_summary)?;

        let min_site = format!(
            "{:18}: {} bp",
            "Min length",
            utils::fmt_num(&self.site.min_sites)
        );
        writeln!(writer, "{}", min_site)?;

        let max_site = format!(
            "{:18}: {} bp",
            "Max length",
            utils::fmt_num(&self.site.max_sites)
        );
        writeln!(writer, "{}", max_site)?;

        let mean_site = format!("{:18}: {:.2} bp", "Mean length", &self.site.mean_sites);

        writeln!(writer, "{}", mean_site)?;

        Ok(())
    }

    fn write_tax_sum<W: Write>(&self, writer: &mut W) -> Result<()> {
        let taxon_summary = "Taxon Summary";
        writeln!(writer, "\n{}", taxon_summary)?;

        let min_tax = format!("{:18}: {}", "Min taxa", utils::fmt_num(&self.chars.min_tax));
        writeln!(writer, "{}", min_tax)?;

        let max_tax = format!("{:18}: {}", "Max taxa", utils::fmt_num(&self.chars.max_tax));
        writeln!(writer, "{}", max_tax)?;

        let mean_tax = format!("{:18}: {:.2}", "Mean taxa", self.chars.mean_tax);
        writeln!(writer, "{}", mean_tax)?;

        Ok(())
    }

    fn write_char_count<W: Write>(&self, writer: &mut W) -> Result<()> {
        let char_count = "Character Count";
        writeln!(writer, "\n{}", char_count)?;
        let alphabet = self.match_alphabet(self.datatype);

        for ch in alphabet.chars() {
            if let Some(count) = self.chars.chars.get(&ch) {
                let count = format!("{:18}: {}", ch, utils::fmt_num(count));
                writeln!(writer, "{}", count)?;
            }
        }
        Ok(())
    }

    fn write_matrix_comp<W: Write>(&self, writer: &mut W) -> Result<()> {
        let matrix_summary = "Matrix Completeness";
        writeln!(writer, "\n{}", matrix_summary)?;
        self.complete
            .completeness
            .iter()
            .for_each(|(percent, ntax)| {
                let percent_str = format!("{}% taxa", percent);
                let percent_str = format!("{:18}: {}", percent_str, utils::fmt_num(ntax));
                writeln!(writer, "{}", percent_str).expect("Failed to write matrix completeness");
            });
        Ok(())
    }

    fn write_cons_seq<W: Write>(&self, writer: &mut W) -> Result<()> {
        let conserved_seq = "Conserved Sequences";
        writeln!(writer, "\n{}", conserved_seq)?;

        let conserved_loci = format!(
            "{:18}: {}",
            "Conserved loci",
            utils::fmt_num(&self.site.cons_loci)
        );

        writeln!(writer, "{}", conserved_loci)?;

        let prop_cons_loci = format!(
            "{:18}: {:.2}%",
            "%Con. loci",
            self.site.prop_cons_loci * 100.0
        );

        writeln!(writer, "{}", prop_cons_loci)?;

        let conserved_sites = format!(
            "{:18}: {}",
            "Conserved sites",
            utils::fmt_num(&self.site.total_cons_site)
        );

        writeln!(writer, "{}", conserved_sites)?;

        let prop_cons_sites = format!(
            "{:18}: {:.2}%",
            "%Con. sites",
            &self.site.prop_cons_site * 100.0
        );

        writeln!(writer, "{}", prop_cons_sites)?;

        let min_cons_site = format!(
            "{:18}: {}",
            "Min con. sites",
            utils::fmt_num(&self.site.min_cons_site)
        );

        writeln!(writer, "{}", min_cons_site)?;

        let max_cons_site = format!(
            "{:18}: {}",
            "Max con. sites",
            utils::fmt_num(&self.site.max_cons_site)
        );

        writeln!(writer, "{}", max_cons_site)?;

        let mean_cons_site = format!("{:18}: {:.2}", "Mean con. sites", &self.site.mean_cons_site);

        writeln!(writer, "{}", mean_cons_site)?;

        Ok(())
    }

    fn write_var_seq<W: Write>(&self, writer: &mut W) -> Result<()> {
        let var_seq = "Variable Sequences";

        writeln!(writer, "\n{}", var_seq)?;

        let var_loci = format!(
            "{:18}: {}",
            "Var. loci",
            utils::fmt_num(&self.site.var_loci)
        );

        writeln!(writer, "{}", var_loci)?;

        let prop_var_loci = format!(
            "{:18}: {:.2}%",
            "%Var. loci",
            self.site.prop_var_loci * 100.0
        );

        writeln!(writer, "{}", prop_var_loci)?;

        let var_sites = format!(
            "{:18}: {}",
            "Var. sites",
            utils::fmt_num(&self.site.total_var_site)
        );

        writeln!(writer, "{}", var_sites)?;

        let prop_var_sites = format!(
            "{:18}: {:.2}%",
            "%Var. sites",
            &self.site.prop_var_site * 100.0
        );

        writeln!(writer, "{}", prop_var_sites)?;

        let min_var_site = format!(
            "{:18}: {}",
            "Min var. sites",
            utils::fmt_num(&self.site.min_var_site)
        );

        writeln!(writer, "{}", min_var_site)?;

        let max_var_site = format!(
            "{:18}: {}",
            "Max var. sites",
            utils::fmt_num(&self.site.max_var_site)
        );

        writeln!(writer, "{}", max_var_site)?;

        let mean_var_site = format!("{:18}: {:.2}", "Mean var. sites", &self.site.mean_var_site);

        writeln!(writer, "{}", mean_var_site)?;

        Ok(())
    }

    fn write_pars_inf<W: Write>(&self, writer: &mut W) -> Result<()> {
        let pars_inf = "Parsimony Informative";
        writeln!(writer, "\n{}", pars_inf)?;

        let inf_loci = format!(
            "{:18}: {}",
            "Inf. loci",
            utils::fmt_num(&self.site.inf_loci)
        );

        writeln!(writer, "{}", inf_loci)?;

        let prop_inf_loci = format!(
            "{:18}: {:.2}%",
            "%Inf. loci",
            self.site.prop_inf_loci * 100.0
        );

        writeln!(writer, "{}", prop_inf_loci)?;

        let inf_sites = format!(
            "{:18}: {}",
            "Inf. sites",
            utils::fmt_num(&self.site.total_inf_site)
        );

        writeln!(writer, "{}", inf_sites)?;

        let prop_inf_sites = format!(
            "{:18}: {:.2}%",
            "%Inf. sites",
            self.site.prop_inf_site * 100.0
        );

        writeln!(writer, "{}", prop_inf_sites)?;

        let min_inf_site = format!(
            "{:18}: {}",
            "Min inf. sites",
            utils::fmt_num(&self.site.min_inf_site)
        );

        writeln!(writer, "{}", min_inf_site)?;

        let max_inf_site = format!(
            "{:18}: {}",
            "Max inf. sites",
            utils::fmt_num(&self.site.max_inf_site)
        );

        writeln!(writer, "{}", max_inf_site)?;

        let mean_inf_site = format!("{:18}: {:.2}", "Mean inf. sites", &self.site.mean_inf_site);
        writeln!(writer, "{}", mean_inf_site)?;

        log::info!("");

        Ok(())
    }
}

use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufWriter, Write};
// use std::path::PathBuf;

use ansi_term::Colour::Yellow;
use anyhow::{Context, Result};

use crate::helper::alphabet;
use crate::helper::stats::{CharSummary, Chars, Completeness, SiteSummary, Sites};
use crate::helper::types::DataType;
use crate::helper::utils;

pub fn print_stats(site: &Sites, dna: &Chars) {
    log::info!("{}", Yellow.paint("Alignment"));
    log::info!("{:18}: {}", "Taxa", utils::fmt_num(&dna.ntax));
    log::info!("{:18}: {}", "Length", utils::fmt_num(&site.counts));

    log::info!("\n{}", Yellow.paint("Characters"));
    log::info!("{:18}: {}", "Total", utils::fmt_num(&dna.total_chars));
    log::info!(
        "{:18}: {}",
        "Missing data",
        utils::fmt_num(&dna.missing_data)
    );
    log::info!(
        "{:18}: {:.2}%",
        "Prop. missing",
        &dna.prop_missing_data * 100.0
    );

    dna.chars
        .iter()
        .for_each(|(ch, count)| log::info!("{:18}: {}", ch, utils::fmt_num(count)));

    log::info!("\n{}", Yellow.paint("Sites"));
    log::info!("{:18}: {}", "Conserved", utils::fmt_num(&site.conserved));
    log::info!("{:18}: {}", "Variable", utils::fmt_num(&site.variable));
    log::info!(
        "{:18}: {}",
        "Parsimony inf.",
        utils::fmt_num(&site.pars_inf)
    );
    log::info!("{:18}: {:.2}%", "%Conserved", site.prop_cons * 100.0);
    log::info!("{:18}: {:.2}%", "%Variable", site.prop_var * 100.0);
    log::info!("{:18}: {:.2}%\n", "%Pars. inf.", site.prop_var * 100.0);
}

trait Alphabet {
    fn get_alphabet(&self, datatype: &DataType) -> &str {
        match datatype {
            DataType::Dna => alphabet::DNA_STR_UPPERCASE,
            DataType::Aa => alphabet::AA_STR_UPPERCASE,
            _ => unreachable!(),
        }
    }
}

impl Alphabet for CsvWriter<'_> {}

pub struct CsvWriter<'a> {
    output: String,
    datatype: &'a DataType,
}

impl<'a> CsvWriter<'a> {
    pub fn new(output: &str, datatype: &'a DataType) -> Self {
        Self {
            output: String::from(output),
            datatype,
        }
    }

    pub fn write_summary_dir(&mut self, stats: &[(Sites, Chars)]) -> Result<()> {
        self.get_ouput_fname();
        let file = File::create(&self.output)
            .with_context(|| format!("Failed creating file {}", self.output))?;
        let mut writer = BufWriter::new(file);
        let alphabet = self.get_alphabet(self.datatype);
        self.write_csv_header(&mut writer, alphabet)?;
        stats.iter().for_each(|(site, chars)| {
            self.write_csv_content(&mut writer, site, chars, alphabet)
                .unwrap();
        });

        log::info!("\n{}", Yellow.paint("Output Files"));
        log::info!("{:18}: {}", "Alignment summary", self.output);

        Ok(())
    }

    pub fn write_summary_file(&mut self, site: &Sites, chars: &Chars) -> Result<()> {
        self.get_ouput_fname();
        let file = File::create(&self.output)
            .with_context(|| format!("Failed creating file {}", self.output))?;
        let mut writer = BufWriter::new(file);
        let alphabet = self.get_alphabet(self.datatype);
        self.write_csv_header(&mut writer, alphabet)?;
        self.write_csv_content(&mut writer, site, chars, alphabet)
            .unwrap();

        Ok(())
    }

    fn get_ouput_fname(&mut self) {
        self.output.push_str("_per_locus.csv")
    }

    fn write_csv_header<W: Write>(&self, writer: &mut W, alphabet: &str) -> Result<()> {
        write!(
            writer,
            "path,\
            locus,\
            ntaxa,\
            chars_count,\
            site_count,\
            conserved_sites,\
            proportion_cons_sites,\
            variable_sites,\
            proportion_var_sites,\
            parsimony_informative_sites,\
            proportion_pars_inf_sites,\
            missing_data,\
            proportion_missing_data,\
            gc_content,\
            at_content,\
        "
        )?;
        alphabet
            .chars()
            .for_each(|ch| write!(writer, "{},", ch).unwrap());
        writeln!(writer)?;
        Ok(())
    }

    fn write_csv_content<W: Write>(
        &self,
        writer: &mut W,
        site: &Sites,
        chars: &Chars,
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
            chars.total_chars
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
        write!(writer, "{},", chars.missing_data)?;
        write!(writer, "{},", chars.prop_missing_data)?;

        // GC content
        write!(
            writer,
            "{},",
            chars.gc_count as f64 / chars.total_chars as f64
        )?;

        // AT content
        write!(
            writer,
            "{},",
            chars.at_count as f64 / chars.total_chars as f64
        )?;

        // Characters
        alphabet.chars().for_each(|ch| {
            let count = chars.chars.get(&ch);
            match count {
                // Some(_) => (),
                Some(count) => write!(writer, "{},", count).unwrap(),
                None => write!(writer, "0,").unwrap(),
            }
        });
        writeln!(writer)?;

        writer.flush()?;
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
        log::info!("{}", Yellow.paint("General Summmary"));
        self.write_gen_sum();
        log::info!("{}", Yellow.paint("\nAlignment Summmary"));
        self.write_aln_sum();
        log::info!("{}", Yellow.paint("\nTaxon Summmary"));
        self.write_tax_sum();

        log::info!("{}", Yellow.paint("\nCharacter Count"));
        self.write_char_count();

        log::info!("{}", Yellow.paint("\nData Matrix Completeness"));
        self.write_matrix_comp();

        log::info!("{}", Yellow.paint("\nConserved Sequences"));
        self.write_cons_seq();

        log::info!("{}", Yellow.paint("\nVariable Sequences"));
        self.write_var_seq();

        log::info!("{}", Yellow.paint("\nParsimony Informative"));
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
            DataType::Aa => (),
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
            "{:18}: {}",
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
        log::info!("{:18}: {:.2} bp", "Mean length", &self.site.mean_sites);
    }

    fn write_tax_sum(&self) {
        log::info!("{:18}: {}", "Min taxa", utils::fmt_num(&self.chars.min_tax));
        log::info!("{:18}: {}", "Max taxa", utils::fmt_num(&self.chars.max_tax));
        log::info!("{:18}: {:.2}", "Mean taxa", self.chars.mean_tax);
    }

    fn write_char_count(&self) {
        let alphabet = self.get_alphabet(self.datatype);
        alphabet.chars().for_each(|ch| {
            if let Some(count) = self.chars.chars.get(&ch) {
                log::info!("{:18}: {}", ch, utils::fmt_num(count));
            }
        });
    }

    fn write_matrix_comp(&self) {
        self.complete
            .completeness
            .iter()
            .for_each(|(percent, ntax)| {
                let percent_str = format!("{}% taxa", percent);
                log::info!("{:18}: {}", percent_str, utils::fmt_num(ntax))
            });
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
        log::info!("{:18}: {:.2}", "Mean con. sites", &self.site.mean_cons_site);
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
        log::info!("{:18}: {:.2}", "Mean var. sites", &self.site.mean_var_site);
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
        log::info!("{:18}: {:.2}", "Mean inf. sites", &self.site.mean_inf_site);
    }
}

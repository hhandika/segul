use std::ffi::OsStr;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};

use crate::core::summary::{CharSummary, Chars, Completeness, SiteSummary, Sites};
use crate::helper::common::DataType;
use crate::helper::utils;

pub fn print_stats(site: &Sites, dna: &Chars) -> Result<()> {
    let io = io::stdout();
    let mut writer = BufWriter::new(io);

    writeln!(writer, "\x1b[0;33mAlignment\x1b[0m")?;
    writeln!(writer, "Taxa\t\t: {}", utils::fmt_num(&dna.ntax))?;
    writeln!(writer, "Length\t\t: {}\n", utils::fmt_num(&site.counts))?;

    writeln!(writer, "\x1b[0;33mSites\x1b[0m")?;
    writeln!(writer, "Conserved\t: {}", utils::fmt_num(&site.conserved))?;
    writeln!(writer, "Variable\t: {}", utils::fmt_num(&site.variable))?;
    writeln!(
        writer,
        "Parsimony inf.\t: {}\n",
        utils::fmt_num(&site.pars_inf)
    )?;
    writeln!(writer, "%Conserved\t: {:.2}%", site.prop_cons * 100.0)?;
    writeln!(writer, "%Variable\t: {:.2}%", site.prop_var * 100.0)?;
    writeln!(writer, "%Pars. inf.\t: {:.2}%\n", site.prop_var * 100.0)?;

    writeln!(writer, "\x1b[0;33mCharacters\x1b[0m")?;
    writeln!(writer, "Total\t\t: {}", utils::fmt_num(&dna.total_chars))?;
    writeln!(
        writer,
        "Missing data\t: {}",
        utils::fmt_num(&dna.missing_data)
    )?;
    writeln!(
        writer,
        "Prop. missing \t: {:.2}%",
        &dna.prop_missing_data * 100.0
    )?;

    dna.chars.iter().for_each(|(ch, count)| {
        writeln!(writer, "{}\t\t: {}", ch, utils::fmt_num(&count)).unwrap()
    });
    writer.flush()?;
    Ok(())
}

trait Alphabet {
    fn get_alphabet(&self, datatype: &DataType) -> &str {
        match datatype {
            DataType::Dna => "-?ACGTNRYSWKMBDHV.",
            DataType::Aa => "?-ARNDCQEGHILKMFPSTWYVYXBZJU.~*",
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

        Ok(())
    }

    pub fn write_summary_file(&mut self, site: &Sites, chars: &Chars) -> Result<()> {
        self.get_ouput_fname();
        let file = File::create(&self.output)
            .with_context(|| format!("Failed creating file {}", self.output))?;
        let mut writer = BufWriter::new(file);
        let alphabet = self.get_alphabet(&self.datatype);
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
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "\x1b[0;33mGeneral Summmary\x1b[0m")?;
        self.write_gen_sum(&mut writer)?;
        writeln!(writer, "\x1b[0;33mAlignment Summmary\x1b[0m")?;
        self.write_aln_sum(&mut writer)?;
        writeln!(writer, "\x1b[0;33mTaxon Summmary\x1b[0m")?;
        self.write_tax_sum(&mut writer)?;

        writeln!(writer, "\x1b[0;33mCharacter Count\x1b[0m")?;
        self.write_char_count(&mut writer)?;

        writeln!(writer, "\x1b[0;33mData Matrix Completeness\x1b[0m")?;
        self.write_matrix_comp(&mut writer)?;

        writeln!(writer, "\x1b[0;33mConserved Sequences\x1b[0m")?;
        self.write_cons_seq(&mut writer)?;

        writeln!(writer, "\x1b[0;33mVariable Sequences\x1b[0m")?;
        self.write_var_seq(&mut writer)?;

        writeln!(writer, "\x1b[0;33mParsimony Informative\x1b[0m")?;
        self.write_pars_inf(&mut writer)?;
        writeln!(writer)?;
        writer.flush()?;
        Ok(())
    }

    pub fn write_sum_to_file(&self, output: &str) -> Result<()> {
        let fname = self.get_output_fname(output);
        let file = File::create(&fname)
            .with_context(|| format!("Failed creating file {}", fname.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "General Summmary")?;
        self.write_gen_sum(&mut writer)?;
        writeln!(writer, "Alignment Summmary")?;
        self.write_aln_sum(&mut writer)?;
        writeln!(writer, "Taxon Summmary")?;
        self.write_tax_sum(&mut writer)?;

        writeln!(writer, "Character Count")?;
        self.write_char_count(&mut writer)?;

        writeln!(writer, "Data Matrix Completeness")?;
        self.write_matrix_comp(&mut writer)?;

        writeln!(writer, "Conserved Sequences")?;
        self.write_cons_seq(&mut writer)?;

        writeln!(writer, "Variable Sequences")?;
        self.write_var_seq(&mut writer)?;

        writeln!(writer, "Parsimony Informative")?;
        self.write_pars_inf(&mut writer)?;
        writeln!(writer)?;
        writer.flush()?;
        Ok(())
    }

    fn get_output_fname(&self, output: &str) -> PathBuf {
        let fname = format!("{}_summary.txt", output);
        PathBuf::from(fname)
    }

    fn write_gen_sum<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(
            writer,
            "Total taxa\t: {}",
            utils::fmt_num(&self.complete.total_tax)
        )?;
        writeln!(
            writer,
            "Total loci\t: {}",
            utils::fmt_num(&self.site.total_loci)
        )?;
        writeln!(
            writer,
            "Total sites\t: {}",
            utils::fmt_num(&self.site.total_sites)
        )?;
        writeln!(
            writer,
            "Missing data\t: {}",
            utils::fmt_num(&self.chars.missing_data)
        )?;
        writeln!(
            writer,
            "%Missing data\t: {:.2}%",
            &self.chars.prop_missing_data * 100.0
        )?;

        match self.datatype {
            DataType::Dna => self.write_dna_sum(writer)?,
            DataType::Aa => (),
            _ => panic!("Please specify datatype"),
        }
        writeln!(writer)?;

        Ok(())
    }

    fn write_dna_sum<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(writer, "GC content\t: {:.2}", self.chars.gc_content)?;
        writeln!(writer, "AT content\t: {:.2}", self.chars.at_content)?;
        writeln!(
            writer,
            "Characters\t: {}",
            utils::fmt_num(&self.chars.total_chars)
        )?;
        writeln!(
            writer,
            "Nucleotides\t: {}",
            utils::fmt_num(&self.chars.total_nucleotides)
        )?;
        Ok(())
    }

    fn write_aln_sum<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(
            writer,
            "Min length\t: {} bp",
            utils::fmt_num(&self.site.min_sites)
        )?;
        writeln!(
            writer,
            "Max length\t: {} bp",
            utils::fmt_num(&self.site.max_sites)
        )?;
        writeln!(writer, "Mean length\t: {:.2} bp\n", &self.site.mean_sites)?;
        Ok(())
    }

    fn write_tax_sum<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(
            writer,
            "Min taxa\t: {}",
            utils::fmt_num(&self.chars.min_tax)
        )?;
        writeln!(
            writer,
            "Max taxa\t: {}",
            utils::fmt_num(&self.chars.max_tax)
        )?;
        writeln!(writer, "Mean taxa\t: {:.2}\n", self.chars.mean_tax)?;

        Ok(())
    }

    fn write_char_count<W: Write>(&self, writer: &mut W) -> Result<()> {
        let alphabet = self.get_alphabet(self.datatype);
        alphabet.chars().for_each(|ch| {
            if let Some(count) = self.chars.chars.get(&ch) {
                writeln!(writer, "{}\t\t: {}", ch, utils::fmt_num(&count)).unwrap();
            }
        });
        writeln!(writer)?;
        Ok(())
    }

    fn write_matrix_comp<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.complete
            .completeness
            .iter()
            .for_each(|(percent, ntax)| {
                self.write_matrix_comp_content(writer, percent, ntax)
                    .expect("Failed printing data matrix completeness to stdout")
            });
        writeln!(writer)?;

        Ok(())
    }

    fn write_matrix_comp_content<W: Write>(
        &self,
        writer: &mut W,
        percent: &usize,
        ntax: &usize,
    ) -> Result<()> {
        if *percent < 10 {
            writeln!(writer, "{}% taxa\t\t: {}", percent, utils::fmt_num(ntax))?;
        } else {
            writeln!(writer, "{}% taxa\t: {}", percent, utils::fmt_num(ntax))?;
        }

        Ok(())
    }

    fn write_cons_seq<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(
            writer,
            "Con. loci\t: {}",
            utils::fmt_num(&self.site.cons_loci)
        )?;
        writeln!(
            writer,
            "%Con. loci\t: {:.2}%",
            self.site.prop_cons_loci * 100.0
        )?;
        writeln!(
            writer,
            "Con. sites\t: {}",
            utils::fmt_num(&self.site.total_cons_site)
        )?;
        writeln!(writer, "%Con. sites\t: {:.2}%", &self.site.prop_cons_site)?;
        writeln!(
            writer,
            "Min con. sites\t: {}",
            utils::fmt_num(&self.site.min_cons_site)
        )?;
        writeln!(
            writer,
            "Max con. sites\t: {}",
            utils::fmt_num(&self.site.max_cons_site)
        )?;
        writeln!(
            writer,
            "Mean con. sites\t: {:.2}\n",
            &self.site.mean_cons_site
        )?;

        Ok(())
    }

    fn write_var_seq<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(
            writer,
            "Var. loci\t: {}",
            utils::fmt_num(&self.site.var_loci)
        )?;
        writeln!(
            writer,
            "%Var. loci\t: {:.2}%",
            self.site.prop_var_loci * 100.0
        )?;
        writeln!(
            writer,
            "Var. sites\t: {}",
            utils::fmt_num(&self.site.total_var_site)
        )?;
        writeln!(writer, "%Var. sites\t: {:.2}%", &self.site.prop_var_site)?;
        writeln!(
            writer,
            "Min var. sites\t: {}",
            utils::fmt_num(&self.site.min_var_site)
        )?;
        writeln!(
            writer,
            "Max var. sites\t: {}",
            utils::fmt_num(&self.site.max_var_site)
        )?;
        writeln!(
            writer,
            "Mean var. sites\t: {:.2}\n",
            &self.site.mean_var_site
        )?;
        Ok(())
    }

    fn write_pars_inf<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(
            writer,
            "Inf. loci\t: {}",
            utils::fmt_num(&self.site.inf_loci)
        )?;
        writeln!(
            writer,
            "%Inf. loci\t: {:.2}%",
            self.site.prop_inf_loci * 100.0
        )?;
        writeln!(
            writer,
            "Inf. sites\t: {}",
            utils::fmt_num(&self.site.total_inf_site)
        )?;
        writeln!(writer, "%Inf. sites\t: {:.2}%", &self.site.prop_inf_site)?;
        writeln!(
            writer,
            "Min inf. sites\t: {}",
            utils::fmt_num(&self.site.min_inf_site)
        )?;
        writeln!(
            writer,
            "Max inf. sites\t: {}",
            utils::fmt_num(&self.site.max_inf_site)
        )?;
        writeln!(writer, "Mean inf. sites\t: {:.2}", &self.site.mean_inf_site)?;

        Ok(())
    }
}

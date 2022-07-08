use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::io::Write;
use std::path::{Path, PathBuf};

use ahash::AHashMap as HashMap;
use ansi_term::Colour::Yellow;
use anyhow::{Context, Result};
use indexmap::IndexSet;

use crate::helper::alphabet;
use crate::helper::stats::{CharSummary, Chars, Completeness, SiteSummary, Sites, Taxa};
use crate::helper::types::DataType;
use crate::helper::utils;
use crate::writer::FileWriter;

trait Alphabet {
    fn get_alphabet(&self, datatype: &DataType) -> &str {
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
    stats: &'a [(Sites, Chars, Taxa)],
}

impl<'a> CsvWriter<'a> {
    pub fn new(
        output: &'a Path,
        prefix: &'a Option<String>,
        datatype: &'a DataType,
        stats: &'a [(Sites, Chars, Taxa)],
    ) -> Self {
        Self {
            output,
            prefix,
            datatype,
            stats,
        }
    }

    // pub fn write_per_locus_summary(&self) -> Result<()> {
        // Ok(())
    // }

    pub fn write_taxon_summary(&self, ids: &IndexSet<String>) -> Result<()> {
        let default_prefix = "taxon_summary";
        let output = self.create_output_fnames(default_prefix);
        let mut writer = self.create_output_file(&output)?;
        write!(writer, "taxon,locus_counts")?;
        let alphabet = self.get_alphabet(self.datatype);
        alphabet
            .chars()
            .for_each(|ch| write!(writer, ",{}", ch).unwrap());
        writeln!(writer)?;
        let taxon_summary = self.summarize_taxa(ids);
        taxon_summary.iter().for_each(|(taxon, counts)| {
            write!(writer, "{},{}", taxon, counts.locus_counts)
                .expect("Failed taxon summary stats");
            alphabet.chars().for_each(|ch| {
                write!(writer, ",{}", counts.char_counts.get(&ch).unwrap_or(&0))
                    .expect("Failed taxon summary stats");
            });
            writeln!(writer).expect("Failed taxon summary stats");
        });
        Ok(())
    }

    pub fn write_locus_summary(&self) -> Result<()> {
        let default_prefix = "locus_summary";
        let output = self.create_output_fnames(default_prefix);
        let mut writer = self.create_output_file(&output)?;
        let alphabet = self.get_alphabet(self.datatype);
        self.write_locus_header(&mut writer, alphabet)?;
        self.stats.iter().for_each(|(site, chars, _)| {
            self.write_locus_content(&mut writer, site, chars, alphabet)
                .unwrap();
        });

        log::info!("{}", Yellow.paint("Output"));
        log::info!("{:18}: {}", "Output dir", self.output.display());

        Ok(())
    }

    fn create_output_fnames(&self, default_prefix: &str) -> PathBuf {
        match self.prefix {
            Some(fname) => {
                let out_name = format!("{}_{}", fname, default_prefix);
                self.output.join(out_name).with_extension("csv")
            },
            None => self.output.join(default_prefix).with_extension("csv"),
        }
    }

    fn summarize_taxa(&self, ids: &IndexSet<String>) -> BTreeMap<String, TaxonRecords> {
        let mut taxon_summary: BTreeMap<String, TaxonRecords> = BTreeMap::new();
        self.stats.iter().for_each(|(_, _, taxa)| {
            ids.iter().for_each(|id| {
                if let Some(char_counts) = taxa.records.get(id) {
                    match taxon_summary.get_mut(id) {
                        Some(taxon) => {
                            char_counts.iter().for_each(|(c, count)| {
                                *taxon.char_counts.entry(*c).or_insert(0) += count;
                            });
                            taxon.locus_counts += 1;
                        }
                        None => {
                            let mut taxon = TaxonRecords::new();
                            taxon.char_counts = char_counts.clone();
                            taxon.locus_counts = 1;
                            taxon_summary.insert(id.to_string(), taxon);
                        }
                    }
                }
            });
        });

        taxon_summary
    }

    fn write_locus_header<W: Write>(&self, writer: &mut W, alphabet: &str) -> Result<()> {
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
        self.write_alphabet_header(writer, alphabet)?;
        Ok(())
    }

    fn write_alphabet_header<W: Write>(&self, writer: &mut W, alphabet: &str) -> Result<()> {
        alphabet
            .chars()
            .for_each(|ch| write!(writer, "{},", ch).unwrap());
        writeln!(writer)?;
        Ok(())
    }

    fn write_locus_content<W: Write>(
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
            write!(writer, "{},", chars.chars.get(&ch).unwrap_or(&0)).unwrap();
        });
        writeln!(writer)?;
        writer.flush()?;
        Ok(())
    }
}

struct TaxonRecords {
    char_counts: HashMap<char, usize>,
    locus_counts: usize,
}

impl TaxonRecords {
    fn new() -> Self {
        Self {
            char_counts: HashMap::new(),
            locus_counts: 0,
        }
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
        log::info!("{}", Yellow.paint("Alignment Summmary"));
        self.write_aln_sum();
        log::info!("{}", Yellow.paint("Taxon Summmary"));
        self.write_tax_sum();

        log::info!("{}", Yellow.paint("Character Count"));
        self.write_char_count();

        log::info!("{}", Yellow.paint("Data Matrix Completeness"));
        self.write_matrix_comp();

        log::info!("{}", Yellow.paint("Conserved Sequences"));
        self.write_cons_seq();

        log::info!("{}", Yellow.paint("Variable Sequences"));
        self.write_var_seq();

        log::info!("{}", Yellow.paint("Parsimony Informative"));
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
        let alphabet = self.get_alphabet(self.datatype);
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

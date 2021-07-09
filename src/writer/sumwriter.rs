use std::fs::File;
use std::io::{self, BufWriter, Result, Write};
use std::path::PathBuf;

use crate::core::stats::{Completeness, Dna, DnaSummary, SiteSummary, Sites};
use crate::utils;

pub fn display_stats(site: &Sites, dna: &Dna) -> Result<()> {
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
    writeln!(writer, "Total\t: {}", utils::fmt_num(&dna.total_chars))?;
    writeln!(writer, "A\t: {}", utils::fmt_num(&dna.a_count))?;
    writeln!(writer, "C\t: {}", utils::fmt_num(&dna.c_count))?;
    writeln!(writer, "G\t: {}", utils::fmt_num(&dna.g_count))?;
    writeln!(writer, "T\t: {}", utils::fmt_num(&dna.t_count))?;
    writeln!(writer, "N\t: {}", utils::fmt_num(&dna.n_count))?;
    writeln!(writer, "?\t: {}", utils::fmt_num(&dna.missings))?;
    writeln!(writer, "-\t: {}", utils::fmt_num(&dna.gaps))?;
    writer.flush()?;
    Ok(())
}

pub struct CsvWriter {
    output: String,
}

impl CsvWriter {
    pub fn new(output: &str) -> Self {
        Self {
            output: String::from(output),
        }
    }

    pub fn write_summary_dir(&mut self, stats: &[(Sites, Dna)]) -> Result<()> {
        self.get_ouput_fname();
        let file = File::create(&self.output).expect("CANNOT WRITE THE STAT RESULTS");
        let mut writer = BufWriter::new(file);
        self.write_csv_header(&mut writer)?;
        stats.iter().for_each(|(site, dna)| {
            self.write_csv_content(&mut writer, site, dna).unwrap();
        });

        Ok(())
    }

    pub fn write_summary_file(&mut self, site: &Sites, dna: &Dna) -> Result<()> {
        self.get_ouput_fname();
        let file = File::create(&self.output).expect("CANNOT WRITE THE STAT RESULTS");
        let mut writer = BufWriter::new(file);
        self.write_csv_header(&mut writer)?;
        self.write_csv_content(&mut writer, site, dna).unwrap();

        Ok(())
    }

    fn get_ouput_fname(&mut self) {
        self.output.push_str("_per_locus.csv")
    }

    fn write_csv_header<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(
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
        AT_content,\
        GC_content,\
        A_counts,\
        T_counts,\
        G_counts,\
        C_counts,\
        gap_counts,\
        missing_counts,\
    "
        )?;
        Ok(())
    }

    fn write_csv_content<W: Write>(&self, writer: &mut W, site: &Sites, dna: &Dna) -> Result<()> {
        write!(
            writer,
            "{},{},{},{},",
            site.path.display(),
            site.path.file_stem().unwrap().to_string_lossy(),
            dna.ntax,
            dna.total_chars
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

        // GC content
        write!(
            writer,
            "{},",
            (dna.g_count as f64 + dna.c_count as f64) / dna.total_chars as f64
        )?;

        // AT content
        write!(
            writer,
            "{},",
            (dna.a_count as f64 + dna.t_count as f64) / dna.total_chars as f64
        )?;

        // Characters
        writeln!(
            writer,
            "{},{},{},{},{},{}",
            dna.a_count, dna.t_count, dna.g_count, dna.c_count, dna.gaps, dna.missings
        )?;

        writer.flush()?;
        Ok(())
    }
}

pub struct SummaryWriter<'s> {
    site: &'s SiteSummary,
    dna: &'s DnaSummary,
    complete: &'s Completeness,
}

impl<'s> SummaryWriter<'s> {
    pub fn new(site: &'s SiteSummary, dna: &'s DnaSummary, complete: &'s Completeness) -> Self {
        Self {
            site,
            dna,
            complete,
        }
    }

    pub fn display_summary(&self) -> Result<()> {
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

        writeln!(writer, "\x1b[0;33mTaxon Completeness\x1b[0m")?;
        self.write_tax_comp(&mut writer)?;

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
        let file = File::create(fname)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "General Summmary")?;
        self.write_gen_sum(&mut writer)?;
        writeln!(writer, "Alignment Summmary")?;
        self.write_aln_sum(&mut writer)?;
        writeln!(writer, "Taxon Summmary")?;
        self.write_tax_sum(&mut writer)?;

        writeln!(writer, "Character Count")?;
        self.write_char_count(&mut writer)?;

        writeln!(writer, "Taxon Completeness")?;
        self.write_tax_comp(&mut writer)?;

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
        writeln!(writer, "GC content\t: {:.2}", self.dna.gc_content)?;
        writeln!(writer, "AT content\t: {:.2}", self.dna.at_content)?;
        writeln!(
            writer,
            "Characters\t: {}",
            utils::fmt_num(&self.dna.total_chars)
        )?;
        writeln!(
            writer,
            "Nucleotides\t: {}",
            utils::fmt_num(&self.dna.total_nucleotides)
        )?;
        writeln!(
            writer,
            "Missing data\t: {}",
            utils::fmt_num(&self.dna.missing_data)
        )?;
        writeln!(
            writer,
            "%Missing data\t: {:.2}%\n",
            &self.dna.prop_missing_data * 100.0
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
        writeln!(writer, "Min taxa\t: {}", utils::fmt_num(&self.dna.min_tax))?;
        writeln!(writer, "Max taxa\t: {}", utils::fmt_num(&self.dna.max_tax))?;
        writeln!(writer, "Mean taxa\t: {:.2}\n", self.dna.mean_tax)?;

        Ok(())
    }

    fn write_char_count<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(writer, "A\t\t: {}", utils::fmt_num(&self.dna.total_a))?;
        writeln!(writer, "C\t\t: {}", utils::fmt_num(&self.dna.total_c))?;
        writeln!(writer, "G\t\t: {}", utils::fmt_num(&self.dna.total_g))?;
        writeln!(writer, "T\t\t: {}", utils::fmt_num(&self.dna.total_t))?;
        writeln!(writer, "N\t\t: {}", utils::fmt_num(&self.dna.total_n))?;
        writeln!(
            writer,
            "?\t\t: {}",
            utils::fmt_num(&self.dna.total_missings)
        )?;
        writeln!(writer, "-\t\t: {}", utils::fmt_num(&self.dna.total_gaps))?;
        writeln!(
            writer,
            "Undetermined\t: {}\n",
            utils::fmt_num(&self.dna.total_undetermined)
        )?;

        Ok(())
    }

    fn write_tax_comp<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(
            writer,
            "95% taxa\t: {}",
            utils::fmt_num(&self.complete.ntax_95)
        )?;
        writeln!(
            writer,
            "90% taxa\t: {}",
            utils::fmt_num(&self.complete.ntax_90)
        )?;
        writeln!(
            writer,
            "85% taxa\t: {}",
            utils::fmt_num(&self.complete.ntax_85)
        )?;
        writeln!(
            writer,
            "80% taxa\t: {}",
            utils::fmt_num(&self.complete.ntax_80)
        )?;
        writeln!(
            writer,
            "75% taxa\t: {}",
            utils::fmt_num(&self.complete.ntax_75)
        )?;
        writeln!(
            writer,
            "70% taxa\t: {}",
            utils::fmt_num(&self.complete.ntax_70)
        )?;
        writeln!(
            writer,
            "65% taxa\t: {}",
            utils::fmt_num(&self.complete.ntax_65)
        )?;
        writeln!(
            writer,
            "60% taxa\t: {}",
            utils::fmt_num(&self.complete.ntax_60)
        )?;
        writeln!(
            writer,
            "55% taxa\t: {}",
            utils::fmt_num(&self.complete.ntax_55)
        )?;
        writeln!(
            writer,
            "50% taxa\t: {}\n",
            utils::fmt_num(&self.complete.ntax_50)
        )?;

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

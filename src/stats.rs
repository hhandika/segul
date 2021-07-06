//! A module for sequence statistics.
use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::File;
use std::io::{self, BufWriter, Result, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc::channel;

use indexmap::IndexMap;
use rayon::prelude::*;

use crate::alignment::Alignment;
use crate::common::SeqFormat;
use crate::finder::IDs;
use crate::utils;

pub struct SeqStats<'a> {
    input_format: &'a SeqFormat,
    output: &'a str,
    ntax: usize,
}

impl<'a> SeqStats<'a> {
    pub fn new(input_format: &'a SeqFormat, output: &'a str) -> Self {
        Self {
            input_format,
            output,
            ntax: 0,
        }
    }

    pub fn get_seq_stats_file(&self, path: &Path) {
        let spin = utils::set_spinner();
        spin.set_message("Getting alignments...");
        let (site, dna) = self.get_stats(path);
        spin.finish_with_message("DONE!\n");
        self.display_stats(&site, &dna).unwrap();
    }

    pub fn get_stats_dir(&mut self, files: &[PathBuf]) {
        let spin = utils::set_spinner();
        spin.set_message("Indexing alignments...");
        self.get_ntax(files);
        spin.set_message("Getting summary stats...");
        let mut stats: Vec<(Sites, Dna)> = self.par_get_stats(files);
        stats.sort_by(|a, b| alphanumeric_sort::compare_path(&a.0.path, &b.0.path));
        let (sites, dna, complete) = self.get_summary_dna(&stats);
        spin.set_message("Writing results...");
        CsvWriter::new(self.output)
            .write_locus_summary(&stats)
            .expect("CANNOT WRITE PER LOCUS SUMMARY STATS");
        let sum = SummaryWriter::new(&sites, &dna, &complete);
        sum.write_sum_to_file(Path::new("SEGUL-summary"))
            .expect("CANNOT CREATE FILE FOR SUMMARY OUPUT");
        spin.finish_with_message("DONE!\n");
        sum.display_summary()
            .expect("CANNOT WRITE SUMMARY TO STDOUT");
    }

    fn get_ntax(&mut self, files: &[PathBuf]) {
        self.ntax = IDs::new(files, self.input_format).get_id_all().len();
    }

    fn par_get_stats(&self, files: &[PathBuf]) -> Vec<(Sites, Dna)> {
        let (send, rec) = channel();
        files.par_iter().for_each_with(send, |s, file| {
            s.send(self.get_stats(file)).unwrap();
        });
        rec.iter().collect()
    }

    fn get_stats(&self, path: &Path) -> (Sites, Dna) {
        let mut aln = Alignment::new();
        aln.get_aln_any(path, self.input_format);
        let mut dna = Dna::new();
        dna.count_chars(&aln);
        let mut sites = Sites::new(path);
        sites.get_stats(&aln.alignment);

        (sites, dna)
    }

    fn get_summary_dna(&self, stats: &[(Sites, Dna)]) -> (SiteSummary, DnaSummary, Completeness) {
        let (sites, dna): (Vec<Sites>, Vec<Dna>) =
            stats.par_iter().map(|p| (p.0.clone(), p.1.clone())).unzip();
        let mut sum_sites = SiteSummary::new();
        sum_sites.get_summary(&sites);
        let mut sum_dna = DnaSummary::new();
        sum_dna.get_summary(&dna);
        let mut ntax_comp = Completeness::new(&self.ntax);
        ntax_comp.get_ntax_completeness(&dna);
        (sum_sites, sum_dna, ntax_comp)
    }

    fn display_stats(&self, site: &Sites, dna: &Dna) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);

        writeln!(writer, "\x1b[0;33mAlignment\x1b[0m")?;
        writeln!(writer, "Taxa\t\t: {}", utils::fmt_num(&dna.ntax))?;
        writeln!(writer, "Length\t\t: {}\n", utils::fmt_num(&dna.total_chars))?;

        writeln!(writer, "\x1b[0;33mSites\x1b[0m")?;
        writeln!(writer, "Count\t\t: {}", utils::fmt_num(&site.counts))?;
        writeln!(writer, "Conserved\t: {}", utils::fmt_num(&site.conserved))?;
        writeln!(writer, "Variable\t: {}", utils::fmt_num(&site.variable))?;
        writeln!(
            writer,
            "Parsimony inf.\t: {}\n",
            utils::fmt_num(&site.pars_inf)
        )?;
        writeln!(writer, "Prop. conserved\t: {:.2}%", site.prop_cons * 100.0)?;
        writeln!(writer, "Prop. variable\t: {:.2}%", site.prop_var * 100.0)?;
        writeln!(writer, "Prop. p. inf.\t: {:.2}%\n", site.prop_var * 100.0)?;

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
}

struct CsvWriter {
    output: String,
}

impl CsvWriter {
    fn new(output: &str) -> Self {
        Self {
            output: String::from(output),
        }
    }

    fn write_locus_summary(&mut self, stats: &[(Sites, Dna)]) -> Result<()> {
        self.get_ouput_fname();
        let file = File::create(&self.output).expect("CANNOT WRITE THE STAT RESULTS");
        let mut writer = BufWriter::new(file);
        self.write_csv_header(&mut writer)?;
        stats.iter().for_each(|(site, dna)| {
            self.write_csv_content(&mut writer, site, dna).unwrap();
        });

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

struct SummaryWriter<'s> {
    site: &'s SiteSummary,
    dna: &'s DnaSummary,
    complete: &'s Completeness,
}

impl<'s> SummaryWriter<'s> {
    fn new(site: &'s SiteSummary, dna: &'s DnaSummary, complete: &'s Completeness) -> Self {
        Self {
            site,
            dna,
            complete,
        }
    }

    fn display_summary(&self) -> Result<()> {
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

    fn write_sum_to_file(&self, output: &Path) -> Result<()> {
        let fname = output.with_extension("txt");
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

struct SiteSummary {
    // General site summary
    total_loci: usize,
    total_sites: usize,
    min_sites: usize,
    max_sites: usize,
    mean_sites: f64,

    // Conserved site summary
    cons_loci: usize,
    prop_cons_loci: f64,
    total_cons_site: usize,
    min_cons_site: usize,
    max_cons_site: usize,
    mean_cons_site: f64,

    // Variable site summary
    var_loci: usize,
    prop_var_loci: f64,
    total_var_site: usize,
    min_var_site: usize,
    max_var_site: usize,
    mean_var_site: f64,

    // Parsimony inf site summary
    inf_loci: usize,
    prop_inf_loci: f64,
    total_inf_site: usize,
    min_inf_site: usize,
    max_inf_site: usize,
    mean_inf_site: f64,
}

impl SiteSummary {
    fn new() -> Self {
        Self {
            total_sites: 0,
            total_loci: 0,
            min_sites: 0,
            max_sites: 0,
            mean_sites: 0.0,
            var_loci: 0,
            prop_var_loci: 0.0,
            total_var_site: 0,
            cons_loci: 0,
            prop_cons_loci: 0.0,
            total_cons_site: 0,
            min_cons_site: 0,
            max_cons_site: 0,
            mean_cons_site: 0.0,
            min_var_site: 0,
            max_var_site: 0,
            mean_var_site: 0.0,
            inf_loci: 0,
            prop_inf_loci: 0.0,
            total_inf_site: 0,
            min_inf_site: 0,
            max_inf_site: 0,
            mean_inf_site: 0.0,
        }
    }

    fn get_summary(&mut self, sites: &[Sites]) {
        self.total_loci = sites.len();
        self.total_sites = sites.iter().map(|s| s.counts).sum();
        self.min_sites = sites.iter().map(|s| s.counts).min().unwrap();
        self.max_sites = sites.iter().map(|s| s.counts).max().unwrap();
        self.mean_sites = self.total_sites as f64 / self.total_loci as f64;
        self.count_cons_sites(&sites);
        self.count_var_sites(&sites);
        self.count_inf_sites(&sites);
    }

    fn count_cons_sites(&mut self, sites: &[Sites]) {
        self.cons_loci = sites.iter().filter(|s| s.variable == 0).count();
        self.prop_cons_loci = self.cons_loci as f64 / self.total_loci as f64;
        self.total_cons_site = sites.iter().map(|s| s.conserved).sum();
        self.min_cons_site = sites.iter().map(|s| s.conserved).min().unwrap();
        self.max_cons_site = sites.iter().map(|s| s.conserved).max().unwrap();
        self.mean_cons_site = self.total_cons_site as f64 / self.total_sites as f64;
    }

    fn count_var_sites(&mut self, sites: &[Sites]) {
        self.var_loci = sites.iter().filter(|s| s.variable > 0).count();
        self.prop_var_loci = self.var_loci as f64 / self.total_loci as f64;
        self.total_var_site = sites.iter().map(|s| s.variable).sum();
        self.min_var_site = sites.iter().map(|s| s.variable).min().unwrap();
        self.max_var_site = sites.iter().map(|s| s.variable).max().unwrap();
        self.mean_var_site = self.total_var_site as f64 / self.total_sites as f64;
    }

    fn count_inf_sites(&mut self, sites: &[Sites]) {
        self.inf_loci = sites.iter().filter(|s| s.pars_inf > 0).count();
        self.prop_inf_loci = self.inf_loci as f64 / self.total_loci as f64;
        self.total_inf_site = sites.iter().map(|s| s.pars_inf).sum();
        self.min_inf_site = sites.iter().map(|s| s.pars_inf).min().unwrap();
        self.max_inf_site = sites.iter().map(|s| s.pars_inf).max().unwrap();
        self.mean_inf_site = self.total_inf_site as f64 / self.total_sites as f64;
    }
}

struct DnaSummary {
    min_tax: usize,
    max_tax: usize,
    mean_tax: f64,
    gc_content: f64,
    at_content: f64,
    missing_data: usize,
    prop_missing_data: f64,
    total_chars: usize,
    total_nucleotides: usize,
    total_a: usize,
    total_c: usize,
    total_g: usize,
    total_t: usize,
    total_n: usize,
    total_missings: usize,
    total_gaps: usize,
    total_undetermined: usize,
}

impl DnaSummary {
    fn new() -> Self {
        Self {
            total_chars: 0,
            min_tax: 0,
            max_tax: 0,
            mean_tax: 0.0,
            gc_content: 0.0,
            at_content: 0.0,
            missing_data: 0,
            prop_missing_data: 0.0,
            total_nucleotides: 0,
            total_a: 0,
            total_c: 0,
            total_g: 0,
            total_t: 0,
            total_n: 0,
            total_missings: 0,
            total_gaps: 0,
            total_undetermined: 0,
        }
    }

    fn get_summary(&mut self, dna: &[Dna]) {
        self.min_tax = dna.iter().map(|d| d.ntax).min().unwrap();
        self.max_tax = dna.iter().map(|d| d.ntax).max().unwrap();
        let sum_tax: usize = dna.iter().map(|d| d.ntax).sum();
        self.mean_tax = sum_tax as f64 / dna.len() as f64;
        self.total_chars = dna.iter().map(|d| d.total_chars).sum();
        self.count_chars(dna);
        self.get_total_nucleotides();
        self.count_gc_at_content();
        self.count_missing_data();
    }

    fn count_chars(&mut self, dna: &[Dna]) {
        self.total_a = dna.iter().map(|d| d.a_count).sum();
        self.total_t = dna.iter().map(|d| d.t_count).sum();
        self.total_g = dna.iter().map(|d| d.g_count).sum();
        self.total_c = dna.iter().map(|d| d.c_count).sum();
        self.total_n = dna.iter().map(|d| d.n_count).sum();
        self.total_missings = dna.iter().map(|d| d.missings).sum();
        self.total_gaps = dna.iter().map(|d| d.gaps).sum();
        self.total_undetermined = dna.iter().map(|d| d.undetermined).sum();
    }

    fn get_total_nucleotides(&mut self) {
        self.total_nucleotides = self.total_a + self.total_t + self.total_g + self.total_c
    }

    fn count_gc_at_content(&mut self) {
        self.gc_content = (self.total_g + self.total_c) as f64 / self.total_chars as f64;
        self.at_content = (self.total_g + self.total_c) as f64 / self.total_chars as f64;
    }

    fn count_missing_data(&mut self) {
        self.missing_data = self.total_missings + self.total_gaps + self.total_n;
        self.prop_missing_data = self.missing_data as f64 / self.total_chars as f64;
    }
}

// #[derive(Debug, Send, Sync)]
struct Completeness {
    ntax_95: usize,
    ntax_90: usize,
    ntax_85: usize,
    ntax_80: usize,
    ntax_75: usize,
    ntax_70: usize,
    ntax_65: usize,
    ntax_60: usize,
    ntax_55: usize,
    ntax_50: usize,
    total_tax: usize,
}

impl Completeness {
    fn new(total_tax: &usize) -> Self {
        Self {
            ntax_95: 0,
            ntax_90: 0,
            ntax_85: 0,
            ntax_80: 0,
            ntax_75: 0,
            ntax_70: 0,
            ntax_65: 0,
            ntax_60: 0,
            ntax_55: 0,
            ntax_50: 0,
            total_tax: *total_tax,
        }
    }

    fn get_ntax_completeness(&mut self, dna: &[Dna]) {
        let ntax: Vec<usize> = dna.iter().map(|d| d.ntax).collect();
        self.ntax_95 = self.count_min_tax(&ntax, 0.95);
        self.ntax_90 = self.count_min_tax(&ntax, 0.9);
        self.ntax_85 = self.count_min_tax(&ntax, 0.85);
        self.ntax_80 = self.count_min_tax(&ntax, 0.8);
        self.ntax_75 = self.count_min_tax(&ntax, 0.75);
        self.ntax_70 = self.count_min_tax(&ntax, 0.7);
        self.ntax_65 = self.count_min_tax(&ntax, 0.65);
        self.ntax_60 = self.count_min_tax(&ntax, 0.6);
        self.ntax_55 = self.count_min_tax(&ntax, 0.55);
        self.ntax_50 = self.count_min_tax(&ntax, 0.5);
    }

    fn count_min_tax(&self, ntax: &[usize], percent: f64) -> usize {
        ntax.iter()
            .filter(|&n| n >= &self.compute_min_taxa(percent))
            .count()
    }

    fn compute_min_taxa(&self, percent: f64) -> usize {
        (self.total_tax as f64 * percent).floor() as usize
    }
}

#[derive(Debug, Clone)]
struct Sites {
    path: PathBuf,
    conserved: usize,
    variable: usize,
    pars_inf: usize,
    counts: usize,
    prop_var: f64,
    prop_cons: f64,
    prop_pinf: f64,
}

impl Sites {
    pub fn new(path: &Path) -> Self {
        Self {
            path: PathBuf::from(path),
            conserved: 0,
            variable: 0,
            pars_inf: 0,
            counts: 0,
            prop_var: 0.0,
            prop_cons: 0.0,
            prop_pinf: 0.0,
        }
    }

    fn get_stats(&mut self, matrix: &IndexMap<String, String>) {
        let site_matrix = self.index_sites(matrix);
        self.get_site_stats(&site_matrix);
        self.count_sites();
        self.get_proportion();
    }

    fn index_sites(&mut self, matrix: &IndexMap<String, String>) -> HashMap<usize, Vec<u8>> {
        let mut site_matrix: HashMap<usize, Vec<u8>> = HashMap::new();
        matrix.values().for_each(|seq| {
            seq.bytes().enumerate().for_each(|(idx, dna)| {
                match site_matrix.get_mut(&idx) {
                    Some(value) => match dna {
                        b'a' | b'g' | b't' | b'c' | b'A' | b'G' | b'T' | b'C' => value.push(dna),
                        _ => (), // ignore ambigous characters
                    },
                    None => match dna {
                        b'a' | b'g' | b't' | b'c' | b'A' | b'G' | b'T' | b'C' => {
                            site_matrix.insert(idx, vec![dna]);
                        }
                        _ => (),
                    },
                }
            })
        });

        site_matrix
    }

    fn get_site_stats(&mut self, site_matrix: &HashMap<usize, Vec<u8>>) {
        site_matrix.values().for_each(|site| {
            let n_patterns = self.get_patterns(site);
            if n_patterns >= 2 {
                self.pars_inf += 1;
            }
        });
    }

    fn get_patterns(&mut self, site: &[u8]) -> usize {
        let mut uniques: Vec<u8> = site.to_vec();
        uniques.sort_unstable();
        uniques.dedup();

        // We consider variable sites
        // when the characters not all the same
        let mut n_patterns = 0;
        if uniques.len() >= 2 {
            self.variable += 1;
            uniques.iter().for_each(|ch| {
                let patterns = site.iter().filter(|&site| site == ch).count();
                if patterns >= 2 {
                    n_patterns += 1;
                }
            });
        } else {
            self.conserved += 1;
        }

        n_patterns
    }

    fn count_sites(&mut self) {
        self.counts = self.conserved + self.variable;
    }

    fn get_proportion(&mut self) {
        self.prop_cons = self.conserved as f64 / self.counts as f64;
        self.prop_var = self.variable as f64 / self.counts as f64;
        self.prop_pinf = self.pars_inf as f64 / self.counts as f64;
    }
}

#[derive(Debug, Clone)]
struct Dna {
    a_count: usize,
    c_count: usize,
    g_count: usize,
    t_count: usize,
    n_count: usize,
    missings: usize,
    gaps: usize,
    undetermined: usize, // alignment length
    total_chars: usize,  // All characters count
    ntax: usize,
}

impl Dna {
    fn new() -> Self {
        Self {
            a_count: 0,
            c_count: 0,
            g_count: 0,
            t_count: 0,
            n_count: 0,
            missings: 0,
            gaps: 0,
            undetermined: 0,
            total_chars: 0,
            ntax: 0,
        }
    }

    fn count_chars(&mut self, aln: &Alignment) {
        self.ntax = aln.header.ntax;
        self.total_chars = aln.header.nchar * self.ntax;
        aln.alignment.values().for_each(|seqs| {
            seqs.bytes().for_each(|ch| match ch {
                b'a' | b'A' => self.a_count += 1,
                b'c' | b'C' => self.c_count += 1,
                b'g' | b'G' => self.g_count += 1,
                b't' | b'T' => self.t_count += 1,
                b'n' | b'N' => self.n_count += 1,
                b'?' | b'.' | b'~' => self.missings += 1,
                b'O' | b'o' | b'X' | b'x' => self.missings += 1, // Following iqtree treatments
                b'-' => self.gaps += 1,
                _ => self.undetermined += 1,
            })
        })
    }
}

#[cfg(test)]
mod test {
    // use indexmap::IndexMap;
    use super::*;

    fn get_matrix(id: &[&str], seq: &[&str]) -> IndexMap<String, String> {
        let mut matrix = IndexMap::new();
        id.iter().zip(seq.iter()).for_each(|(i, s)| {
            matrix.insert(i.to_string(), s.to_string());
        });

        matrix
    }

    #[test]
    fn pattern_count_test() {
        let site = b"AATT";
        let site_2 = b"AATTGG";
        let pattern = Sites::new(Path::new(".")).get_patterns(site);
        let pattern_2 = Sites::new(Path::new(".")).get_patterns(site_2);
        assert_eq!(2, pattern);
        assert_eq!(3, pattern_2);
    }

    #[test]
    fn count_parsimony_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::new(Path::new("."));
        let smat = site.index_sites(&mat);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
    }

    #[test]
    fn count_variable_sites_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT", "ATTA", "ATGC", "ATGA"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::new(Path::new("."));
        let smat = site.index_sites(&mat);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
    }

    #[test]
    fn count_parsimony_gap_test() {
        let id = ["ABC", "ABE", "ABF", "ABD"];
        let seq = ["AATT---", "ATTA---", "ATGC---", "ATGA---"];
        let mat = get_matrix(&id, &seq);
        let mut site = Sites::new(Path::new("."));
        let smat = site.index_sites(&mat);
        site.get_site_stats(&smat);
        assert_eq!(1, site.pars_inf);
        assert_eq!(3, site.variable);
    }

    #[test]
    fn get_site_stats_test() {
        let path = Path::new("test_files/concat.fasta");
        let input_format = SeqFormat::Fasta;
        let mut aln = Alignment::new();
        aln.get_aln_any(path, &input_format);
        let mut site = Sites::new(Path::new("."));
        let smat = site.index_sites(&aln.alignment);
        site.get_site_stats(&smat);
        assert_eq!(18, site.conserved);
        assert_eq!(8, site.variable);
        assert_eq!(2, site.pars_inf);
    }

    #[test]
    fn filter_min_tax_test() {
        let ntax = vec![10, 8, 20, 30, 60];
        let comp = Completeness::new(&60);
        assert_eq!(2, comp.count_min_tax(&ntax, 0.5))
    }
}

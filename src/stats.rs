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
use crate::common::{Header, SeqFormat};
// use crate::finder::IDs;
use crate::utils;

pub fn get_seq_stats(path: &Path, input_format: &SeqFormat) {
    let spin = utils::set_spinner();
    spin.set_message("Getting alignments...");
    let mut aln = Alignment::new();
    aln.get_aln_any(path, input_format);
    spin.set_message("Counting characters...");
    let mut dna = Dna::new();
    dna.count_chars(&aln);
    spin.set_message("Getting summary stats...");
    let mut sites = Sites::new(path);
    sites.get_stats(&aln.alignment);
    spin.finish_with_message("DONE!\n");
    display_stats(&sites, &dna, &aln.header).unwrap();
}

pub fn get_stats_dir(files: &[PathBuf], input_format: &SeqFormat) {
    let spin = utils::set_spinner();
    spin.set_message("Processing alignments...");
    let (send, rec) = channel();
    files.par_iter().for_each_with(send, |s, file| {
        s.send(get_stats(file, input_format)).unwrap();
    });

    let mut stats: Vec<(Sites, Dna)> = rec.iter().collect();

    stats.sort_by(|a, b| alphanumeric_sort::compare_path(&a.0.path, &b.0.path));
    spin.set_message("Getting summary stats...");
    let (sites, dna) = get_summary_dna(&stats);
    spin.set_message("Writing results...");
    write_aln_stats(&stats).unwrap();
    spin.finish_with_message("DONE!\n");
    display_summary(&sites, &dna).unwrap();
}

fn get_stats(path: &Path, input_format: &SeqFormat) -> (Sites, Dna) {
    let mut aln = Alignment::new();
    aln.get_aln_any(path, input_format);
    let mut dna = Dna::new();
    dna.count_chars(&aln);
    let mut sites = Sites::new(path);
    sites.get_stats(&aln.alignment);

    (sites, dna)
}

fn get_summary_dna(stats: &[(Sites, Dna)]) -> (SiteSummary, DnaSummary) {
    let (sites, dna): (Vec<Sites>, Vec<Dna>) =
        stats.par_iter().map(|p| (p.0.clone(), p.1.clone())).unzip();
    let mut sum_sites = SiteSummary::new();
    sum_sites.get_summary(&sites);
    let mut sum_dna = DnaSummary::new();
    sum_dna.get_summary(&dna);

    (sum_sites, sum_dna)
}

fn display_summary(site: &SiteSummary, dna: &DnaSummary) -> Result<()> {
    let io = io::stdout();
    let mut writer = BufWriter::new(io);
    writeln!(writer, "Total loci\t: {}", utils::fmt_num(&site.total_loci))?;
    writeln!(
        writer,
        "Pars. inf. loci\t: {}",
        utils::fmt_num(&site.inf_loci)
    )?;
    writeln!(
        writer,
        "Total sites\t: {}",
        utils::fmt_num(&site.total_sites)
    )?;
    writeln!(
        writer,
        "Total chars\t: {}",
        utils::fmt_num(&dna.total_chars)
    )?;
    writeln!(writer, "Min taxa\t: {}", utils::fmt_num(&dna.min_tax))?;
    writeln!(writer, "Max taxa\t: {}", utils::fmt_num(&dna.max_tax))?;
    writeln!(writer, "Mean taxa\t: {:.2}", dna.mean_tax)?;
    writeln!(writer)?;
    writer.flush()?;
    Ok(())
}

fn write_aln_stats(stats: &[(Sites, Dna)]) -> Result<()> {
    let fname = "SEGUL-stats_per_locus.csv";
    let file = File::create(fname).expect("CANNOT WRITE THE STAT RESULTS");
    let mut writer = BufWriter::new(file);
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
    stats.iter().for_each(|(site, dna)| {
        write_content(&mut writer, site, dna).unwrap();
    });

    Ok(())
}

fn write_content<W: Write>(writer: &mut W, site: &Sites, dna: &Dna) -> Result<()> {
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

fn display_stats(site: &Sites, dna: &Dna, aln: &Header) -> Result<()> {
    let io = io::stdout();
    let mut writer = BufWriter::new(io);

    writeln!(writer, "\x1b[0;33mAlignment\x1b[0m")?;
    writeln!(writer, "Taxa\t\t: {}", utils::fmt_num(&aln.ntax))?;
    writeln!(writer, "Length\t\t: {}\n", utils::fmt_num(&aln.nchar))?;

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

struct SiteSummary {
    total_loci: usize,
    total_sites: usize,
    mean_sites: usize,
    inf_loci: usize,
}

impl SiteSummary {
    fn new() -> Self {
        Self {
            total_sites: 0,
            total_loci: 0,
            mean_sites: 0,
            inf_loci: 0,
        }
    }

    fn get_summary(&mut self, sites: &[Sites]) {
        self.total_loci = sites.len();
        self.total_sites = sites.iter().map(|s| s.counts).sum();
        self.inf_loci = sites.iter().filter(|s| s.pars_inf > 0).count();
        self.mean_sites = self.total_sites / self.total_loci;
    }
}

struct DnaSummary {
    min_tax: usize,
    max_tax: usize,
    mean_tax: f64,
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
                b'?' | b'.' => self.missings += 1,
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
}

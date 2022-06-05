use std::path::{Path, PathBuf};

use ansi_term::Colour::Yellow;
use clap::ArgMatches;

use crate::cli::{ConcatCli, InputCli, InputPrint, OutputCli};
use crate::handler::filter::{Params, SeqFilter};
use crate::helper::finder::IDs;
use crate::helper::types::{DataType, InputFmt, OutputFmt, PartitionFmt};
use crate::helper::{filenames, utils};
use crate::parser::txt;

impl InputCli for FilterParser<'_> {}
impl InputPrint for FilterParser<'_> {}
impl OutputCli for FilterParser<'_> {}
impl ConcatCli for FilterParser<'_> {}

pub(in crate::cli) struct FilterParser<'a> {
    matches: &'a ArgMatches,
    input_fmt: InputFmt,
    input_dir: Option<PathBuf>,
    output_dir: PathBuf,
    files: Vec<PathBuf>,
    params: Params,
    ntax: usize,
    percent: f64,
    datatype: DataType,
}

impl<'a> FilterParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches) -> Self {
        Self {
            matches,
            input_fmt: InputFmt::Fasta,
            input_dir: None,
            output_dir: PathBuf::new(),
            files: Vec::new(),
            params: Params::MinTax(0),
            ntax: 0,
            percent: 0.0,
            datatype: DataType::Dna,
        }
    }

    pub(in crate::cli) fn filter(&mut self) {
        self.input_fmt = self.parse_input_fmt(self.matches);
        self.datatype = self.parse_datatype(self.matches);
        let task_desc = "Alignment filtering";
        self.files = if self.matches.is_present("dir") {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &self.input_fmt)
        } else {
            self.parse_input(self.matches)
        };
        self.print_input(
            &self.input_dir,
            task_desc,
            self.files.len(),
            &self.input_fmt,
            &self.datatype,
        );

        if self.is_npercent() {
            self.get_min_taxa_npercent();
        } else {
            self.get_params();
            self.set_output_path();
            self.filter_aln();
        }
    }

    fn get_min_taxa_npercent(&mut self) {
        let npercent = self.parse_npercent();
        self.count_ntax();
        npercent.iter().for_each(|&np| {
            self.percent = np;
            let min_taxa = self.count_min_tax();
            self.params = Params::MinTax(min_taxa);
            self.set_multi_output_path();
            self.filter_aln();
            utils::print_divider();
        });
    }

    fn filter_aln(&self) {
        let is_overwrite = self.parse_overwrite_opts(self.matches);
        self.check_output_dir_exist(&self.output_dir, is_overwrite);
        self.print_params();
        let mut filter = SeqFilter::new(
            &self.files,
            &self.input_fmt,
            &self.datatype,
            &self.output_dir,
            &self.params,
        );
        match self.check_concat() {
            Some(part_fmt) => {
                let output_fmt = if self.matches.is_present("output-format") {
                    self.parse_output_fmt(self.matches)
                } else {
                    OutputFmt::Nexus
                };
                let prefix = self.parse_prefix(self.matches, &self.output_dir);
                let output = filenames::create_output_fname(&self.output_dir, &prefix, &output_fmt);
                filter.set_concat(&output, &output_fmt, &part_fmt);
                filter.filter_aln();
            }
            None => filter.filter_aln(),
        }
    }

    fn get_params(&mut self) {
        self.params = match self.matches {
            m if m.is_present("percent") => {
                self.percent = self.get_percent();
                self.count_ntax();
                let min_taxa = self.count_min_tax();
                Params::MinTax(min_taxa)
            }
            m if m.is_present("aln-len") => Params::AlnLen(self.parse_aln_len()),
            m if m.is_present("pars-inf") => Params::ParsInf(self.parse_pars_inf()),
            m if m.is_present("percent-inf") => Params::PercInf(self.count_percent_inf()),
            m if m.is_present("taxon-all") => Params::TaxonAll(self.parse_taxon_id()),
            _ => unreachable!("Invalid parameters!"),
        }
    }

    fn parse_taxon_id(&self) -> Vec<String> {
        let id_path = Path::new(
            self.matches
                .value_of("taxon-id")
                .expect("Failed to parse taxon-id"),
        );
        txt::parse_text_file(id_path)
    }

    fn count_percent_inf(&self) -> f64 {
        let perc_inf = self
            .matches
            .value_of("percent-inf")
            .expect("Failed parsing percent informative values");
        perc_inf
            .parse::<f64>()
            .expect("Failed parsing percent inf to floating points")
    }

    fn parse_aln_len(&self) -> usize {
        let len = self
            .matches
            .value_of("aln-len")
            .expect("Failed parsing an alignment length value");
        len.parse::<usize>()
            .expect("Failed parsing an alignment value to integer")
    }

    fn parse_pars_inf(&self) -> usize {
        let len = self
            .matches
            .value_of("pars-inf")
            .expect("Failed parsing a parsimony informative value");
        len.parse::<usize>()
            .expect("Failed parsing a parsimony informative value to in integer")
    }

    fn count_ntax(&mut self) {
        if self.matches.is_present("ntax") {
            self.ntax = self.parse_ntax();
        } else {
            let spin = utils::set_spinner();
            spin.set_message("Counting the number of taxa...");
            self.ntax = IDs::new(&self.files, &self.input_fmt, &self.datatype)
                .id_unique()
                .len();
            spin.finish_with_message("Finished counting the number of taxa!\n");
        };
    }

    fn count_min_tax(&self) -> usize {
        (self.ntax as f64 * self.percent).floor() as usize
    }

    fn parse_npercent(&self) -> Vec<f64> {
        self.matches
            .values_of("npercent")
            .expect("Failed parsing npercent")
            .map(|np| self.parse_percent(np))
            .collect()
    }

    fn get_percent(&mut self) -> f64 {
        let percent = self
            .matches
            .value_of("percent")
            .expect("Failed parsing a percentage value");
        self.parse_percent(percent)
    }

    fn parse_percent(&self, percent: &str) -> f64 {
        percent
            .parse::<f64>()
            .expect("Failed parsing a percentage value to a floating point number")
    }

    fn parse_ntax(&self) -> usize {
        let ntax = self
            .matches
            .value_of("ntax")
            .expect("Failed parsing a ntax value");
        ntax.parse::<usize>()
            .expect("Failed parsing a ntax value to integer")
    }

    fn check_concat(&self) -> Option<PartitionFmt> {
        if self.matches.is_present("concat") {
            Some(self.get_part_fmt())
        } else {
            None
        }
    }

    fn get_part_fmt(&self) -> PartitionFmt {
        if self.matches.is_present("partition") {
            self.parse_partition_fmt(self.matches)
        } else {
            PartitionFmt::Nexus
        }
    }

    fn set_output_path(&mut self) {
        if self.matches.is_present("output") {
            self.output_dir = self.parse_output(self.matches);
        } else {
            match self.input_dir.as_ref() {
                Some(dir) => self.output_dir = self.fmt_output_path(dir),
                None => panic!("Please, define an output directory!"),
            }
        }
    }

    fn set_multi_output_path(&mut self) {
        if self.matches.is_present("output") {
            let output_dir = self.parse_output(self.matches);
            self.output_dir = self.fmt_output_path(&output_dir)
        } else {
            match self.input_dir.as_ref() {
                Some(dir) => self.output_dir = self.fmt_output_path(dir),
                None => panic!("Please, define an output directory!"),
            }
        }
    }

    fn fmt_output_path(&self, dir: &Path) -> PathBuf {
        let parent = dir.parent().expect("Failed parsing input directory");
        let last: String = match dir.file_name() {
            Some(fname) => fname.to_string_lossy().to_string(),
            None => String::from("segul-filter"),
        };
        let output_dir = match self.params {
            Params::MinTax(_) => format!("{}_{}p", last, self.percent * 100.0),
            Params::AlnLen(len) => format!("{}_{}bp", last, len),
            Params::ParsInf(inf) => format!("{}_{}inf", last, inf),
            Params::PercInf(perc_inf) => format!("{}_{}percent_inf", last, perc_inf * 100.0),
            Params::TaxonAll(_) => format!("{}_taxon_id", last),
        };
        parent.join(output_dir)
    }

    fn print_params(&self) {
        log::info!("{}", Yellow.paint("Params"));
        match &self.params {
            Params::MinTax(min_taxa) => {
                log::info!("{:18}: {}", "Taxon count", self.ntax);
                log::info!("{:18}: {}%", "Percent", self.percent * 100.0);
                log::info!("{:18}: {}\n", "Min tax", min_taxa);
            }
            Params::AlnLen(len) => log::info!("{:18}: {} bp\n", "Min aln len", len),
            Params::ParsInf(inf) => log::info!("{:18}: {}\n", "Min pars. inf", inf),
            Params::PercInf(perc_inf) => {
                log::info!("{:18}: {}%\n", "% pars. inf", perc_inf * 100.0)
            }
            Params::TaxonAll(taxon_id) => {
                log::info!("{:18}: {} taxa\n", "Taxon id", taxon_id.len())
            }
        }
    }

    fn is_npercent(&mut self) -> bool {
        self.matches.is_present("npercent")
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use clap::{Arg, Command};

    #[test]
    fn test_min_taxa_output_dir() {
        let arg = Command::new("segul-test")
            .arg(Arg::new("test"))
            .get_matches();
        let mut min_taxa = FilterParser::new(&arg);
        let dir = "./test_taxa/";
        min_taxa.percent = 0.75;
        let res = PathBuf::from("./test_taxa_75p");
        let output = min_taxa.fmt_output_path(Path::new(dir));
        assert_eq!(res, output);
    }

    #[test]
    fn test_min_taxa() {
        let arg = Command::new("segul-test")
            .arg(Arg::new("filter-test"))
            .get_matches();
        let mut filter = FilterParser::new(&arg);
        filter.percent = 0.65;
        filter.ntax = 10;
        assert_eq!(6, filter.count_min_tax());
    }
}

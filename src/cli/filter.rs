use std::path::PathBuf;

use ansi_term::Colour::Yellow;
use clap::ArgMatches;

use crate::cli::*;
use crate::core::filter;
use crate::helper::types::{DataType, InputFmt, OutputFmt};

impl InputCli for FilterParser<'_> {}
impl InputPrint for FilterParser<'_> {}
impl OutputCli for FilterParser<'_> {}
impl PartCLi for FilterParser<'_> {}

pub(in crate::cli) struct FilterParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_fmt: InputFmt,
    input_dir: Option<PathBuf>,
    output_dir: PathBuf,
    files: Vec<PathBuf>,
    params: filter::Params,
    ntax: usize,
    percent: f64,
    datatype: DataType,
}

impl<'a> FilterParser<'a> {
    pub(in crate::cli) fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_fmt: InputFmt::Fasta,
            input_dir: None,
            output_dir: PathBuf::new(),
            files: Vec::new(),
            params: filter::Params::MinTax(0),
            ntax: 0,
            percent: 0.0,
            datatype: DataType::Dna,
        }
    }

    pub(in crate::cli) fn filter(&mut self) {
        self.input_fmt = self.parse_input_fmt(self.matches);
        self.datatype = self.parse_datatype(self.matches);
        let task_desc = "Alignment filtering";
        self.files = if self.is_input_dir() {
            let dir = self.parse_dir_input(self.matches);
            self.input_dir = Some(PathBuf::from(dir));
            self.get_files(dir, &self.input_fmt)
        } else {
            self.parse_input_wcard(self.matches)
        };

        if self.is_npercent() {
            self.get_min_taxa_npercent(task_desc);
        } else {
            self.get_params();
            self.set_output_path();
            self.filter_aln(task_desc);
        }
    }

    fn get_min_taxa_npercent(&mut self, task_desc: &str) {
        let npercent = self.parse_npercent();
        npercent.iter().for_each(|&np| {
            self.percent = np;
            let min_tax = self.get_min_taxa();
            self.params = filter::Params::MinTax(min_tax);
            self.set_multi_output_path();
            self.filter_aln(task_desc);
            utils::print_divider();
        });
    }

    fn is_input_dir(&self) -> bool {
        self.matches.is_present("dir")
    }

    fn filter_aln(&self, task_desc: &str) {
        self.print_input_multi(
            &self.input_dir,
            task_desc,
            self.files.len(),
            &self.input_fmt,
        );
        self.print_params();
        let mut filter = filter::SeqFilter::new(
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
                filter.set_concat(&output_fmt, &part_fmt);
                filter.filter_aln();
            }
            None => filter.filter_aln(),
        }
    }

    fn get_params(&mut self) {
        self.params = match self.matches {
            m if m.is_present("percent") => {
                self.percent = self.get_percent();
                filter::Params::MinTax(self.get_min_taxa())
            }
            m if m.is_present("aln-len") => filter::Params::AlnLen(self.get_aln_len()),
            m if m.is_present("pars-inf") => filter::Params::ParsInf(self.get_pars_inf()),
            _ => unreachable!("Invalid parameters!"),
        }
    }

    fn get_min_taxa(&mut self) -> usize {
        self.get_ntax();
        self.count_min_tax()
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

    fn get_aln_len(&self) -> usize {
        let len = self
            .matches
            .value_of("aln-len")
            .expect("CANNOT GET ALIGNMENT LENGTH VALUES");
        len.parse::<usize>()
            .expect("CANNOT PARSE ALIGNMENT LENGTH VALUES TO INTEGERS")
    }

    fn get_pars_inf(&self) -> usize {
        let len = self
            .matches
            .value_of("pars-inf")
            .expect("CANNOT GET PARSIMONY INFORMATIVE VALUES");
        len.parse::<usize>()
            .expect("CANNOT PARSE PARSIMONY INFORMATIVE VALUES TO INTEGERS")
    }

    fn get_ntax(&mut self) {
        self.ntax = if self.matches.is_present("ntax") {
            self.parse_ntax()
        } else {
            IDs::new(&self.files, &self.input_fmt, &self.datatype)
                .get_id_all()
                .len()
        };
    }

    fn count_min_tax(&self) -> usize {
        (self.ntax as f64 * self.percent).floor() as usize
    }

    fn parse_npercent(&self) -> Vec<f64> {
        self.matches
            .values_of("npercent")
            .expect("FAILED PARSING npercent")
            .map(|np| self.parse_percent(np))
            .collect()
    }

    fn is_npercent(&mut self) -> bool {
        self.matches.is_present("npercent")
    }

    fn get_percent(&mut self) -> f64 {
        let percent = self
            .matches
            .value_of("percent")
            .expect("CANNOT GET PERCENTAGE VALUES");
        self.parse_percent(percent)
    }

    fn parse_percent(&self, percent: &str) -> f64 {
        percent
            .parse::<f64>()
            .expect("CANNOT PARSE PERCENTAGE VALUES TO FLOATING POINTS")
    }

    fn parse_ntax(&self) -> usize {
        let ntax = self
            .matches
            .value_of("ntax")
            .expect("CANNOT GET NTAX VALUES");
        ntax.parse::<usize>()
            .expect("CANNOT PARSE NTAX VALUES TO INTEGERS")
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
            filter::Params::MinTax(_) => format!("{}_{}p", last, self.percent * 100.0),
            filter::Params::AlnLen(len) => format!("{}_{}bp", last, len),
            filter::Params::ParsInf(inf) => format!("{}_{}inf", last, inf),
        };
        parent.join(output_dir)
    }

    fn print_params(&self) {
        log::info!("{}", Yellow.paint("Parameters"));
        match self.params {
            filter::Params::MinTax(min_taxa) => {
                log::info!("{:18}: {}", "Taxon count", self.ntax);
                log::info!("{:18}: {}%", "Percent", self.percent * 100.0);
                log::info!("{:18}: {}\n", "Min tax", min_taxa);
            }
            filter::Params::AlnLen(len) => log::info!("{:18}: {}bp\n", "Min aln len", len),
            filter::Params::ParsInf(inf) => log::info!("{:18}: {}\n", "Min pars. inf", inf),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use clap::{App, Arg};

    #[test]
    fn min_taxa_output_dir_test() {
        let arg = App::new("segul-test")
            .arg(Arg::with_name("test"))
            .get_matches();
        let mut min_taxa = FilterParser::new(&arg);
        let dir = "./test_taxa/";
        min_taxa.percent = 0.75;
        let res = PathBuf::from("./test_taxa_75p");
        let output = min_taxa.fmt_output_path(Path::new(dir));
        assert_eq!(res, output);
    }

    #[test]
    fn min_taxa_test() {
        let arg = App::new("segul-test")
            .arg(Arg::with_name("filter-test"))
            .get_matches();
        let mut filter = FilterParser::new(&arg);
        filter.percent = 0.65;
        filter.ntax = 10;
        assert_eq!(6, filter.count_min_tax());
    }
}

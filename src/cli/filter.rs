use std::path::PathBuf;

use colored::Colorize;

use crate::cli::{AlignSeqInput, ConcatCli, InputCli, OutputCli};
use crate::handler::align::filter::{Params, SeqFilter};
use crate::helper::finder::IDs;
use crate::helper::types::{DataType, InputFmt, PartitionFmt};
use crate::helper::{files, utils};
use crate::parser::txt;

use super::args::AlignFilterArgs;
use super::{collect_paths, AlignSeqPrint};

impl InputCli for FilterParser<'_> {}
impl OutputCli for FilterParser<'_> {}
impl ConcatCli for FilterParser<'_> {}
impl AlignSeqInput for FilterParser<'_> {}

pub(in crate::cli) struct FilterParser<'a> {
    args: &'a AlignFilterArgs,
    input_dir: Option<PathBuf>,
    input_fmt: InputFmt,
    output_dir: PathBuf,
    files: Vec<PathBuf>,
    params: Params,
    ntax: usize,
    percent: f64,
    datatype: DataType,
}

impl<'a> FilterParser<'a> {
    pub(in crate::cli) fn new(args: &'a AlignFilterArgs) -> Self {
        Self {
            args,
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
        self.input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        self.datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let task_desc = "Alignment filtering";
        let dir = &self.args.io.dir;
        let input_fmt = &self.input_fmt; // Binding to satisfy the macro
        self.files = collect_paths!(self, dir, input_fmt);
        AlignSeqPrint::new(
            &self.input_dir,
            input_fmt,
            &self.datatype,
            task_desc,
            self.files.len(),
        )
        .print();

        if let Some(npercent) = &self.args.npercent {
            self.filter_min_taxa_npercent(npercent);
        } else {
            self.parse_params();
            self.fmt_output_path();
            self.filter_aln();
        }
    }

    fn filter_min_taxa_npercent(&mut self, npercent: &[f64]) {
        self.count_ntax();
        npercent.iter().for_each(|np| {
            self.percent = *np;
            let min_taxa = self.count_min_tax();
            self.params = Params::MinTax(min_taxa);
            self.fmt_output_path();
            self.filter_aln();
            utils::print_divider();
        });
    }

    fn filter_aln(&self) {
        self.check_output_dir_exist(&self.output_dir, self.args.io.force);
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
                let output_fmt = self.parse_output_fmt(&self.args.out_fmt.output_fmt);
                let prefix = self.parse_prefix(&self.args.partition.prefix, &self.output_dir);
                let output = files::create_output_fname(&self.output_dir, &prefix, &output_fmt);
                filter.set_concat(&output, &output_fmt, &part_fmt);
                filter.filter_aln();
            }
            None => filter.filter_aln(),
        }
    }

    fn parse_params(&mut self) {
        self.params = match self.args {
            m if m.percent.is_some() => self.parse_percent(),
            m if m.len.is_some() => Params::AlnLen(self.parse_aln_len()),
            m if m.pinf.is_some() => Params::ParsInf(self.parse_pars_inf()),
            m if m.percent_inf.is_some() => Params::PercInf(self.count_percent_inf()),
            m if m.ids.is_some() => Params::TaxonAll(self.parse_taxon_id()),
            _ => unreachable!("Invalid parameters!"),
        }
    }

    fn parse_percent(&mut self) -> Params {
        match self.args.percent {
            Some(percent) => {
                self.percent = percent;
                self.count_ntax();
                let min_taxa = self.count_min_tax();
                Params::MinTax(min_taxa)
            }
            None => unreachable!("Invalid parameters!"),
        }
    }

    fn parse_taxon_id(&self) -> Vec<String> {
        match &self.args.ids {
            Some(path) => txt::parse_text_file(path),
            None => unreachable!("Invalid parameters!"),
        }
    }

    fn count_percent_inf(&self) -> f64 {
        match self.args.percent_inf {
            Some(p) => p,
            None => unreachable!("Invalid parameters!"),
        }
    }

    fn parse_aln_len(&self) -> usize {
        match self.args.len {
            Some(len) => len,
            None => unreachable!("Invalid parameters!"),
        }
    }

    fn parse_pars_inf(&self) -> usize {
        match self.args.pinf {
            Some(pinf) => pinf,
            None => unreachable!("Invalid parameters!"),
        }
    }

    fn count_ntax(&mut self) {
        if let Some(ntax) = self.args.ntax {
            self.ntax = ntax;
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

    fn check_concat(&self) -> Option<PartitionFmt> {
        if self.args.concat {
            Some(self.parse_partition_fmt(&self.args.partition.part_fmt, self.args.partition.codon))
        } else {
            None
        }
    }

    fn fmt_output_path(&mut self) {
        let output_dir = match self.params {
            Params::MinTax(_) => {
                format!("{}_{}p", self.args.output, self.percent * 100.0)
            }
            Params::AlnLen(len) => format!("{}_{}bp", self.args.output, len),
            Params::ParsInf(inf) => format!("{}_{}inf", self.args.output, inf),
            Params::PercInf(perc_inf) => {
                format!("{}_{}percent_inf", self.args.output, perc_inf * 100.0)
            }
            Params::TaxonAll(_) => format!("{}_taxonID", self.args.output),
        };

        self.output_dir = PathBuf::from(output_dir);
    }

    fn print_params(&self) {
        log::info!("{}", "Params".yellow());
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
}

#[cfg(test)]
mod test {
    use crate::cli::args::{CommonConcatArgs, CommonSeqInput, CommonSeqOutput, IOArgs};

    use super::*;

    macro_rules! args {
        ($args: ident) => {
            let $args = AlignFilterArgs {
                io: IOArgs {
                    input: None,
                    dir: Some("./test_taxa/".to_string()),
                    force: false,
                },
                percent: Some(0.75),
                ntax: None,
                len: None,
                pinf: None,
                percent_inf: None,
                ids: None,
                concat: false,
                partition: CommonConcatArgs {
                    part_fmt: "raxml".to_string(),
                    codon: false,
                    prefix: None,
                },
                in_fmt: CommonSeqInput {
                    input_fmt: "phylip".to_string(),

                    datatype: "dna".to_string(),
                },
                out_fmt: CommonSeqOutput {
                    output_fmt: "phylip".to_string(),
                },
                output: "SEGUL-Filter".to_string(),
                npercent: None,
            };
        };
    }

    #[test]
    fn test_min_taxa_output_dir() {
        args!(args);
        let mut min_taxa = FilterParser::new(&args);
        let res = PathBuf::from("SEGUL-Filter_75p");
        min_taxa.parse_params();
        min_taxa.fmt_output_path();
        assert_eq!(res, min_taxa.output_dir);
    }

    #[test]
    fn test_min_taxa() {
        args!(args);
        let mut filter = FilterParser::new(&args);
        filter.percent = 0.65;
        filter.ntax = 10;
        assert_eq!(6, filter.count_min_tax());
    }
}

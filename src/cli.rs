use std::fs::File;
use std::io::{self, BufWriter, Result, Write};
use std::path::{Path, PathBuf};

use clap::{App, AppSettings, Arg, ArgMatches};
use indexmap::IndexSet;
use rayon::prelude::*;

use crate::core::converter::Converter;
use crate::core::filter;
use crate::core::msa;
use crate::core::summary::SeqStats;
use crate::helper::common::{PartitionFormat, SeqFormat};
use crate::helper::finder::{Files, IDs};
use crate::helper::utils;

fn get_args(version: &str) -> ArgMatches {
    App::new("segul")
        .version(version)
        .about("An ultra-fast and efficient alignment manipulation tool")
        .author("Heru Handika")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .subcommand(
            App::new("convert")
                .about("Converts sequence formats")
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .long("input")
                        .help("Convert a fasta file")
                        .takes_value(true)
                        .required_unless("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT FILE"),
                )
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Convert multiple fasta files inside a dir")
                        .takes_value(true)
                        .required_unless("input")
                        .conflicts_with("input")
                        .value_name("DIR"),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Sets target directory or use a costume file name for a single input")
                        .takes_value(true)
                        .required_unless("input")
                        .value_name("OUTPUT"),
                )
                .arg(
                    Arg::with_name("format")
                        .short("f")
                        .long("format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("FORMAT")
                        .possible_values(&[
                            "fasta",
                            "nexus",
                            "phylip",
                            "fasta-int",
                            "nexus-int",
                            "phylip-int",
                        ]),
                )
                .arg(
                    Arg::with_name("output-format")
                        .short("F")
                        .long("output-format")
                        .help("Sets target output format")
                        .takes_value(true)
                        .default_value("nexus")
                        .value_name("FORMAT")
                        .possible_values(&[
                            "nexus",
                            "phylip",
                            "fasta",
                            "fasta-int",
                            "nexus-int",
                            "phylip-int",
                        ]),
                )
                .arg(
                    Arg::with_name("sort")
                        .long("sort")
                        .help("Sorts the alignments")
                        .takes_value(false),
                ),
        )
        .subcommand(
            App::new("concat")
                .about("Concatenates alignments")
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Inputs alignment dir")
                        .takes_value(true)
                        .required(true)
                        .value_name("DIR"),
                )
                .arg(
                    Arg::with_name("format")
                        .short("f")
                        .long("format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("FORMAT")
                        .possible_values(&[
                            "nexus",
                            "phylip",
                            "fasta",
                            "fasta-int",
                            "nexus-int",
                            "phylip-int",
                        ]),
                )
                .arg(
                    Arg::with_name("partition")
                        .short("-p")
                        .long("part")
                        .help("Sets partition format")
                        .takes_value(true)
                        .required(true)
                        .default_value("charset")
                        .value_name("FORMAT")
                        .possible_values(&["charset", "nexus", "raxml"]),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Uses a costume output filename")
                        .takes_value(true)
                        .required(true)
                        .default_value("concat")
                        .value_name("OUTPUT"),
                )
                .arg(
                    Arg::with_name("output-format")
                        .short("F")
                        .long("output-format")
                        .help("Sets target output format")
                        .takes_value(true)
                        .default_value("nexus")
                        .value_name("FORMAT")
                        .possible_values(&[
                            "nexus",
                            "phylip",
                            "fasta",
                            "fasta-int",
                            "nexus-int",
                            "phylip-int",
                        ]),
                )
                .arg(
                    Arg::with_name("codon")
                        .long("codon")
                        .help("Sets codon model partition format")
                        .takes_value(false),
                ),
        )
        .subcommand(
            App::new("filter")
                .about("Picks alignments with specified min taxa")
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Inputs a path to alignment dir")
                        .takes_value(true)
                        .required(true)
                        .value_name("DIR"),
                )
                .arg(
                    Arg::with_name("format")
                        .short("f")
                        .long("format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("FORMAT")
                        .possible_values(&[
                            "fasta",
                            "nexus",
                            "phylip",
                            "fasta-int",
                            "nexus-int",
                            "phylip-int",
                        ]),
                )
                .arg(
                    Arg::with_name("percent")
                        .long("percent")
                        .help("Sets percentage of minimal taxa")
                        .takes_value(true)
                        .required_unless_all(&["npercent", "aln-len", "pars-inf"])
                        .conflicts_with_all(&["npercent", "aln-len", "pars-inf"])
                        .value_name("FORMAT"),
                )
                .arg(
                    Arg::with_name("npercent")
                        .long("npercent")
                        .help("Sets minimal taxa in multiple percentages")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "aln-len", "pars-inf"])
                        .multiple(true)
                        .value_name("FORMAT"),
                )
                .arg(
                    Arg::with_name("ntax")
                        .long("ntax")
                        .help("Inputs the total number of taxa")
                        .takes_value(true)
                        .conflicts_with_all(&["aln-len", "pars-inf"])
                        .value_name("TAXON-COUNT"),
                )
                .arg(
                    Arg::with_name("aln-len")
                        .long("len")
                        .help("Sets minimal alignment length")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "npercent", "pars-inf"])
                        .value_name("FORMAT"),
                )
                .arg(
                    Arg::with_name("pars-inf")
                        .long("pinf")
                        .help("Sets minimal alignment length")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "npercent", "aln-len"])
                        .value_name("FORMAT"),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Sets an output directory")
                        .takes_value(true)
                        .value_name("OUTPUT"),
                )
                .arg(
                    Arg::with_name("output-format")
                        .short("F")
                        .long("output-format")
                        .help("Sets output format if concat")
                        .takes_value(true)
                        .value_name("FORMAT")
                        .possible_values(&[
                            "fasta",
                            "nexus",
                            "phylip",
                            "fasta-int",
                            "nexus-int",
                            "phylip-int",
                        ]),
                )
                .arg(
                    Arg::with_name("partition")
                        .short("-p")
                        .long("part")
                        .help("Sets partition format if concat")
                        .takes_value(true)
                        .value_name("FORMAT")
                        .possible_values(&["charset", "nexus", "raxml"]),
                )
                .arg(
                    Arg::with_name("concat")
                        .long("concat")
                        .help("Concats the final results")
                        .required_ifs(&[
                            ("filter", "partition"),
                            ("filter", "codon"),
                            ("filter", "output-format"),
                        ])
                        .takes_value(false),
                )
                .arg(
                    Arg::with_name("codon")
                        .long("codon")
                        .help("Sets codon model partition format")
                        .takes_value(false),
                ),
        )
        .subcommand(
            App::new("id")
                .about("Gets sample ids from multiple alignments")
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Inputs dir with alignment files")
                        .takes_value(true)
                        .value_name("INPUT FILE"),
                )
                .arg(
                    Arg::with_name("format")
                        .short("f")
                        .long("format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("FORMAT")
                        .possible_values(&[
                            "fasta",
                            "nexus",
                            "phylip",
                            "fasta-int",
                            "nexus-int",
                            "phylip-int",
                        ]),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Uses a costume output filename")
                        .takes_value(true)
                        .required(true)
                        .default_value("id")
                        .value_name("OUTPUT"),
                ),
        )
        .subcommand(
            App::new("summary")
                .about("Gets alignment summary stats")
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .long("input")
                        .help("Convert a fasta file")
                        .takes_value(true)
                        .required_unless("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT FILE"),
                )
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Inputs dir with alignment files")
                        .takes_value(true)
                        .conflicts_with("input")
                        .value_name("INPUT FILE"),
                )
                .arg(
                    Arg::with_name("format")
                        .short("f")
                        .long("format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("FORMAT")
                        .possible_values(&[
                            "fasta",
                            "nexus",
                            "phylip",
                            "fasta-int",
                            "nexus-int",
                            "phylip-int",
                        ]),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Uses a costume output filename")
                        .takes_value(true)
                        .required(true)
                        .default_value("SEGUL-stats")
                        .value_name("OUTPUT"),
                )
                .arg(
                    Arg::with_name("decrement")
                        .long("decrement")
                        .help("Sets a custom percentage decrement value for counting taxon completeness")
                        .takes_value(true)
                        .value_name("DECREMENT")
                        .default_value("5")
                        .possible_values(&["1", "2", "5", "10"]),
                ),
        )
        .get_matches()
}

pub fn parse_cli(version: &str) {
    let args = get_args(version);
    let text = format!("SEGUL v{}", version);
    utils::print_title(&text);
    match args.subcommand() {
        ("convert", Some(convert_matches)) => ConvertParser::new(convert_matches).convert(),
        ("concat", Some(concat_matches)) => ConcatParser::new(concat_matches).concat(),
        ("filter", Some(pick_matches)) => FilterParser::new(pick_matches).filter(),
        ("id", Some(id_matches)) => IdParser::new(id_matches).get_id(),
        ("summary", Some(stats_matches)) => StatsParser::new(stats_matches).show_stats(),
        _ => unreachable!(),
    }
}

trait Cli {
    fn get_file_input<'a>(&self, matches: &'a ArgMatches) -> &'a str {
        matches
            .value_of("input")
            .expect("CANNOT FIND AN INPUT FILE")
    }

    fn get_dir_input<'a>(&self, matches: &'a ArgMatches) -> &'a str {
        matches.value_of("dir").expect("CANNOT READ DIR PATH")
    }

    fn get_files(&self, dir: &str, input_fmt: &SeqFormat) -> Vec<PathBuf> {
        Files::new(dir, input_fmt).get_files()
    }

    fn get_output<'a>(&self, matches: &'a ArgMatches) -> &'a str {
        matches.value_of("output").expect("CANNOT READ OUTPUT PATH")
    }

    fn get_output_path(&self, matches: &ArgMatches) -> PathBuf {
        if matches.is_present("output") {
            let output = self.get_output(matches);
            PathBuf::from(output)
        } else {
            PathBuf::from(".")
        }
    }

    fn get_input_fmt(&self, matches: &ArgMatches) -> SeqFormat {
        let input_fmt = matches
            .value_of("format")
            .expect("CANNOT READ FORMAT INPUT");
        match input_fmt {
            "fasta" | "fasta-int" => SeqFormat::Fasta,
            "nexus" | "nexus-int" => SeqFormat::Nexus,
            "phylip" => SeqFormat::Phylip,
            "phylip-int" => SeqFormat::PhylipInt,
            _ => panic!(
                "UNSUPPORTED FORMAT. \
        THE PROGRAM ONLY ACCEPT fasta, fasta-int, nexus, nexus-int, phylip, and phylip-int. ALL IN lowercase. \
        YOUR INPUT: {} ",
                input_fmt
            ),
        }
    }

    fn get_output_fmt(&self, matches: &ArgMatches) -> SeqFormat {
        let output_fmt = matches
            .value_of("output-format")
            .expect("CANNOT READ FORMAT INPUT");
        match output_fmt {
            "nexus" => SeqFormat::Nexus,
            "phylip" => SeqFormat::Phylip,
            "fasta" => SeqFormat::Fasta,
            "nexus-int" => SeqFormat::NexusInt,
            "fasta-int" => SeqFormat::FastaInt,
            "phylip-int" => SeqFormat::PhylipInt,
            _ => panic!(
                "UNSUPPORTED FORMAT. \
        THE PROGRAM ONLY ACCEPT fasta, fasta-int, nexus, nexus-int, phylip, and phylip-int. ALL IN lowercase. \
        YOUR INPUT: {} ",
                output_fmt
            ),
        }
    }
}

impl Cli for ConvertParser<'_> {}

struct ConvertParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_fmt: SeqFormat,
    is_dir: bool,
}

impl<'a> ConvertParser<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_fmt: SeqFormat::Fasta,
            is_dir: false,
        }
    }

    fn convert(&mut self) {
        self.input_fmt = self.get_input_fmt(&self.matches);
        if self.matches.is_present("input") {
            self.convert_file();
        } else {
            self.convert_multiple_fasta();
        }
    }

    fn convert_file(&mut self) {
        let input = Path::new(self.get_file_input(self.matches));
        let output_fmt = self.get_output_fmt(self.matches);
        let output = self.get_output_path(self.matches);
        self.print_input_file(input).unwrap();

        self.convert_any(input, &output, &output_fmt);
    }

    fn convert_multiple_fasta(&mut self) {
        let dir = self.get_dir_input(self.matches);
        let files = self.get_files(dir, &self.input_fmt);
        let output_fmt = self.get_output_fmt(self.matches);
        let output = self.get_output_path(&self.matches);
        self.is_dir = true;
        self.print_input_dir(Path::new(dir), files.len(), &output)
            .unwrap();
        let spin = utils::set_spinner();
        spin.set_message("Converting alignments...");
        files.par_iter().for_each(|file| {
            let output = output.join(file.file_stem().unwrap());
            self.convert_any(file, &output, &output_fmt);
        });
        spin.finish_with_message("DONE!");
    }

    fn convert_any(&self, input: &Path, output: &Path, output_fmt: &SeqFormat) {
        let mut convert = Converter::new(input, output, output_fmt, self.is_dir);
        if self.is_sort() {
            convert.convert_sorted(&self.input_fmt);
        } else {
            convert.convert_unsorted(&self.input_fmt);
        }
    }

    fn is_sort(&self) -> bool {
        self.matches.is_present("sort")
    }

    fn print_input_file(&self, input: &Path) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul convert")?;
        writeln!(writer, "Input\t\t: {}", &input.display())?;
        Ok(())
    }

    fn print_input_dir(&self, input: &Path, nfile: usize, output: &Path) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul convert")?;
        writeln!(writer, "Input dir\t: {}", &input.display())?;
        writeln!(writer, "Total files\t: {}", utils::fmt_num(&nfile))?;
        writeln!(writer, "Output dir \t: {}\n", output.display())?;
        Ok(())
    }
}

trait PartCLi {
    fn parse_partition_fmt(&self, matches: &ArgMatches) -> PartitionFormat {
        let part_fmt = matches
            .value_of("partition")
            .expect("CANNOT READ PARTITION FORMAT");
        if matches.is_present("codon") {
            self.parse_partition_fmt_codon(part_fmt)
        } else {
            self.parse_partition_fmt_std(part_fmt)
        }
    }

    fn parse_partition_fmt_std(&self, part_fmt: &str) -> PartitionFormat {
        match part_fmt {
            "nexus" => PartitionFormat::Nexus,
            "raxml" => PartitionFormat::Raxml,
            "charset" => PartitionFormat::Charset,
            _ => PartitionFormat::Nexus,
        }
    }

    fn parse_partition_fmt_codon(&self, part_fmt: &str) -> PartitionFormat {
        match part_fmt {
            "charset" => PartitionFormat::CharsetCodon,
            "nexus" => PartitionFormat::NexusCodon,
            "raxml" => PartitionFormat::RaxmlCodon,
            _ => PartitionFormat::NexusCodon,
        }
    }

    fn check_partition_format(&self, output_fmt: &SeqFormat, part_fmt: &PartitionFormat) {
        match output_fmt {
            SeqFormat::Nexus | SeqFormat::NexusInt => (),
            _ => {
                if let PartitionFormat::Nexus | PartitionFormat::NexusCodon = part_fmt {
                    panic!(
                        "CANNOT WRITE EMBEDDED-NEXUS PARTITION TO NON-NEXUS OUTPUT. \
                MAYBE YOU MEAN TO WRITE THE PARTITION TO 'charset' INSTEAD."
                    )
                }
            }
        }
    }
}

impl PartCLi for ConcatParser<'_> {}

impl Cli for ConcatParser<'_> {
    fn get_output_path(&self, matches: &ArgMatches) -> PathBuf {
        PathBuf::from(self.get_output(matches))
    }
}

struct ConcatParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_fmt: SeqFormat,
    output_fmt: SeqFormat,
    part_fmt: PartitionFormat,
}

impl<'a> ConcatParser<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_fmt: SeqFormat::Fasta,
            output_fmt: SeqFormat::Nexus,
            part_fmt: PartitionFormat::Charset,
        }
    }

    fn concat(&mut self) {
        self.input_fmt = self.get_input_fmt(self.matches);
        let dir = self.get_dir_input(self.matches);
        let output = self.get_output(self.matches);
        self.output_fmt = self.get_output_fmt(self.matches);
        self.part_fmt = self.parse_partition_fmt(self.matches);
        self.check_partition_format(&self.output_fmt, &self.part_fmt);
        self.print_input_dir(&dir).unwrap();
        let concat =
            msa::MSAlignment::new(&self.input_fmt, output, &self.output_fmt, &self.part_fmt);
        let mut files = self.get_files(dir, &self.input_fmt);
        concat.concat_alignment(&mut files);
    }

    fn print_input_dir(&self, input: &str) -> Result<()> {
        let io = io::stdout();
        let mut writer = io::BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul concat")?;
        writeln!(writer, "Input dir\t: {}\n", input)?;
        Ok(())
    }
}

impl Cli for FilterParser<'_> {}
impl PartCLi for FilterParser<'_> {}

struct FilterParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_fmt: SeqFormat,
    output_dir: PathBuf,
    files: Vec<PathBuf>,
    params: filter::Params,
    ntax: usize,
    percent: f64,
}

impl<'a> FilterParser<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_fmt: SeqFormat::Fasta,
            output_dir: PathBuf::new(),
            files: Vec::new(),
            params: filter::Params::MinTax(0),
            ntax: 0,
            percent: 0.0,
        }
    }

    fn filter(&mut self) {
        self.input_fmt = self.get_input_fmt(self.matches);
        let dir = self.get_dir_input(self.matches);
        self.files = self.get_files(dir, &self.input_fmt);
        if self.is_npercent() {
            self.get_min_taxa_npercent(dir);
        } else {
            self.get_params();
            self.set_output_path(dir);
            self.print_input(dir).expect("CANNOT DISPLAY TO STDOUT");
            self.filter_aln();
        }
    }

    fn get_min_taxa_npercent(&mut self, dir: &str) {
        let npercent = self.parse_npercent();
        npercent.iter().for_each(|&np| {
            self.percent = np;
            let min_tax = self.get_min_taxa();
            self.params = filter::Params::MinTax(min_tax);
            self.set_multi_output_path(dir);
            self.print_input(dir).expect("CANNOT DISPLAY TO STDOUT");
            self.filter_aln();
            utils::print_divider();
        });
    }

    fn filter_aln(&self) {
        let mut filter =
            filter::SeqFilter::new(&self.files, &self.input_fmt, &self.output_dir, &self.params);
        match self.check_concat() {
            Some(part_fmt) => {
                let output_fmt = if self.matches.is_present("output-format") {
                    self.get_output_fmt(self.matches)
                } else {
                    SeqFormat::Nexus
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
            _ => unreachable!(),
        }
    }

    fn get_min_taxa(&mut self) -> usize {
        self.get_ntax();
        self.count_min_tax()
    }

    fn check_concat(&self) -> Option<PartitionFormat> {
        if self.matches.is_present("concat") {
            Some(self.get_part_fmt())
        } else {
            None
        }
    }

    fn get_part_fmt(&self) -> PartitionFormat {
        if self.matches.is_present("partition") {
            self.parse_partition_fmt(self.matches)
        } else {
            PartitionFormat::Charset
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
            IDs::new(&self.files, &self.input_fmt).get_id_all().len()
        };
    }

    fn count_min_tax(&self) -> usize {
        (self.ntax as f64 * self.percent).floor() as usize
    }

    fn parse_npercent(&self) -> Vec<f64> {
        let npercent: Vec<&str> = self.matches.values_of("npercent").unwrap().collect();
        npercent.iter().map(|np| self.parse_percent(np)).collect()
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

    fn set_output_path<P: AsRef<Path>>(&mut self, dir: P) {
        if self.matches.is_present("output") {
            self.output_dir = PathBuf::from(self.get_output(self.matches));
        } else {
            self.output_dir = self.fmt_output_path(dir.as_ref());
        }
    }

    fn set_multi_output_path<P: AsRef<Path>>(&mut self, dir: P) {
        if self.matches.is_present("output") {
            let output_dir = PathBuf::from(self.get_output(self.matches));
            self.output_dir = self.fmt_output_path(&output_dir)
        } else {
            self.output_dir = self.fmt_output_path(dir.as_ref());
        }
    }

    fn fmt_output_path(&self, dir: &Path) -> PathBuf {
        let parent = dir.parent().unwrap();
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

    fn print_input(&self, dir: &str) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "\x1b[0;33mInput\x1b[0m")?;
        writeln!(writer, "Dir\t\t: {}", dir)?;
        writeln!(
            writer,
            "File count\t: {}",
            utils::fmt_num(&self.files.len())
        )?;
        match self.params {
            filter::Params::MinTax(min_taxa) => {
                writeln!(writer, "Taxon count\t: {}", self.ntax)?;
                writeln!(writer, "Percent\t\t: {}%", self.percent * 100.0)?;
                writeln!(writer, "Min tax\t\t: {}", min_taxa)?;
            }
            filter::Params::AlnLen(len) => writeln!(writer, "Min aln len\t: {}bp", len)?,
            filter::Params::ParsInf(inf) => writeln!(writer, "Min pars. inf\t: {}", inf)?,
        }
        writeln!(writer)?;
        Ok(())
    }
}

impl Cli for IdParser<'_> {
    fn get_output_path(&self, matches: &ArgMatches) -> PathBuf {
        if matches.is_present("output") {
            let output = self.get_output(matches);
            PathBuf::from(output).with_extension("txt")
        } else {
            let input = Path::new(self.get_dir_input(matches));
            input.with_extension("txt")
        }
    }
}

struct IdParser<'a> {
    matches: &'a ArgMatches<'a>,
}

impl<'a> IdParser<'a> {
    fn new(matches: &'a ArgMatches) -> Self {
        Self { matches }
    }

    fn get_id(&self) {
        let dir = self.get_dir_input(&self.matches);
        let input_fmt = self.get_input_fmt(&self.matches);
        let files = self.get_files(dir, &input_fmt);
        self.print_input(dir).unwrap();
        let spin = utils::set_spinner();
        spin.set_message("Indexing IDs..");
        let ids = IDs::new(&files, &input_fmt).get_id_all();
        spin.finish_with_message("DONE!");
        self.write_results(&ids);
    }

    fn write_results(&self, ids: &IndexSet<String>) {
        let fname = self.get_output_path(&self.matches);
        let file = File::create(&fname).expect("CANNOT CREATE AN OUTPUT FILE");
        let mut writer = BufWriter::new(file);
        ids.iter().for_each(|id| {
            writeln!(writer, "{}", id).unwrap();
        });
        writer.flush().unwrap();
        self.print_output(&fname, ids.len()).unwrap();
    }

    fn print_input(&self, dir: &str) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Command\t\t\t: segul id")?;
        writeln!(writer, "Input dir\t\t: {}\n", dir)?;

        Ok(())
    }

    fn print_output(&self, output: &Path, ids: usize) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "\nTotal unique IDs\t: {}", ids)?;
        writeln!(writer, "Output\t\t\t: {}", output.display())?;

        Ok(())
    }
}

impl Cli for StatsParser<'_> {}

struct StatsParser<'a> {
    matches: &'a ArgMatches<'a>,
    decrement: usize,
}

impl<'a> StatsParser<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            decrement: 0,
        }
    }

    fn show_stats(&mut self) {
        let input_fmt = self.get_input_fmt(&self.matches);
        self.decrement = self.parse_decrement();
        if self.matches.is_present("dir") {
            self.show_stats_dir(&input_fmt);
        } else {
            self.show_stats_file(&input_fmt);
        }
    }

    fn show_stats_dir(&self, input_fmt: &SeqFormat) {
        let dir = self.get_dir_input(&self.matches);
        let files = self.get_files(dir, input_fmt);
        let output = self.get_output(&self.matches);
        self.print_input_file(Path::new(dir)).unwrap();
        SeqStats::new(input_fmt, output, self.decrement).get_stats_dir(&files);
    }

    fn show_stats_file(&self, input_fmt: &SeqFormat) {
        self.get_input_fmt(&self.matches);
        let input = Path::new(self.get_file_input(self.matches));
        let output = self.get_output(&self.matches);
        self.print_input_file(input).unwrap();
        SeqStats::new(input_fmt, output, self.decrement).get_seq_stats_file(input);
    }

    fn parse_decrement(&self) -> usize {
        let decrement = self
            .matches
            .value_of("decrement")
            .expect("CAN'T GET DECREMENT VALUES");
        decrement
            .parse::<usize>()
            .expect("FAIL PARSING DECREMENT VALUES")
    }

    fn print_input_file(&self, input: &Path) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul summary")?;
        writeln!(writer, "Input\t\t: {}\n", input.display())?;

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;

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
    fn get_id_output_path_test() {
        let arg = App::new("segul-test")
            .arg(Arg::with_name("dir").default_value("./test_dir/"))
            .get_matches();
        let id = IdParser::new(&arg);
        let res = PathBuf::from("./test_dir.txt");
        assert_eq!(res, id.get_output_path(&arg));
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

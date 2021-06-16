use std::io::{self, Result, Write};
use std::path::{Path, PathBuf};

use clap::{App, AppSettings, Arg, ArgMatches};
use glob::glob;
use rayon::prelude::*;

use crate::common::{InputFormat, OutputFormat, PartitionFormat};
use crate::converter::Converter;
use crate::msa;
use crate::utils;

fn get_args(version: &str) -> ArgMatches {
    App::new("segul")
        .version(version)
        .about("A genomic sequence tool")
        .author("Heru Handika")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .subcommand(
            App::new("convert")
                .about("Convert sequence formats")
                .subcommand(
                    App::new("fasta")
                        .about("Convert fasta files")
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
                            Arg::with_name("phylip")
                                .long("phylip")
                                .help("Convert fasta to phylip. Default: nexus.")
                                .takes_value(false),
                        ),
                )
                .subcommand(
                    App::new("nexus")
                        .about("Convert nexus files")
                        .arg(
                            Arg::with_name("input")
                                .short("i")
                                .long("input")
                                .help("Convert a nexus file")
                                .takes_value(true)
                                .required_unless("dir")
                                .conflicts_with("dir")
                                .value_name("INPUT FILE"),
                        )
                        .arg(
                            Arg::with_name("dir")
                                .short("d")
                                .long("dir")
                                .help("Convert multiple nexus files inside a dir")
                                .takes_value(true)
                                .required_unless("input")
                                .conflicts_with("input")
                                .value_name("DIR"),
                        )
                        .arg(
                            Arg::with_name("phylip")
                                .long("phylip")
                                .help("Convert nexus to phylip. Default: fasta")
                                .takes_value(false),
                        )
                        .arg(
                            Arg::with_name("output")
                                .short("o")
                                .long("output")
                                .help("Uses a costume output filename")
                                .takes_value(true)
                                .required_unless("input")
                                .value_name("OUTPUT"),
                        ),
                )
                .subcommand(
                    App::new("phylip")
                        .about("Convert phylip files")
                        .arg(
                            Arg::with_name("input")
                                .short("i")
                                .long("input")
                                .help("Convert a phylip file")
                                .takes_value(true)
                                .required_unless("dir")
                                .conflicts_with("dir")
                                .value_name("INPUT FILE"),
                        )
                        .arg(
                            Arg::with_name("dir")
                                .short("d")
                                .long("dir")
                                .help("Convert multiple phylip files inside a dir")
                                .takes_value(true)
                                .required_unless("input")
                                .conflicts_with("input")
                                .value_name("DIR"),
                        )
                        .arg(
                            Arg::with_name("nexus")
                                .long("phylip")
                                .help("Convert nexus to phylip. Default: fasta")
                                .takes_value(false),
                        ),
                ),
        )
        .subcommand(
            App::new("concat")
                .about("Concat alignments")
                .subcommand(
                    App::new("fasta")
                        .about("Concats fasta alignments")
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
                                .help("Sets output format. Choices: nexus, fasta, phylip")
                                .takes_value(true)
                                .default_value("nexus")
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("partition")
                                .short("-p")
                                .long("part")
                                .help("Sets partition format. Choice: nexus, phylip, nexsep")
                                .takes_value(true)
                                .required(true)
                                .default_value("nexsep")
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("output")
                                .short("o")
                                .long("output")
                                .help("Uses a costume output filename. Default: concat.")
                                .takes_value(true)
                                .required(true)
                                .default_value("concat")
                                .value_name("OUTPUT"),
                        ),
                )
                .subcommand(
                    App::new("nexus")
                        .about("Concats nexus alignments")
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
                                .help("Sets output format. Choices: nexus, fasta, phylip")
                                .takes_value(true)
                                .default_value("nexus")
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("partition")
                                .short("-p")
                                .long("part")
                                .help("Sets partition format. Choice: nexus, phylip, nexsep")
                                .takes_value(true)
                                .required(true)
                                .default_value("nexsep")
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("output")
                                .short("o")
                                .long("output")
                                .help("Uses a costume output filename. Default: concat.")
                                .takes_value(true)
                                .required(true)
                                .default_value("concat")
                                .value_name("OUTPUT"),
                        ),
                )
                .subcommand(
                    App::new("phylip")
                        .about("Concats phylip alignments")
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
                                .help("Sets output format. Choices: nexus, fasta, phylip")
                                .takes_value(true)
                                .default_value("nexus")
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("partition")
                                .short("-p")
                                .long("part")
                                .help("Sets partition format. Choice: nexus, raxml, nexsep")
                                .takes_value(true)
                                .required(true)
                                .default_value("nexsep") 
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("output")
                                .short("o")
                                .long("output")
                                .help("Uses a costume output filename. Default: concat.")
                                .takes_value(true)
                                .required(true)
                                .default_value("concat")
                                .value_name("OUTPUT"),
                        ),
                ),
        )
        .get_matches()
}

pub fn parse_cli(version: &str) {
    let args = get_args(version);
    let text = format!("SEGUL v{}", version);
    utils::print_divider(&text, 50);
    match args.subcommand() {
        ("convert", Some(convert_matches)) => parse_convert_subcommand(convert_matches),
        ("concat", Some(concat_matches)) => parse_concat_subcommand(concat_matches),
        _ => unreachable!(),
    }
}

fn parse_convert_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("fasta", Some(fasta_matches)) => {
            ConvertParser::new(fasta_matches, InputFormat::Fasta).convert()
        }
        ("nexus", Some(nexus_matches)) => {
            ConvertParser::new(nexus_matches, InputFormat::Nexus).convert()
        }
        ("phylip", Some(phylip_matches)) => {
            ConvertParser::new(phylip_matches, InputFormat::Phylip).convert()
        }
        _ => unreachable!(),
    }
}

fn parse_concat_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("fasta", Some(fasta_matches)) => {
            ConcatParser::new(fasta_matches, InputFormat::Fasta).concat()
        }
        ("nexus", Some(nexus_matches)) => {
            ConcatParser::new(nexus_matches, InputFormat::Nexus).concat()
        }
        ("phylip", Some(phylip_matches)) => {
            ConcatParser::new(phylip_matches, InputFormat::Phylip).concat()
        }
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

    fn get_files(&self, pattern: &str) -> Vec<PathBuf> {
        glob(pattern)
            .expect("COULD NOT FIND FILES")
            .filter_map(|ok| ok.ok())
            .collect()
    }

    fn get_output<'a>(&self, matches: &'a ArgMatches) -> &'a str {
        matches.value_of("output").expect("CANNOT READ OUTPUT PATH")
    }

    fn set_output(&self, matches: &ArgMatches) -> PathBuf {
        if matches.is_present("output") {
            let output = self.get_output(matches);
            PathBuf::from(output)
        } else {
            PathBuf::from(".")
        }
    }

    fn get_partition_format(&self, matches: &ArgMatches) -> PartitionFormat {
        let part_format = matches
            .value_of("partition")
            .expect("CANNOT READ PARTITION FORMAT");
        match part_format {
            "nexus" => PartitionFormat::Nexus,
            "raxml" => PartitionFormat::Raxml,
            "nexsep" => PartitionFormat::NexusSeparate,
            _ => PartitionFormat::Nexus,
        }
    }

    fn get_output_format(&self, matches: &ArgMatches) -> OutputFormat {
        let format = matches
            .value_of("format")
            .expect("CANNOT READ FORMAT INPUT");
        match format {
            "nexus" => OutputFormat::Nexus,
            "phylip" => OutputFormat::Phylip,
            "fasta" => OutputFormat::Fasta,
            _ => panic!(
                "UNSUPPORTED FORMAT. \
        THE PROGRAM ONLY ACCEPT nexus, phylip, and fasta. All in lowercase.\
        Your input: {} ",
                format
            ),
        }
    }
}

impl Cli for ConvertParser<'_> {}
impl Cli for ConcatParser<'_> {}

struct ConvertParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_format: InputFormat,
    output: PathBuf,
    is_dir: bool,
}

impl<'a> ConvertParser<'a> {
    fn new(matches: &'a ArgMatches<'a>, input_format: InputFormat) -> Self {
        Self {
            matches,
            input_format,
            output: PathBuf::new(),
            is_dir: false,
        }
    }

    fn convert(&mut self) {
        if self.matches.is_present("input") {
            self.convert_file();
        } else {
            self.convert_multiple_fasta();
        }
    }

    fn convert_file(&mut self) {
        let input = Path::new(self.get_file_input(self.matches));
        self.output = self.set_output(self.matches);
        self.display_input_file(input).unwrap();
        self.convert_any(input);
    }

    fn convert_multiple_fasta(&mut self) {
        let dir = self.get_dir_input(self.matches);
        let pattern = self.get_pattern(dir);
        let files = self.get_files(&pattern);
        self.output = self.set_output(&self.matches);
        self.is_dir = true;
        self.display_input_dir(Path::new(dir), files.len()).unwrap();
        files.par_iter().for_each(|file| {
            self.convert_any(file);
        });
    }

    fn convert_any(&self, input: &Path) {
        match self.input_format {
            InputFormat::Fasta => self.convert_fasta(input),
            InputFormat::Nexus => self.convert_nexus(input),
            InputFormat::Phylip => self.convert_phylip(input),
        }
    }

    fn convert_fasta(&self, input: &Path) {
        if self.matches.is_present("phylip") {
            Converter::new(input, &self.output, &OutputFormat::Phylip, self.is_dir).convert_fasta();
        } else {
            Converter::new(input, &self.output, &OutputFormat::Nexus, self.is_dir).convert_fasta();
        }
    }

    fn convert_nexus(&self, input: &Path) {
        if self.matches.is_present("phylip") {
            Converter::new(input, &self.output, &OutputFormat::Phylip, self.is_dir).convert_nexus();
        } else {
            Converter::new(input, &self.output, &OutputFormat::Fasta, self.is_dir).convert_nexus();
        }
    }

    fn convert_phylip(&self, input: &Path) {
        if self.matches.is_present("nexus") {
            Converter::new(input, &self.output, &OutputFormat::Nexus, self.is_dir).convert_phylip();
        } else {
            Converter::new(input, &self.output, &OutputFormat::Fasta, self.is_dir).convert_phylip();
        }
    }

    fn get_pattern(&self, dir: &str) -> String {
        match self.input_format {
            InputFormat::Fasta => format!("{}/*.fa*", dir),
            InputFormat::Nexus => format!("{}/*.nex*", dir),
            InputFormat::Phylip => format!("{}/*.phy*", dir),
        }
    }

    fn display_input_file(&self, input: &Path) -> Result<()> {
        let io = io::stdout();
        let mut writer = io::BufWriter::new(io);
        writeln!(writer, "Command\t: segul convert")?;
        writeln!(writer, "Input\t: {}", &input.display())?;
        Ok(())
    }

    fn display_input_dir(&self, input: &Path, nfile: usize) -> Result<()> {
        let io = io::stdout();
        let mut writer = io::BufWriter::new(io);
        writeln!(writer, "Command\t: segul convert")?;
        writeln!(writer, "Input dir\t: {}", &input.display())?;
        writeln!(writer, "Total files\t: {}", nfile)?;
        writeln!(writer, "Output dir\t: {}", self.output.display())?;
        Ok(())
    }
}

struct ConcatParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_format: InputFormat,
}

impl<'a> ConcatParser<'a> {
    fn new(matches: &'a ArgMatches<'a>, input_format: InputFormat) -> Self {
        Self {
            matches,
            input_format,
        }
    }

    fn concat(&self) {
        let dir = self.get_dir_input(self.matches);
        let output = self.get_output(self.matches);
        let output_format = self.get_output_format(self.matches);
        let part_format = self.get_partition_format(self.matches);
        let mut concat = msa::MSAlignment::new(dir, output, output_format, part_format);
        self.concat_any(&mut concat)
    }

    fn concat_any(&self, concat: &mut msa::MSAlignment) {
        match self.input_format {
            InputFormat::Fasta => concat.concat_fasta(),
            InputFormat::Nexus => concat.concat_nexus(),
            InputFormat::Phylip => concat.concat_phylip(),
        }
    }
}

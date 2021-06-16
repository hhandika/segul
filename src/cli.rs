use std::path::{Path, PathBuf};

use clap::{App, AppSettings, Arg, ArgMatches};
use glob::glob;

use crate::common::{OutputFormat, PartitionFormat};
use crate::converter::Converter;
use crate::msa;

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
    match args.subcommand() {
        ("convert", Some(convert_matches)) => parse_convert_subcommand(convert_matches),
        ("concat", Some(concat_matches)) => parse_concat_subcommand(concat_matches),
        _ => unreachable!(),
    }
}

fn parse_convert_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("fasta", Some(fasta_matches)) => Fasta::new(fasta_matches).convert_fasta(),
        ("nexus", Some(nexus_matches)) => Nexus::new(nexus_matches).convert_nexus(),
        ("phylip", Some(phylip_matches)) => Phylip::new(phylip_matches).convert_phylip(),
        _ => unreachable!(),
    }
}

fn parse_concat_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("fasta", Some(fasta_matches)) => Fasta::new(fasta_matches).concat_fasta(),
        ("nexus", Some(nexus_matches)) => Nexus::new(nexus_matches).concat_nexus(),
        ("phylip", Some(phylip_matches)) => Phylip::new(phylip_matches).concat_phylip(),
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

impl Cli for Fasta<'_> {}
impl Cli for Nexus<'_> {}
impl Cli for Phylip<'_> {}

struct Fasta<'a> {
    matches: &'a ArgMatches<'a>,
    output: PathBuf,
    is_dir: bool,
}

impl<'a> Fasta<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            output: PathBuf::new(),
            is_dir: false,
        }
    }

    fn convert_fasta(&mut self) {
        if self.matches.is_present("input") {
            self.convert_single_fasta();
        } else {
            self.convert_multiple_fasta();
        }
    }

    fn convert_single_fasta(&mut self) {
        let input = Path::new(self.get_file_input(self.matches));
        self.output = self.set_output(self.matches);
        self.convert(input);
    }

    fn convert_multiple_fasta(&mut self) {
        let dir = self.get_dir_input(self.matches);
        let pattern = format!("{}/*.fa*", dir);
        let files = self.get_files(&pattern);
        self.is_dir = true;
        files.iter().for_each(|file| {
            self.output = self.set_output(&self.matches);
            self.convert(file);
        })
    }

    fn convert(&self, input: &Path) {
        if self.matches.is_present("phylip") {
            Converter::new(input, &self.output, &OutputFormat::Phylip, self.is_dir).convert_fasta();
        } else {
            Converter::new(input, &self.output, &OutputFormat::Nexus, self.is_dir).convert_fasta();
        }
    }

    fn concat_fasta(&self) {
        let dir = self.get_dir_input(self.matches);
        let output = self.get_output(self.matches);
        let output_format = self.get_output_format(self.matches);
        let part_format = self.get_partition_format(self.matches);

        msa::MSAlignment::new(dir, output, output_format, part_format).concat_fasta();
    }
}

struct Nexus<'a> {
    matches: &'a ArgMatches<'a>,
    output: PathBuf,
    is_dir: bool,
}

impl<'a> Nexus<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            output: PathBuf::new(),
            is_dir: false,
        }
    }

    fn convert_nexus(&mut self) {
        if self.matches.is_present("input") {
            self.convert_single_nexus();
        } else {
            self.convert_multiple_nexus();
        }
    }

    fn convert_single_nexus(&mut self) {
        let input = Path::new(self.get_file_input(self.matches));
        self.output = self.set_output(&self.matches);
        self.convert(input);
    }

    fn convert_multiple_nexus(&mut self) {
        let dir = self.get_dir_input(self.matches);
        let pattern = format!("{}/*.nex*", dir);
        let files = self.get_files(&pattern);
        self.is_dir = true;
        files.iter().for_each(|file| {
            self.output = self.set_output(&self.matches);
            self.convert(file);
        })
    }

    fn convert(&self, input: &Path) {
        if self.matches.is_present("phylip") {
            Converter::new(input, &self.output, &OutputFormat::Phylip, self.is_dir).convert_nexus();
        } else {
            Converter::new(input, &self.output, &OutputFormat::Fasta, self.is_dir).convert_nexus();
        }
    }

    fn concat_nexus(&self) {
        let dir = self.get_dir_input(self.matches);
        let output = self.get_output(self.matches);
        let output_format = self.get_output_format(self.matches);
        let part_format = self.get_partition_format(self.matches);

        msa::MSAlignment::new(dir, output, output_format, part_format).concat_nexus();
    }
}

struct Phylip<'a> {
    matches: &'a ArgMatches<'a>,
    output: PathBuf,
}

impl<'a> Phylip<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            output: PathBuf::new(),
        }
    }

    fn convert_phylip(&self) {
        let input = Path::new(self.get_file_input(self.matches));
        if self.matches.is_present("nexus") {
            Converter::new(&input, &self.output, &OutputFormat::Nexus, false).convert_phylip();
        } else {
            Converter::new(&input, &self.output, &OutputFormat::Fasta, false).convert_phylip();
        }
    }

    fn concat_phylip(&self) {
        let dir = self.get_dir_input(self.matches);
        let output = self.get_output(self.matches);
        let output_format = self.get_output_format(self.matches);
        let part_format = self.get_partition_format(self.matches);

        msa::MSAlignment::new(dir, output, output_format, part_format).concat_phylip();
    }
}

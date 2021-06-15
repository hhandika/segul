use clap::{App, AppSettings, Arg, ArgMatches};

use crate::common::{SeqFormat, SeqPartition};
use crate::fasta;
use crate::msa;
use crate::nexus;
use crate::phylip;

fn get_args(version: &str) -> ArgMatches {
    App::new("segul")
        .version(version)
        .about("A genomic sequence tool")
        .author("Heru Handika")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .subcommand(
            App::new("fasta")
                .about("Fasta tools")
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Uses a costume output filename")
                        .takes_value(true)
                        .required(true)
                        .global(true)
                        .default_value("concat")
                        .value_name("OUTPUT"),
                )
                .subcommand(
                    App::new("convert")
                        .about("Convert fasta")
                        .arg(
                            Arg::with_name("input")
                                .short("i")
                                .long("input")
                                .help("Inputs a fasta file path")
                                .takes_value(true)
                                .required(true)
                                .value_name("INPUT FILE"),
                        )
                        .arg(
                            Arg::with_name("phylip")
                                .long("phylip")
                                .help("Convert fasta to phylip")
                                .takes_value(false),
                        ),
                )
                .subcommand(
                    App::new("concat")
                        .about("Concat fasta alignments")
                        .arg(
                            Arg::with_name("dir")
                                .short("d")
                                .long("dir")
                                .help("Inputs dir containing fasta files")
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
                                .required(true)
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
                        ),
                ),
        )
        .subcommand(
            App::new("nexus")
                .about("Nexus Tools")
                .subcommand(
                    App::new("convert")
                        .about("Convert nexus")
                        .arg(
                            Arg::with_name("input")
                                .short("i")
                                .long("input")
                                .help("Inputs a nexus file path")
                                .takes_value(true)
                                .required(true)
                                .value_name("INPUT FILE"),
                        )
                        .arg(
                            Arg::with_name("phylip")
                                .long("phylip")
                                .help("Convert nexus to phylip")
                                .takes_value(false),
                        ),
                )
                .subcommand(
                    App::new("concat")
                        .about("Concat nexus alignments")
                        .arg(
                            Arg::with_name("dir")
                                .short("d")
                                .long("dir")
                                .help("Inputs a dir containing nexus files")
                                .takes_value(true)
                                .required(true)
                                .value_name("DIR"),
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
                            Arg::with_name("format")
                                .short("f")
                                .long("format")
                                .help("Sets an output format. Default to nexus. Choices: nexus, fasta, phylip.")
                                .takes_value(true)
                                .required(true)
                                .default_value("nexus")
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("partition")
                                .short("-p")
                                .long("part")
                                .help("Sets a partition format. Default to nexus separate (nexsep). Choice: nexus, phylip, nexsep.")
                                .takes_value(true)
                                .required(true)
                                .default_value("nexsep")
                                .value_name("FORMAT"),
                        ),
                ),
        )
        .subcommand(
            App::new("phylip").about("Any phylip tools").subcommand(
                App::new("convert")
                    .about("Convert phylip")
                    .arg(
                        Arg::with_name("input")
                            .short("i")
                            .long("input")
                            .help("Inputs file path")
                            .takes_value(true)
                            .required(true)
                            .value_name("INPUT FILE"),
                    )
                    .arg(
                        Arg::with_name("nexus")
                            .long("nexus")
                            .help("Convert phylip to nexus")
                            .takes_value(false),
                    ),
            ),
        )
        .get_matches()
}

pub fn parse_cli(version: &str) {
    let args = get_args(version);
    match args.subcommand() {
        ("fasta", Some(fasta_matches)) => parse_fasta_subcommand(fasta_matches),
        ("nexus", Some(nexus_matches)) => parse_nexus_subcommand(nexus_matches),
        ("phylip", Some(phylip_matches)) => parse_phylip_subcommand(phylip_matches),
        _ => unreachable!(),
    }
}

fn parse_fasta_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("convert", Some(convert_matches)) => Fasta::new(convert_matches).convert_fasta(),
        ("concat", Some(concat_matches)) => Fasta::new(concat_matches).concat_fasta(),
        _ => unreachable!(),
    }
}

fn parse_nexus_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("convert", Some(convert_matches)) => Nexus::new(convert_matches).convert_nexus(),
        ("concat", Some(concat_matches)) => Nexus::new(concat_matches).concat_nexus(),
        _ => unreachable!(),
    }
}

fn parse_phylip_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("convert", Some(convert_matches)) => Phylip::new(convert_matches).convert_phylip(),
        ("concat", Some(concat_matches)) => Phylip::new(concat_matches).concat_phylip(),
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

    fn get_output<'a>(&self, matches: &'a ArgMatches) -> &'a str {
        matches.value_of("output").expect("CANNOT READ OUTPUT PATH")
    }

    fn get_partition_format(&self, matches: &ArgMatches) -> SeqPartition {
        let part_format = matches
            .value_of("partition")
            .expect("CANNOT READ PARTITION FORMAT");
        match part_format {
            "nexus" => SeqPartition::Nexus,
            "phylip" => SeqPartition::Phylip,
            "nexsep" => SeqPartition::NexusSeparate,
            _ => SeqPartition::Nexus,
        }
    }

    fn get_output_format(&self, matches: &ArgMatches) -> SeqFormat {
        let format = matches
            .value_of("format")
            .expect("CANNOT READ FORMAT INPUT");
        match format {
            "nexus" => SeqFormat::Nexus,
            "phylip" => SeqFormat::Phylip,
            "fasta" => SeqFormat::Fasta,
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
}

impl<'a> Fasta<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self { matches }
    }

    fn convert_fasta(&self) {
        let input = self.get_file_input(self.matches);
        let output = self.get_output(self.matches);
        if self.matches.is_present("phylip") {
            fasta::convert_fasta(input, output, SeqFormat::Phylip);
        } else {
            fasta::convert_fasta(input, output, SeqFormat::Nexus);
        }
    }

    fn concat_fasta(&self) {
        let dir = self.get_dir_input(self.matches);
        let output = self.get_output(self.matches);
        let output_format = self.get_output_format(self.matches);
        let part_format = self.get_partition_format(self.matches);

        msa::MSA::new(dir, output, output_format, part_format).concat_fasta();
    }
}

struct Nexus<'a> {
    matches: &'a ArgMatches<'a>,
}

impl<'a> Nexus<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self { matches }
    }

    fn convert_nexus(&self) {
        let input = self.get_file_input(self.matches);
        if self.matches.is_present("phylip") {
            nexus::convert_nexus(input, SeqFormat::Phylip);
        } else {
            nexus::convert_nexus(input, SeqFormat::Fasta);
        }
    }

    fn concat_nexus(&self) {
        let dir = self.get_dir_input(self.matches);
        let output = self.get_output(self.matches);
        let output_format = self.get_output_format(self.matches);
        let part_format = self.get_partition_format(self.matches);

        msa::MSA::new(dir, output, output_format, part_format).concat_nexus();
    }
}

struct Phylip<'a> {
    matches: &'a ArgMatches<'a>,
}

impl<'a> Phylip<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self { matches }
    }

    fn convert_phylip(&self) {
        let input = self.get_file_input(self.matches);
        let output = self.get_output(self.matches);
        if self.matches.is_present("nexus") {
            phylip::convert_phylip(input, SeqFormat::Nexus);
        } else {
            phylip::convert_phylip(input, SeqFormat::Fasta);
        }
    }

    fn concat_phylip(&self) {
        let dir = self.get_dir_input(self.matches);
        let output = self.get_output(self.matches);
        let output_format = self.get_output_format(self.matches);
        let part_format = self.get_partition_format(self.matches);

        msa::MSA::new(dir, output, output_format, part_format).concat_phylip();
    }
}

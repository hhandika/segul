use clap::{App, AppSettings, Arg, ArgMatches};

use crate::alignment;
use crate::common::SeqFormat;
use crate::fasta;
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
                .about("Any fasta tools")
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
                    Arg::with_name("id")
                        .long("id")
                        .takes_value(false)
                        .help("Gets IDs only"),
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
                                .help("Inputs file path")
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
                                .help("Input a dir to the alignment path")
                                .takes_value(true)
                                .required(true)
                                .value_name("DIR"),
                        )
                        .arg(
                            Arg::with_name("output")
                                .short("o")
                                .long("output")
                                .help("Input an output file name")
                                .takes_value(true)
                                .required(true)
                                .default_value("concat")
                                .value_name("OUTPUT"),
                        )
                        .arg(
                            Arg::with_name("format")
                                .long("format")
                                .help("Inputs an output format. Choice: nexus, fasta, phylip.")
                                .takes_value(true)
                                .required(true)
                                .default_value("nexus")
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
        ("fasta", Some(fasta_matches)) => parse_fasta(fasta_matches),
        ("nexus", Some(nexus_matches)) => parse_nexus_subcommand(nexus_matches),
        ("phylip", Some(phylip_matches)) => parse_phylip_subcommand(phylip_matches),
        _ => unreachable!(),
    }
}

fn parse_nexus_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("convert", Some(convert_matches)) => convert_nexus(convert_matches),
        ("concat", Some(concat_matches)) => concat_nexus(concat_matches),
        _ => unreachable!(),
    }
}

fn parse_phylip_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("convert", Some(convert_matches)) => convert_phylip(convert_matches),
        _ => unreachable!(),
    }
}

fn convert_nexus(matches: &ArgMatches) {
    let input = matches
        .value_of("input")
        .expect("CANNOT FIND AN INPUT FILE");
    if matches.is_present("phylip") {
        nexus::convert_nexus(input, SeqFormat::Phylip);
    } else {
        nexus::convert_nexus(input, SeqFormat::Fasta);
    }
}

fn concat_nexus(matches: &ArgMatches) {
    let dir = matches.value_of("dir").expect("CANNOT READ DIR PATH");
    let output = matches.value_of("output").expect("CANNOT READ OUTPUT PATH");
    let format = matches
        .value_of("format")
        .expect("CANNOT READ FORMAT INPUT");
    let filetype = get_file_type(format);

    alignment::concat_nexus(dir, output, filetype);
}

fn get_file_type(format: &str) -> SeqFormat {
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

fn parse_fasta(matches: &ArgMatches) {
    let input = matches
        .value_of("input")
        .expect("CANNOT FIND AN INPUT FILE");
    let id = matches.is_present("id");

    if id {
        fasta::parse_fasta_id(input);
    } else {
        fasta::parse_fasta(input);
    }
}

fn convert_phylip(matches: &ArgMatches) {
    let input = matches
        .value_of("input")
        .expect("CANNOT FIND AN INPUT FILE");
    if matches.is_present("nexus") {
        phylip::convert_phylip(input, SeqFormat::Nexus);
    } else {
        phylip::convert_phylip(input, SeqFormat::Fasta);
    }
}

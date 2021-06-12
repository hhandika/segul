use clap::{App, AppSettings, Arg, ArgMatches};

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
            App::new("nexus").about("Nexus Tools").subcommand(
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
            ),
        )
        .subcommand(
            App::new("phylip")
                .about("Any phylip tools")
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
                    Arg::with_name("convert")
                        .long("convert")
                        .help("Convert nexus to fasta")
                        .takes_value(false),
                ),
        )
        .get_matches()
}

pub fn parse_cli(version: &str) {
    let args = get_args(version);
    match args.subcommand() {
        ("fasta", Some(fasta_matches)) => parse_fasta(fasta_matches),
        ("nexus", Some(nexus_matches)) => parse_nexus_subcommand(nexus_matches),
        ("phylip", Some(phylip_matches)) => parse_phylip(phylip_matches),
        _ => unreachable!(),
    }
}

fn parse_nexus_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("convert", Some(convert_matches)) => convert_nexus(convert_matches),
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

fn parse_phylip(matches: &ArgMatches) {
    let input = matches
        .value_of("input")
        .expect("CANNOT FIND AN INPUT FILE");
    phylip::read_phylip(input);
}

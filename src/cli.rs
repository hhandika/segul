use clap::{App, AppSettings, Arg, ArgMatches};

use crate::fasta;

fn get_args(version: &str) -> ArgMatches {
    App::new("segul")
        .version(version)
        .about("A genomic sequence tool")
        .author("Heru Handika")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .subcommand(
            App::new("fasta").about("Any fasta tools").arg(
                Arg::with_name("input")
                    .short("i")
                    .long("input")
                    .help("Inputs file path")
                    .takes_value(true)
                    .required(true)
                    .value_name("INPUT FILE"),
            ),
        )
        .get_matches()
}

pub fn parse_cli(version: &str) {
    let args = get_args(version);
    match args.subcommand() {
        ("fasta", Some(fasta_matches)) => parse_fasta(fasta_matches),
        _ => unreachable!(),
    }
}

fn parse_fasta(matches: &ArgMatches) {
    let input = matches
        .value_of("input")
        .expect("CANNOT FIND AN INPUT FILE");
    fasta::parse_fasta_id(input);
}

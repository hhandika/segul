mod new_args;

use crate::cli::new_args::{Cli, MainSubcommand};
use clap::Parser;

pub fn parse_cli() {
    let args = Cli::parse();

    match args.subcommand {
        MainSubcommand::RawRead(_) => {
            println!("RawRead");
        }
        MainSubcommand::Contig(_) => {
            println!("Contig");
        }
        MainSubcommand::Alignment(_) => {
            println!("Alignment");
        }
        MainSubcommand::Partition(_) => {
            println!("Partition");
        }
        MainSubcommand::Sequence(_) => {
            println!("Sequence");
        }
    }
}

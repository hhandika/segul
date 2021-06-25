use std::io::{self, Result, Write};
use std::path::{Path, PathBuf};

use clap::{App, AppSettings, Arg, ArgMatches};
use rayon::prelude::*;

use crate::common::{PartitionFormat, SeqFormat};
use crate::converter::Converter;
use crate::finder::Files;
use crate::msa;
use crate::picker::Picker;
use crate::stats::AlnStats;
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
                            Arg::with_name("format")
                                .short("f")
                                .long("format")
                                .help("Sets output format. Choices: nexus, nexus-int, fasta, fasta-int, phylip, phylip-int")
                                .takes_value(true)
                                .default_value("nexus")
                                .value_name("FORMAT"),
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
                            Arg::with_name("format")
                                .short("f")
                                .long("format")
                                .help("Sets output format. Choices: nexus, nexus-int, fasta, fasta-int, phylip, phylip-int")
                                .takes_value(true)
                                .default_value("fasta")
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("output")
                                .short("o")
                                .long("output")
                                .help("Sets target directory or use a costume file name for a single input")
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
                            Arg::with_name("output")
                                .short("o")
                                .long("output")
                                .help("Sets target directory or use a costume file name for a single input")
                                .takes_value(true)
                                .required_unless("input")
                                .value_name("OUTPUT"),
                        )
                        .arg(
                            Arg::with_name("interleave")
                            .long("interleave")
                            .help("Is interleave phylip")
                            .takes_value(false)
                        )
                        .arg(
                            Arg::with_name("format")
                                .short("f")
                                .long("format")
                                .help("Sets output format. Choices: nexus, nexus-int, fasta, fasta-int, phylip, phylip-int")
                                .takes_value(true)
                                .default_value("fasta")
                                .value_name("FORMAT"),
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
                                .help("Sets output format. Choices: nexus, nexus-int, fasta, fasta-int, phylip, phylip-int")
                                .takes_value(true)
                                .default_value("nexus")
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("partition")
                                .short("-p")
                                .long("part")
                                .help("Sets partition format. Choices: nexus, phylip, charset")
                                .takes_value(true)
                                .required(true)
                                .default_value("charset")
                                .value_name("FORMAT"),
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
                            Arg::with_name("codon")
                            .long("codon")
                            .help("Sets codon model partition format")
                            .takes_value(false)
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
                                .help("Sets output format. Choices: nexus, nexus-int, fasta, fasta-int, phylip, phylip-int")
                                .takes_value(true)
                                .default_value("nexus")
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("partition")
                                .short("-p")
                                .long("part")
                                .help("Sets partition format. Choices: nexus, phylip, charset")
                                .takes_value(true)
                                .required(true)
                                .default_value("charset")
                                .value_name("FORMAT"),
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
                            Arg::with_name("codon")
                            .long("codon")
                            .help("Sets codon model partition format")
                            .takes_value(false)
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
                                .help("Sets output format. Choices: nexus, nexus-int, fasta, fasta-int, phylip, phylip-int")
                                .takes_value(true)
                                .default_value("nexus")
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("partition")
                                .short("-p")
                                .long("part")
                                .help("Sets partition format. Choices: nexus, raxml, charset")
                                .takes_value(true)
                                .required(true)
                                .default_value("charset") 
                                .value_name("FORMAT"),
                        )
                        .arg(
                            Arg::with_name("interleave")
                            .long("interleave")
                            .help("Is interleave phylip")
                            .takes_value(false)
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
                            Arg::with_name("codon")
                            .long("codon")
                            .help("Sets codon model partition format")
                            .takes_value(false)
                        ),
                ),
        )
        .subcommand(App::new("pick").about("Gets alignment statistics").arg(
                            Arg::with_name("dir")
                                .short("d")
                                .long("dir")
                                .help("Inputs dir with alignment files")
                                .takes_value(true)
                                .value_name("INPUT FILE"),
                        ))
        .subcommand(App::new("stats").about("Gets alignment statistics").arg(
                            Arg::with_name("input")
                                .short("i")
                                .long("input")
                                .help("Inputs alignment file")
                                .takes_value(true)
                                .required_unless("dir")
                                .conflicts_with("dir")
                                .value_name("INPUT FILE"),
                        ))
        .get_matches()
}

pub fn parse_cli(version: &str) {
    let args = get_args(version);
    let text = format!("SEGUL v{}", version);
    utils::print_divider(&text, 50);
    match args.subcommand() {
        ("convert", Some(convert_matches)) => parse_convert_subcommand(convert_matches),
        ("concat", Some(concat_matches)) => parse_concat_subcommand(concat_matches),
        ("pick", Some(pick_matches)) => {
            PickParser::new(pick_matches, SeqFormat::Nexus).get_min_taxa()
        }
        ("stats", Some(stats_matches)) => {
            StatsParser::new(stats_matches, SeqFormat::Nexus).show_stats()
        }
        _ => unreachable!(),
    }
}

fn parse_convert_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("fasta", Some(fasta_matches)) => {
            ConvertParser::new(fasta_matches, SeqFormat::Fasta).convert()
        }
        ("nexus", Some(nexus_matches)) => {
            ConvertParser::new(nexus_matches, SeqFormat::Nexus).convert()
        }
        ("phylip", Some(phylip_matches)) => {
            ConvertParser::new(phylip_matches, SeqFormat::Phylip).convert()
        }
        _ => unreachable!(),
    }
}

fn parse_concat_subcommand(args: &ArgMatches) {
    match args.subcommand() {
        ("fasta", Some(fasta_matches)) => {
            ConcatParser::new(fasta_matches, SeqFormat::Fasta).concat()
        }
        ("nexus", Some(nexus_matches)) => {
            ConcatParser::new(nexus_matches, SeqFormat::Nexus).concat()
        }
        ("phylip", Some(phylip_matches)) => {
            ConcatParser::new(phylip_matches, SeqFormat::Phylip).concat()
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

    fn get_files(&self, dir: &str, input_format: &SeqFormat) -> Vec<PathBuf> {
        Files::new(dir, input_format).get_files()
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

    fn check_phylip_interleave(&self, matches: &ArgMatches) -> bool {
        matches.is_present("interleave")
    }

    fn get_output_format(&self, matches: &ArgMatches) -> SeqFormat {
        let output_format = matches
            .value_of("format")
            .expect("CANNOT READ FORMAT INPUT");
        match output_format {
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
                output_format
            ),
        }
    }
}

impl Cli for ConvertParser<'_> {}
impl Cli for ConcatParser<'_> {
    fn set_output(&self, matches: &ArgMatches) -> PathBuf {
        let output = self.get_output(matches);
        PathBuf::from(output)
    }
}

impl Cli for PickParser<'_> {}
impl Cli for StatsParser<'_> {}

struct ConvertParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_format: SeqFormat,
    output: PathBuf,
    is_dir: bool,
}

impl<'a> ConvertParser<'a> {
    fn new(matches: &'a ArgMatches<'a>, input_format: SeqFormat) -> Self {
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
        let output_format = self.get_output_format(self.matches);
        self.output = self.set_output(self.matches);
        self.display_input_file(input).unwrap();
        self.convert_any(input, &output_format);
    }

    fn convert_multiple_fasta(&mut self) {
        let dir = self.get_dir_input(self.matches);
        let files = self.get_files(dir, &self.input_format);
        let output_format = self.get_output_format(self.matches);
        self.output = self.set_output(&self.matches);
        self.is_dir = true;
        self.display_input_dir(Path::new(dir), files.len()).unwrap();
        files.par_iter().for_each(|file| {
            self.convert_any(file, &output_format);
        });
    }

    fn convert_any(&self, input: &Path, output_format: &SeqFormat) {
        let mut convert = Converter::new(input, &self.output, output_format, self.is_dir);
        match self.input_format {
            SeqFormat::Fasta => convert.convert_fasta(),
            SeqFormat::Nexus => convert.convert_nexus(),
            SeqFormat::Phylip => {
                let interleave = self.check_phylip_interleave(&self.matches);
                convert.convert_phylip(interleave)
            }
            _ => (),
        }
    }

    fn display_input_file(&self, input: &Path) -> Result<()> {
        let io = io::stdout();
        let mut writer = io::BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul convert")?;
        writeln!(writer, "Input\t\t: {}", &input.display())?;
        Ok(())
    }

    fn display_input_dir(&self, input: &Path, nfile: usize) -> Result<()> {
        let io = io::stdout();
        let mut writer = io::BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul convert")?;
        writeln!(writer, "Input dir\t: {}", &input.display())?;
        writeln!(
            writer,
            "Total files\t: {}",
            utils::format_thousand_sep(&nfile)
        )?;
        writeln!(writer, "Output dir\t: {}", self.output.display())?;
        Ok(())
    }
}

struct ConcatParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_format: SeqFormat,
    part_format: PartitionFormat,
    codon: bool,
}

impl<'a> ConcatParser<'a> {
    fn new(matches: &'a ArgMatches<'a>, input_format: SeqFormat) -> Self {
        Self {
            matches,
            input_format,
            part_format: PartitionFormat::Charset,
            codon: false,
        }
    }

    fn concat(&mut self) {
        let dir = self.get_dir_input(self.matches);
        let output = self.get_output(self.matches);
        let output_format = self.get_output_format(self.matches);
        self.get_partition_format();
        let interleave = self.check_phylip_interleave(&self.matches);
        if interleave {
            self.input_format = SeqFormat::PhylipInt;
        }
        self.display_input_dir(&dir).unwrap();
        let concat =
            msa::MSAlignment::new(&self.input_format, output, output_format, &self.part_format);
        let mut files = self.get_files(dir, &self.input_format);
        concat.concat_alignment(&mut files);
    }

    fn get_partition_format(&mut self) {
        let part_format = self
            .matches
            .value_of("partition")
            .expect("CANNOT READ PARTITION FORMAT");
        if self.matches.is_present("codon") {
            self.codon = true;
            self.get_partition_format_codon(part_format);
        } else {
            self.get_partition_format_std(part_format);
        }
        self.check_partition_format();
    }

    fn get_partition_format_std(&mut self, part_format: &str) {
        self.part_format = match part_format {
            "nexus" => PartitionFormat::Nexus,
            "raxml" => PartitionFormat::Raxml,
            "charset" => PartitionFormat::Charset,
            _ => PartitionFormat::Nexus,
        };
    }

    fn get_partition_format_codon(&mut self, part_format: &str) {
        self.part_format = match part_format {
            "charset" => PartitionFormat::CharsetCodon,
            "nexus" => PartitionFormat::NexusCodon,
            "raxml" => PartitionFormat::RaxmlCodon,
            _ => PartitionFormat::NexusCodon,
        };
    }

    fn display_input_dir(&self, input: &str) -> Result<()> {
        let io = io::stdout();
        let mut writer = io::BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul concat")?;
        writeln!(writer, "Input dir\t: {}\n", input)?;
        Ok(())
    }

    fn check_partition_format(&self) {
        match self.input_format {
            SeqFormat::Nexus => (),
            _ => {
                if let PartitionFormat::Nexus | PartitionFormat::NexusCodon = self.part_format {
                    panic!(
                        "CANNOT WRITE EMBEDDED-NEXUS PARTITION TO NON-NEXUS OUTPUT. \
                MAYBE YOU MEAN TO WRITE THE PARTITION TO 'charset' INSTEAD."
                    )
                }
            }
        }
    }
}

struct PickParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_format: SeqFormat,
    percent: f64,
    output_dir: PathBuf,
}

impl<'a> PickParser<'a> {
    fn new(matches: &'a ArgMatches<'a>, input_format: SeqFormat) -> Self {
        Self {
            matches,
            input_format,
            percent: 0.9,
            output_dir: PathBuf::from("min_taxa_test"),
        }
    }

    fn get_min_taxa(&self) {
        let dir = self.get_dir_input(self.matches);
        let mut files = self.get_files(dir, &self.input_format);
        let pick = Picker::new(
            &mut files,
            &self.input_format,
            &self.output_dir,
            self.percent,
        );
        pick.get_min_taxa();
    }
}

struct StatsParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_format: SeqFormat,
}

impl<'a> StatsParser<'a> {
    fn new(matches: &'a ArgMatches<'a>, input_format: SeqFormat) -> Self {
        Self {
            matches,
            input_format,
        }
    }

    fn show_stats(&mut self) {
        let input = Path::new(self.get_file_input(self.matches));
        let interleave = self.check_phylip_interleave(&self.matches);
        if interleave {
            self.input_format = SeqFormat::PhylipInt;
        }
        AlnStats::new().get_stats(input, &self.input_format);
    }
}

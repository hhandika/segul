use std::io::{self, BufWriter, Result, Write};
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
                        .help("Sets input format. Choices: fasta, nexus, phylip, fasta-int, nexus-int, phylip-int")
                        .takes_value(true)
                        .required(true)
                        .value_name("FORMAT"),
                )
                .arg(
                    Arg::with_name("target")
                        .short("t")
                        .long("target")
                        .help("Sets target output format. Choices: fasta, nexus, phylip, fasta-int, nexus-int, phylip-int")
                        .takes_value(true)
                        .default_value("nexus")
                        .value_name("FORMAT"),
                )
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
                        .help("Sets input format. Choices: fasta, nexus, phylip, fasta-int, nexus-int, phylip-int")
                        .takes_value(true)
                        .required(true)
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
                    Arg::with_name("target")
                        .short("t")
                        .long("target")
                        .help("Sets target output format. Choices: fasta, nexus, phylip, fasta-int, nexus-int, phylip-int")
                        .takes_value(true)
                        .default_value("nexus")
                        .value_name("FORMAT"),
                )
                .arg(
                    Arg::with_name("codon")
                    .long("codon")
                    .help("Sets codon model partition format")
                    .takes_value(false)
                )
        )
        .subcommand(
            App::new("pick")
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
                        .help("Sets input format. Choices: fasta, nexus, phylip, fasta-int, nexus-int, phylip-int")
                        .takes_value(true)
                        .required(true)
                        .value_name("FORMAT"),
                )
                .arg(
                    Arg::with_name("percent")
                        .long("percent")
                        .help("Sets percentage of minimal taxa")
                        .takes_value(true)
                        .required_unless("npercent")
                        .default_value("0.75")
                        .value_name("FORMAT"),
                )
                .arg(
                    Arg::with_name("npercent")
                        .long("npercent")
                        .help("Sets minimal taxa in multiple percentages")
                        .takes_value(true)
                        .multiple(true)
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
                    Arg::with_name("concat")
                    .long("codon")
                    .help("Concats the final results")
                    .takes_value(false)
                ),
        )
        .subcommand(App::new("summary").about("Gets alignment summary stats").arg(
                            Arg::with_name("dir")
                                .short("d")
                                .long("dir")
                                .help("Inputs dir with alignment files")
                                .takes_value(true)
                                .value_name("INPUT FILE"),
                        ))
        .get_matches()
}

pub fn parse_cli(version: &str) {
    let args = get_args(version);
    let text = format!("SEGUL v{}", version);
    utils::print_title(&text);
    match args.subcommand() {
        ("convert", Some(convert_matches)) => ConvertParser::new(convert_matches).convert(),
        ("concat", Some(concat_matches)) => ConcatParser::new(concat_matches).concat(),
        ("pick", Some(pick_matches)) => PickParser::new(pick_matches).get_min_taxa(),
        ("stats", Some(stats_matches)) => StatsParser::new(stats_matches).show_stats(),
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

    fn get_output_path(&self, matches: &ArgMatches) -> PathBuf {
        if matches.is_present("output") {
            let output = self.get_output(matches);
            PathBuf::from(output)
        } else {
            PathBuf::from(".")
        }
    }

    fn get_input_format(&self, matches: &ArgMatches) -> SeqFormat {
        let input_format = matches
            .value_of("format")
            .expect("CANNOT READ FORMAT INPUT");
        match input_format {
            "fasta" | "fasta-int" => SeqFormat::Fasta,
            "nexus" | "nexus-int" => SeqFormat::Nexus,
            "phylip" => SeqFormat::Phylip,
            "phylip-int" => SeqFormat::PhylipInt,
            _ => panic!(
                "UNSUPPORTED FORMAT. \
        THE PROGRAM ONLY ACCEPT fasta, fasta-int, nexus, nexus-int, phylip, and phylip-int. ALL IN lowercase. \
        YOUR INPUT: {} ",
                input_format
            ),
        }
    }

    fn get_output_format(&self, matches: &ArgMatches) -> SeqFormat {
        let output_format = matches
            .value_of("target")
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
    fn get_output_path(&self, matches: &ArgMatches) -> PathBuf {
        PathBuf::from(self.get_output(matches))
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
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_format: SeqFormat::Fasta,
            output: PathBuf::new(),
            is_dir: false,
        }
    }

    fn convert(&mut self) {
        self.input_format = self.get_input_format(&self.matches);
        if self.matches.is_present("input") {
            self.convert_file();
        } else {
            self.convert_multiple_fasta();
        }
    }

    fn convert_file(&mut self) {
        let input = Path::new(self.get_file_input(self.matches));
        let output_format = self.get_output_format(self.matches);
        self.output = self.get_output_path(self.matches);
        self.display_input_file(input).unwrap();
        self.convert_any(input, &output_format);
    }

    fn convert_multiple_fasta(&mut self) {
        let dir = self.get_dir_input(self.matches);
        let files = self.get_files(dir, &self.input_format);
        let output_format = self.get_output_format(self.matches);
        self.output = self.get_output_path(&self.matches);
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
            SeqFormat::Phylip => convert.convert_phylip(false),
            SeqFormat::PhylipInt => convert.convert_phylip(true),
            _ => (),
        }
    }

    fn display_input_file(&self, input: &Path) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul convert")?;
        writeln!(writer, "Input\t\t: {}", &input.display())?;
        Ok(())
    }

    fn display_input_dir(&self, input: &Path, nfile: usize) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "Command\t\t: segul convert")?;
        writeln!(writer, "Input dir\t: {}", &input.display())?;
        writeln!(writer, "Total files\t: {}", utils::fmt_num(&nfile))?;
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
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_format: SeqFormat::Fasta,
            part_format: PartitionFormat::Charset,
            codon: false,
        }
    }

    fn concat(&mut self) {
        self.input_format = self.get_input_format(self.matches);
        let dir = self.get_dir_input(self.matches);
        let output = self.get_output(self.matches);
        let output_format = self.get_output_format(self.matches);
        self.get_partition_format();
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
    output_dir: PathBuf,
}

impl<'a> PickParser<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_format: SeqFormat::Fasta,
            output_dir: PathBuf::new(),
        }
    }

    fn get_min_taxa(&mut self) {
        self.input_format = self.get_input_format(self.matches);
        let dir = self.get_dir_input(self.matches);
        let mut files = self.get_files(dir, &self.input_format);
        if self.is_npercent() {
            self.get_min_taxa_npercent(dir, &mut files);
        } else {
            let percent = self.get_percent();
            self.set_output_path(dir, &percent);
            self.display_input(dir).expect("CANNOT DISPLAY TO STDOUT");
            self.get_min_taxa_percent(&mut files, percent);
        }
    }

    fn get_min_taxa_percent(&mut self, files: &mut [PathBuf], percent: f64) {
        let mut pick = Picker::new(files, &self.input_format, &self.output_dir, percent);
        pick.get_min_taxa();
    }

    fn get_min_taxa_npercent(&mut self, dir: &str, files: &mut [PathBuf]) {
        let npercent = self.get_npercent();
        npercent.iter().for_each(|np| {
            self.set_multi_output_path(dir, np);
            self.display_input(dir).expect("CANNOT DISPLAY TO STDOUT");
            self.get_min_taxa_percent(files, *np);
            utils::print_divider();
        });
    }

    fn get_npercent(&self) -> Vec<f64> {
        let npercent: Vec<&str> = self.matches.values_of("npercent").unwrap().collect();
        npercent
            .iter()
            .map(|np| self.parse_percentage(np))
            .collect()
    }

    fn is_npercent(&mut self) -> bool {
        self.matches.is_present("npercent")
    }

    fn get_percent(&self) -> f64 {
        let percent = self
            .matches
            .value_of("percent")
            .expect("CANNOT GET PERCENTAGE VALUES");
        self.parse_percentage(percent)
    }

    fn parse_percentage(&self, percent: &str) -> f64 {
        percent
            .parse::<f64>()
            .expect("CANNOT PARSE PERCENTAGE VALUES TO FLOATING POINTS")
    }

    fn set_output_path<P: AsRef<Path>>(&mut self, dir: P, percent: &f64) {
        if self.matches.is_present("output") {
            self.output_dir = PathBuf::from(self.get_output(self.matches));
        } else {
            self.output_dir = self.get_formatted_output(dir.as_ref(), percent);
        }
    }

    fn set_multi_output_path<P: AsRef<Path>>(&mut self, dir: P, percent: &f64) {
        if self.matches.is_present("output") {
            let output_dir = PathBuf::from(self.get_output(self.matches));
            self.output_dir = self.get_formatted_output(&output_dir, percent)
        } else {
            self.output_dir = self.get_formatted_output(dir.as_ref(), percent);
        }
    }

    fn get_formatted_output(&self, dir: &Path, percent: &f64) -> PathBuf {
        let parent = dir.parent().unwrap();
        let last = dir.file_name().unwrap().to_string_lossy();
        let output_dir = format!("{}_{}p", last, percent * 100.0);
        parent.join(output_dir)
    }

    fn display_input(&self, dir: &str) -> Result<()> {
        let io = io::stdout();
        let mut writer = BufWriter::new(io);
        writeln!(writer, "\x1b[0;33mInput\x1b[0m")?;
        writeln!(writer, "Dir\t\t: {}", dir)?;

        Ok(())
    }
}

struct StatsParser<'a> {
    matches: &'a ArgMatches<'a>,
    input_format: SeqFormat,
}

impl<'a> StatsParser<'a> {
    fn new(matches: &'a ArgMatches<'a>) -> Self {
        Self {
            matches,
            input_format: SeqFormat::Nexus,
        }
    }

    fn show_stats(&mut self) {
        self.get_input_format(&self.matches);
        let input = Path::new(self.get_file_input(self.matches));
        AlnStats::new().get_stats(input, &self.input_format);
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
        let min_taxa = PickParser::new(&arg);
        let dir = "./test_taxa/";
        let percent = 0.75;
        let res = PathBuf::from("./test_taxa_75p");
        let output = min_taxa.get_formatted_output(Path::new(dir), &percent);
        assert_eq!(res, output);
    }
}

use clap::{App, AppSettings, Arg, ArgMatches};

pub fn get_args(version: &str) -> ArgMatches {
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
                        .conflicts_with_all(&["dir", "wildcard"])
                        .value_name("INPUT FILE"),
                )
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Convert multiple fasta files inside a dir")
                        .takes_value(true)
                        .required_unless("input")
                        .conflicts_with_all(&["input", "wildcard"])
                        .value_name("DIR"),
                )
                .arg(
                    Arg::with_name("wildcard")
                        .short("c")
                        .long("wcard")
                        .help("Convert multiple fasta files using wildcard as an input")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless("input")
                        .conflicts_with_all(&["input", "dir"])
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
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
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
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
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
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
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
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
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
                        .help("Gets summary from a file")
                        .takes_value(true)
                        .required_unless("dir")
                        .conflicts_with_all(&["dir", "wildcard"])
                        .value_name("INPUT FILE"),
                )
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Gets summary from alignment files")
                        .takes_value(true)
                        .conflicts_with_all(&["input", "wildcard"])
                        .value_name("INPUT FILE"),
                )
                .arg(
                    Arg::with_name("wildcard")
                        .short("c")
                        .long("wcard")
                        .help("Convert multiple fasta files using wildcard as an input")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless("input")
                        .conflicts_with_all(&["input", "dir"])
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
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
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

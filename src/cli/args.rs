use clap::{crate_description, crate_name, App, AppSettings, Arg, ArgMatches};

pub fn get_args(version: &str) -> ArgMatches {
    App::new(crate_name!())
        .version(version)
        .about(crate_description!())
        .author("Heru Handika")
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .subcommand(
            App::new("convert")
                .about("Converts sequence formats")
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Input a directory")
                        .takes_value(true)
                        .required_unless("input")
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .long("input")
                        .help("Input path with wildcard support")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Sets target directory or use a costume file name for a single input")
                        .takes_value(true)
                        .required(true)
                        .default_value("SEGUL-Convert")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("input-format")
                        .short("f")
                        .long("input-format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
                            "fasta",
                            "nexus",
                            "phylip",
                        ]),
                )
                .arg(
                    Arg::with_name("output-format")
                        .short("F")
                        .long("output-format")
                        .help("Sets target output format")
                        .takes_value(true)
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
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
                    Arg::with_name("datatype")
                        .long("datatype")
                        .help("Sets data type")
                        .takes_value(true)
                        .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
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
                        .required_unless("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .long("input")
                        .help("Input path with wildcard support")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::with_name("input-format")
                        .short("f")
                        .long("input-format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
                            "nexus",
                            "phylip",
                            "fasta",
                        ]),
                )
                .arg(
                    Arg::with_name("partition")
                        .short("-p")
                        .long("part")
                        .help("Sets partition format")
                        .takes_value(true)
                        .required(true)
                        .default_value("nexus")
                        .value_name("PART-FORMAT")
                        .possible_values(&["charset", "nexus", "raxml"]),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Uses a costume output directory")
                        .takes_value(true)
                        .required(true)
                        .default_value("SEGUL-Concat")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("output-format")
                        .short("F")
                        .long("output-format")
                        .help("Sets target output format")
                        .takes_value(true)
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
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
                    Arg::with_name("prefix")
                        .long("prefix")
                        .help("Uses a costume output filename")
                        .takes_value(true)
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("datatype")
                        .long("datatype")
                        .help("Sets data type")
                        .takes_value(true)
                        .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
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
                .about("Filter alignments based on min taxon completeness, alignment length, and parsimony informative sites")
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Inputs a path to alignment dir")
                        .takes_value(true)
                        .required_unless("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .long("input")
                        .help("Input path with wildcard support")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::with_name("input-format")
                        .short("f")
                        .long("input-format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
                            "fasta",
                            "nexus",
                            "phylip",
                        ]),
                )
                .arg(
                    Arg::with_name("datatype")
                        .long("datatype")
                        .help("Sets data type")
                        .takes_value(true)
                        .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::with_name("percent")
                        .long("percent")
                        .help("Sets percentage of minimal taxa")
                        .takes_value(true)
                        .required_unless_all(&["npercent", "aln-len", "pars-inf", "percent-inf"])
                        .conflicts_with_all(&["npercent", "aln-len", "pars-inf", "percent-inf"])
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::with_name("npercent")
                        .long("npercent")
                        .help("Inputs percentages of minimal taxa (allow multiple values)")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "aln-len", "pars-inf"])
                        .multiple(true)
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::with_name("ntax")
                        .long("ntax")
                        .help("Inputs the total number of taxa")
                        .takes_value(true)
                        .conflicts_with_all(&["aln-len", "pars-inf"])
                        .value_name("INTEGER"),
                )
                .arg(
                    Arg::with_name("aln-len")
                        .long("len")
                        .help("Inputs minimal alignment length")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "npercent", "pars-inf"])
                        .value_name("INTEGER"),
                )
                .arg(
                    Arg::with_name("pars-inf")
                        .long("pinf")
                        .help("Inputs minimal parsimony informative sites")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "npercent", "aln-len"])
                        .value_name("INTEGER"),
                )
                .arg(
                    Arg::with_name("percent-inf")
                        .long("percent-inf")
                        .help("Inputs percent parsimony informative sites")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "npercent", "aln-len", "pars-inf"])
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Sets an output directory")
                        .takes_value(true)
                        .required_unless("dir")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("output-format")
                        .short("F")
                        .long("output-format")
                        .help("Sets output format if concat")
                        .takes_value(true)
                        .value_name("SEQ-FORMAT")
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
                    Arg::with_name("prefix")
                        .long("prefix")
                        .help("Specifies prefix for output filename")
                        .takes_value(true)
                        .requires("concat")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("partition")
                        .short("-p")
                        .long("part")
                        .help("Sets partition format if concat")
                        .takes_value(true)
                        .requires("concat")
                        .value_name("PART-FORMAT")
                        .possible_values(&["charset", "nexus", "raxml"]),
                )
                .arg(
                    Arg::with_name("concat")
                        .long("concat")
                        .help("Concats the final results")
                        .requires("partition")
                        .takes_value(false),
                )
                .arg(
                    Arg::with_name("codon")
                        .long("codon")
                        .requires("concat")
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
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .long("input")
                        .help("Input path with wildcard support")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless("input")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::with_name("input-format")
                        .short("f")
                        .long("input-format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
                            "fasta",
                            "nexus",
                            "phylip",
                        ]),
                )
                .arg(
                    Arg::with_name("datatype")
                        .long("datatype")
                        .help("Sets data type")
                        .takes_value(true)
                        .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Uses a costume output filename")
                        .takes_value(true)
                        .required(true)
                        .default_value("id")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("map")
                        .short("m")
                        .long("map")
                        .help("Map ID across all alignments")
                        .takes_value(false),
                ),
        )
        .subcommand(
            App::new("rename")
            .about("Batch renaming sequence IDs in multiple alignments")
                .arg(
                    Arg::with_name("dir")
                    .short("d")
                    .long("dir")
                    .help("Inputs dir with alignment files")
                    .takes_value(true)
                    .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("input")
                            .short("i")
                            .long("input")
                            .help("Input path with wildcard support")
                            .takes_value(true)
                            .multiple(true)
                            .required_unless("input")
                            .conflicts_with_all(&["input", "dir"])
                            .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::with_name("input-format")
                        .short("f")
                        .long("input-format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
                            "fasta",
                            "nexus",
                            "phylip",
                        ]),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Uses a costume output filename")
                        .takes_value(true)
                        .required(true)
                        .default_value("SEGUL-Rename")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("output-format")
                        .short("F")
                        .long("output-format")
                        .help("Sets target output format")
                        .takes_value(true)
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
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
                    Arg::with_name("datatype")
                        .long("datatype")
                        .help("Sets data type")
                        .takes_value(true)
                        .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::with_name("names")
                        .short("n")
                        .long("names")
                        .help("Inputs IDs to rename")
                        .takes_value(true)
                        .required(true)
                        .value_name("PATH"),

                )
                .arg(
                    Arg::with_name("dry-run")
                        .long("dry-run")
                        .help("Check if the program can parse the input ids correctly")
                        .takes_value(false),
                ),
        )
        .subcommand(
            App::new("summary")
                .about("Gets alignment summary stats")
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Gets summary from alignment files")
                        .takes_value(true)
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .long("input")
                        .help("Input path with wildcard support")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::with_name("input-format")
                        .short("f")
                        .long("input-format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
                            "fasta",
                            "nexus",
                            "phylip",
                        ]),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Uses a costume output filename")
                        .takes_value(true)
                        .required(true)
                        .default_value("SEGUL-Stats")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("datatype")
                        .long("datatype")
                        .help("Sets data type")
                        .takes_value(true)
                        .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::with_name("percent-interval")
                        .long("interval")
                        .help("Sets a custom percentage interval value for counting data matrix completeness")
                        .takes_value(true)
                        .value_name("INTEGER")
                        .default_value("5")
                        .possible_values(&["1", "2", "5", "10"]),
                ),
        )
        .subcommand(
            App::new("extract")
                .about("Extract sequences from a collection of alignments")
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Inputs a directory path to alignments")
                        .takes_value(true)
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .long("input")
                        .help("Input path with wildcard support")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::with_name("input-format")
                        .short("f")
                        .long("input-format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
                            "fasta",
                            "nexus",
                            "phylip",
                        ]),
                )
                .arg(
                    Arg::with_name("regex")
                        .long("re")
                        .help("Extract sequence IDs that match regular expression")
                        .conflicts_with_all(&["id", "file"])
                        .takes_value(true)
                        .require_equals(true)
                        .value_name("REGEX")       
                )
                .arg(
                    Arg::with_name("file")
                        .long("file")
                        .help("Inputs sequence IDs in a file")
                        .conflicts_with_all(&["id", "regex"])
                        .takes_value(true)
                        .value_name("A-TEXT-FILE-PATH")       
                )
                .arg(
                    Arg::with_name("id")
                        .long("id")
                        .help("Extract sequence using a list of IDs")
                        .conflicts_with_all(&["regex", "file"])
                        .required_unless_all(&["regex", "file"])
                        .takes_value(true)
                        .multiple(true)
                        .value_name("STRING")       
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Specifies a directory name")
                        .takes_value(true)
                        .required(true)
                        .default_value("SEGUL-Extract")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("output-format")
                        .short("F")
                        .long("output-format")
                        .help("Sets target output format")
                        .takes_value(true)
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
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
                    Arg::with_name("datatype")
                        .long("datatype")
                        .help("Sets data type")
                        .takes_value(true)
                        .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
        )
        .subcommand(
            App::new("translate")
                .about("Translate dna sequences to amino acid sequences")
                .arg(
                    Arg::with_name("dir")
                        .short("d")
                        .long("dir")
                        .help("Inputs a directory path to alignments")
                        .takes_value(true)
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .long("input")
                        .help("Input path with wildcard support")
                        .takes_value(true)
                        .multiple(true)
                        .required_unless_one(&["dir", "show-tables"])
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::with_name("input-format")
                        .short("f")
                        .long("input-format")
                        .help("Sets input format")
                        .takes_value(true)
                        .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .possible_values(&[
                            "auto",
                            "fasta",
                            "nexus",
                            "phylip",
                        ]),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .help("Specifies a directory name")
                        .takes_value(true)
                        .required(true)
                        .default_value("SEGUL-translation")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::with_name("output-format")
                        .short("F")
                        .long("output-format")
                        .help("Sets target output format")
                        .takes_value(true)
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
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
                    Arg::with_name("datatype")
                        .long("datatype")
                        .help("Sets data type")
                        .takes_value(true)
                        .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::with_name("table")
                        .long("table")
                        .help("Sets the ncbi translation table")
                        .takes_value(true)
                        .default_value("1")
                        .value_name("INTEGER")
                        .possible_values(&[
                            "1",
                            "2",
                            "3",
                            "4",
                            "5",
                            "6",
                            "9",
                            "10",
                            "11",
                            "12",
                            "13",
                            "14",
                            "16",
                            "21",
                            "22",
                            "23",
                            "24",
                            "25",
                            "26",
                            "29",
                            "30",
                            "33",
                        ]),
                )
                .arg(
                    Arg::with_name("reading-frame")
                        .long("rf")
                        .help("Sets the translation reading frame")
                        .takes_value(true)
                        .default_value("1")
                        .value_name("INTEGER")
                        .possible_values(&[
                            "1",
                            "2",
                            "3",
                            "auto",
                        ]),
                )
                .arg(
                    Arg::with_name("show-tables")
                        .long("show-tables")
                        .help("Show supported NCBI Genetic Code tables")
                        .takes_value(false)
                )
        )
        .get_matches()
}

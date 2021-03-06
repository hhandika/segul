use clap::{crate_description, crate_name, Arg, ArgMatches, Command};

#[cfg(not(tarpaulin_include))]
pub fn get_args(version: &str) -> ArgMatches {
    Command::new(crate_name!())
        .version(version)
        .about(crate_description!())
        .author("Heru Handika")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(
            Command::new("convert")
                .about("Convert sequence formats")
                .arg(
                    Arg::new("dir")
                        .short('d')
                        .long("dir")
                        .help("Input a directory path")
                        .takes_value(true)
                        .required_unless_present("input")
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        .takes_value(true)
                        .multiple_values(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Specify an output directory")
                        .takes_value(true)
                        // // .required(true)
                        .default_value("SEGUL-Convert")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .takes_value(true)
                        // // .required(true)
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
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
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
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        // // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("sort")
                        .long("sort")
                        .help("Sort the alignments")
                        .takes_value(false),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        .takes_value(false)
                ),
        )
        .subcommand(
            Command::new("concat")
                .about("Concatenate alignments")
                .arg(
                    Arg::new("dir")
                        .short('d')
                        .long("dir")
                        .help("Input alignment dir")
                        .takes_value(true)
                        .required_unless_present("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        .takes_value(true)
                        .multiple_values(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .takes_value(true)
                        // // .required(true)
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
                    Arg::new("partition")
                        .short('p')
                        .long("part")
                        .help("Specify partition format")
                        .takes_value(true)
                        // // .required(true)
                        .default_value("nexus")
                        .value_name("PART-FORMAT")
                        .possible_values(&["charset", "nexus", "raxml"]),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Use a custom output directory")
                        .takes_value(true)
                        // // .required(true)
                        .default_value("SEGUL-Concat")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
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
                    Arg::new("prefix")
                        .long("prefix")
                        .help("Use a custom output filename")
                        .takes_value(true)
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        // // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("codon")
                        .long("codon")
                        .help("Specify codon model partition format")
                        .takes_value(false),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        .takes_value(false)
                ),
        )
        .subcommand(
            Command::new("filter")
                .about("Filter alignments based on min taxon completeness, alignment length, and parsimony informative sites")
                .arg(
                    Arg::new("dir")
                        .short('d')
                        .long("dir")
                        .help("Input a path to alignment dir")
                        .takes_value(true)
                        .required_unless_present("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        .takes_value(true)
                        .multiple_values(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .takes_value(true)
                        // // .required(true)
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
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        // // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("percent")
                        .long("percent")
                        .help("Filtered based on percentage of minimal taxa")
                        .takes_value(true)
                        .required_unless_present_any(&["npercent", "aln-len", "pars-inf", "percent-inf", "taxon-all"])
                        .conflicts_with_all(&["npercent", "aln-len", "pars-inf", "percent-inf", "taxon-all"])
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("npercent")
                        .long("npercent")
                        .help("Filtered based on percentages of minimal taxa (allow multiple values)")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "aln-len", "pars-inf", "taxon-all"])
                        .multiple_values(true)
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("ntax")
                        .long("ntax")
                        .help("Input the total number of taxa")
                        .takes_value(true)
                        .conflicts_with_all(&["aln-len", "pars-inf", "taxon-all"])
                        .value_name("INTEGER"),
                )
                .arg(
                    Arg::new("aln-len")
                        .long("len")
                        .help("Filter based on minimal alignment length")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "npercent", "pars-inf", "taxon-all"])
                        .value_name("INTEGER"),
                )
                .arg(
                    Arg::new("pars-inf")
                        .long("pinf")
                        .help("Filter based on minimal parsimony informative sites")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "npercent", "aln-len", "taxon-all"])
                        .value_name("INTEGER"),
                )
                .arg(
                    Arg::new("percent-inf")
                        .long("percent-inf")
                        .help("Filter based on percent parsimony informative sites")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "npercent", "aln-len", "pars-inf", "taxon-all"])
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("taxon-all")
                        .long("taxon-all")
                        .help("Filter based on taxon ID")
                        .takes_value(true)
                        .conflicts_with_all(&["percent", "npercent", "aln-len", "pars-inf"])
                        .value_name("PATH-TO-taxon-all-FILE"),

                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Specify an output directory")
                        .takes_value(true)
                        .required_unless_present("dir")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format if concat")
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
                    Arg::new("prefix")
                        .long("prefix")
                        .help("Specify prefix for output filename")
                        .takes_value(true)
                        .requires("concat")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("partition")
                        .short('p')
                        .long("part")
                        .help("Specify a partition format if concat")
                        .takes_value(true)
                        .requires("concat")
                        .value_name("PART-FORMAT")
                        .possible_values(&["charset", "nexus", "raxml"]),
                )
                .arg(
                    Arg::new("concat")
                        .long("concat")
                        .help("Concat filtered alignments")
                        .requires("partition")
                        .takes_value(false),
                )
                .arg(
                    Arg::new("codon")
                        .long("codon")
                        .requires("concat")
                        .help("Specify a codon model partition format")
                        .takes_value(false),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        .takes_value(false)
                ),
        )
        .subcommand(
            Command::new("id")
                .about("Get sample ids from multiple alignments")
                .arg(
                    Arg::new("dir")
                        .short('d')
                        .long("dir")
                        .help("Input a directory path")
                        .takes_value(true)
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        .takes_value(true)
                        .multiple_values(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .takes_value(true)
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
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Use a custom output filename")
                        .takes_value(true)
                        .default_value("id")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("map")
                        .long("map")
                        .help("Map ID across all alignments")
                        .takes_value(false),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        .takes_value(false)
                ),
        )
        .subcommand(
            Command::new("remove")
            .about("Remove sequence based on user-defined IDs")
                .arg(
                    Arg::new("dir")
                    .short('d')
                    .long("dir")
                    .help("Input a directory path")
                    .takes_value(true)
                    .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                            .short('i')
                            .long("input")
                            .help("Input path (include wildcard support)")
                            .takes_value(true)
                            .multiple_values(true)
                            .required_unless_present("dir")
                            .conflicts_with("dir")
                            .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .takes_value(true)
                        // .required(true)
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
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Use a custom output filename")
                        .takes_value(true)
                        // .required(true)
                        .default_value("SEGUL-Remove")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
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
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("re")
                        .long("re")
                        .help("Remove sequence that match regular expression")
                        .conflicts_with("id")
                        .takes_value(true)
                        .require_equals(true)
                        .value_name("REGEX")       
                )
                .arg(
                    Arg::new("id")
                        .long("id")
                        .help("Input sequence IDs using terminal commands (STDIN)")
                        .required_unless_present("re")
                        .takes_value(true)
                        .multiple_values(true)
                        .value_name("STRING")       
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        .takes_value(false)
                ),
    )
        .subcommand(
            Command::new("rename")
            .about("Batch renaming sequence IDs in multiple_values alignments")
                .arg(
                    Arg::new("dir")
                    .short('d')
                    .long("dir")
                    .help("Input a directory path")
                    .takes_value(true)
                    .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                            .short('i')
                            .long("input")
                            .help("Input path (include wildcard support)")
                            .takes_value(true)
                            .multiple_values(true)
                            .required_unless_present("dir")
                            .conflicts_with("dir")
                            .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .takes_value(true)
                        // .required(true)
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
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Use a custom output filename")
                        .takes_value(true)
                        // .required(true)
                        .default_value("SEGUL-Rename")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
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
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("replace-id")
                        .long("replace")
                        .help("Rename using input IDs in a file")
                        .takes_value(true)
                        .required_unless_present_any(&["remove", "remove-re", "replace-from", "replace-from-re"])
                        .conflicts_with_all(&["remove", "remove-re", "remove-re-all", "replace-from", "replace-to"])
                        .value_name("PATH"),

                )
                .arg(
                    Arg::new("remove")
                        .long("remove")
                        .help("Remove matching input string")
                        .takes_value(true)
                        .conflicts_with_all(&["remove-re", "remove-re-all"])
                        .require_equals(true)
                        .value_name("STRING"),

                )
                .arg(
                    Arg::new("remove-re")
                        .long("remove-re")
                        .help("Remove first found matching regular expression")
                        .takes_value(true)
                        .conflicts_with_all(&["remove-re-all"])
                        .require_equals(true)
                        .value_name("REGEX"),

                )
                .arg(
                    Arg::new("remove-re-all")
                        .long("remove-re-all")
                        .help("Remove all matching regular expression")
                        .takes_value(true)
                        .conflicts_with_all(&["replace-id"])
                        .require_equals(true)
                        .value_name("REGEX"),

                )
                .arg(
                    Arg::new("replace-from")
                        .long("replace-from")
                        .help("Replace matching input string with the output string")
                        .takes_value(true)
                        .conflicts_with_all(&["remove", "remove-re", "remove-re-all"])
                        .require_equals(true)
                        .requires("replace-to")
                        .value_name("INPUT-STRING"),

                )
                .arg(
                    Arg::new("replace-from-re")
                        .long("replace-from-re")
                        .help("Replace matching input regex with the output string/regex")
                        .takes_value(true)
                        .conflicts_with_all(&["remove", "remove-re", "remove-re-all"])
                        .require_equals(true)
                        .requires("replace-to")
                        .value_name("INPUT-REGEX"),

                )
                .arg(
                    Arg::new("replace-to")
                        .long("replace-to")
                        .help("Replace matching input string with the output string")
                        .takes_value(true)
                        .conflicts_with_all(&["remove", "remove-re", "remove-re-all"])
                        .require_equals(true)
                        .value_name("OUTPUT-STRING"),

                )
                .arg(
                    Arg::new("dry-run")
                        .long("dry-run")
                        .help("Check if the program can parse the input ids correctly")
                        .takes_value(false),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        .takes_value(false)
                ),
        )
        .subcommand(
            Command::new("summary")
                .about("Compute alignment summary statistics")
                .arg(
                    Arg::new("dir")
                        .short('d')
                        .long("dir")
                        .help("Input a directory path")
                        .takes_value(true)
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        .takes_value(true)
                        .multiple_values(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .takes_value(true)
                        // .required(true)
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
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Use a custom output directory name")
                        .takes_value(true)
                        // .required(true)
                        .default_value("SEGUL-Summary")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("prefix")
                        .long("prefix")
                        .help("Use a costum output filename")
                        .takes_value(true)
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("percent-interval")
                        .long("interval")
                        .help("Specify a custom percentage interval value for counting data matrix completeness")
                        .takes_value(true)
                        .value_name("INTEGER")
                        .default_value("5")
                        .possible_values(&["1", "2", "5", "10"]),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        .takes_value(false)
                ),
        )
        .subcommand(
            Command::new("partition")
                .about("Convert partition formats")
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path")
                        .takes_value(true)
                        .multiple_values(true)
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("partition")
                        .short('p')
                        .long("input-part")
                        .help("Specify input partition format")
                        .takes_value(true)
                        .value_name("PART-FORMAT")
                        .possible_values(&["nexus", "raxml", "charset"]),
                )
                .arg(
                    Arg::new("output-partition")
                        .short('P')
                        .long("output-part")
                        .help("Specify output partition format")
                        .takes_value(true)
                        .value_name("PART-FORMAT")
                        // .required(true)
                        .default_value("nexus")
                        .possible_values(&["nexus", "raxml"]),
                )
                .arg(
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("codon")
                        .long("codon")
                        .help("Specify codon model partition format")
                        .takes_value(false),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        .takes_value(false)
                )
                .arg(
                    Arg::new("uncheck")
                        .long("uncheck")
                        .help("Skip checking partition formats")
                        .takes_value(false)
                ),
        )
        .subcommand(
            Command::new("extract")
                .about("Extract sequences from a collection of alignments")
                .arg(
                    Arg::new("dir")
                        .short('d')
                        .long("dir")
                        .help("Input a directory path")
                        .takes_value(true)
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        .takes_value(true)
                        .multiple_values(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .takes_value(true)
                        // .required(true)
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
                    Arg::new("re")
                        .long("re")
                        .help("Extract sequence IDs that match regular expression")
                        .conflicts_with_all(&["id", "file"])
                        .takes_value(true)
                        .require_equals(true)
                        .value_name("REGEX")       
                )
                .arg(
                    Arg::new("file")
                        .long("file")
                        .help("Input sequence IDs in a file")
                        .conflicts_with_all(&["id", "re"])
                        .takes_value(true)
                        .value_name("A-TEXT-FILE-PATH")       
                )
                .arg(
                    Arg::new("id")
                        .long("id")
                        .help("Input sequence IDs using terminal commands (STDIN)")
                        .conflicts_with_all(&["re", "file"])
                        .required_unless_present_any(&["re", "file"])
                        .takes_value(true)
                        .multiple_values(true)
                        .value_name("STRING")       
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Specify an output directory name")
                        .takes_value(true)
                        // .required(true)
                        .default_value("SEGUL-Extract")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
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
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        .takes_value(false)
                )
        )
        .subcommand(
            Command::new("translate")
                .about("Translate DNA sequences to amino acid sequences")
                .arg(
                    Arg::new("dir")
                        .short('d')
                        .long("dir")
                        .help("Input a directory path")
                        .takes_value(true)
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        .takes_value(true)
                        .multiple_values(true)
                        .required_unless_present_any(&["dir", "show-tables"])
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .takes_value(true)
                        // .required(true)
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
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Specify an output directory name")
                        .takes_value(true)
                        // .required(true)
                        .default_value("SEGUL-Translation")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
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
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("table")
                        .long("table")
                        .help("Specify the NCBI translation table")
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
                    Arg::new("reading-frame")
                        .long("rf")
                        .help("Specify the translation reading frame")
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
                    Arg::new("show-tables")
                        .long("show-tables")
                        .help("Show supported NCBI translation tables")
                        .takes_value(false)
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrites the output file/directory if exist")
                        .takes_value(false)
                )
        )
        .subcommand(
            Command::new("split")
                .about("Split alignments by partitions")
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input partition path")
                        .takes_value(true)
                        .multiple_values(false)
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Specify an output directory")
                        .takes_value(true)
                        .default_value("SEGUL-Split")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .takes_value(true)
                        // .required(true)
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
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
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
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        .takes_value(true)
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .possible_values(&["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("input-partition")
                        .short('I')
                        .long("input-part")
                        .help("Input sequence partitions")
                        .required(true)
                        .takes_value(true)
                        .value_name("PART-PATH"),
                )
                .arg(
                    Arg::new("partition")
                        .short('p')
                        .long("part")
                        .help("Specify partition format")
                        .takes_value(true)
                        .value_name("PART-FORMAT")
                        .possible_values(&["nexus", "raxml"]),
                )
                .arg(
                    Arg::new("prefix")
                        .long("prefix")
                        .help("Add prefix to output file names")
                        .takes_value(true)
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        .takes_value(false)
                )
                .arg(
                    Arg::new("uncheck")
                        .long("uncheck")
                        .help("Skip checking partition formats")
                        .takes_value(false)
                ),
        )
        .get_matches()
}

#[cfg(test)]
mod test {
    use std::path::Path;

    use assert_cmd::Command;
    use tempdir::TempDir;

    pub fn segul(dir: &Path) -> Command {
        let mut cmd = Command::cargo_bin("segul").unwrap();
        cmd.current_dir(dir);
        cmd
    }

    macro_rules! test_subcommand {
        ($func: ident, $subcommand: expr) => {
            #[test]
            fn $func() {
                let tmp_dir = TempDir::new("temp").unwrap();
                segul(tmp_dir.path())
                    .arg($subcommand)
                    .arg("-h")
                    .assert()
                    .success();
                tmp_dir.close().unwrap();
            }
        };
    }

    test_subcommand! {test_concat_subcommand, "concat"}
    test_subcommand! {test_convert_subcommand, "convert"}
    test_subcommand! {test_extract_subcommand, "extract"}
    test_subcommand! {test_filter_subcommand, "filter"}
    test_subcommand! {test_id_subcommand, "id"}
    test_subcommand! {test_partition_subcommand, "partition"}
    test_subcommand! {test_remove_subcommand, "remove"}
    test_subcommand! {test_rename_subcommand, "rename"}
    test_subcommand! {test_split_subcommand, "split"}
    test_subcommand! {test_summary_subcommand, "summary"}
    test_subcommand! {test_translate_subcommand, "translate"}
}

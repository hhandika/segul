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
                        
                        .required_unless_present("input")
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        .multiple(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Specify an output directory")
                        
                        // // .required(true)
                        .default_value("SEGUL-Convert")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        
                        // // .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .value_parser([
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
                        
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
                        .value_parser([
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
                        
                        // // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("sort")
                        .long("sort")
                        .help("Sort the alignments")
                        ,
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        
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
                        
                        .required_unless_present("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        
                        .multiple(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .value_parser([
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
                        
                        .default_value("nexus")
                        .value_name("PART-FORMAT")
                        .value_parser(["charset", "nexus", "raxml"]),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Use a custom output directory")
                        
                        .default_value("SEGUL-Concat")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
                        
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
                        .value_parser([
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
                        
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("codon")
                        .long("codon")
                        .help("Specify codon model partition format")
                        ,
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        
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
                        
                        .required_unless_present("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        
                        .multiple(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        
                        // // .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .value_parser([
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
                        
                        // // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("percent")
                        .long("percent")
                        .help("Filtered based on percentage of minimal taxa")
                        
                        .required_unless_present_any(["npercent", "aln-len", "pars-inf", "percent-inf", "taxon-all"])
                        .conflicts_with_all(["npercent", "aln-len", "pars-inf", "percent-inf", "taxon-all"])
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("npercent")
                        .long("npercent")
                        .help("Filtered based on percentages of minimal taxa (allow multiple values)")
                        
                        .conflicts_with_all(["percent", "aln-len", "pars-inf", "taxon-all"])
                        .multiple(true)
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("ntax")
                        .long("ntax")
                        .help("Input the total number of taxa")
                        
                        .conflicts_with_all(["aln-len", "pars-inf", "taxon-all"])
                        .value_name("INTEGER"),
                )
                .arg(
                    Arg::new("aln-len")
                        .long("len")
                        .help("Filter based on minimal alignment length")
                        
                        .conflicts_with_all(["percent", "npercent", "pars-inf", "taxon-all"])
                        .value_name("INTEGER"),
                )
                .arg(
                    Arg::new("pars-inf")
                        .long("pinf")
                        .help("Filter based on minimal parsimony informative sites")
                        
                        .conflicts_with_all(["percent", "npercent", "aln-len", "taxon-all"])
                        .value_name("INTEGER"),
                )
                .arg(
                    Arg::new("percent-inf")
                        .long("percent-inf")
                        .help("Filter based on percent parsimony informative sites")
                        
                        .conflicts_with_all(["percent", "npercent", "aln-len", "pars-inf", "taxon-all"])
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("taxon-all")
                        .long("taxon-all")
                        .help("Filter based on taxon ID")
                        
                        .conflicts_with_all(["percent", "npercent", "aln-len", "pars-inf"])
                        .value_name("PATH-TO-taxon-all-FILE"),

                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Specify an output directory")
                        
                        .required_unless_present("dir")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format if concat")
                        
                        .value_name("SEQ-FORMAT")
                        .value_parser([
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
                        
                        .requires("concat")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("partition")
                        .short('p')
                        .long("part")
                        .help("Specify a partition format if concat")
                        
                        .requires("concat")
                        .value_name("PART-FORMAT")
                        .value_parser(["charset", "nexus", "raxml"]),
                )
                .arg(
                    Arg::new("concat")
                        .long("concat")
                        .help("Concat filtered alignments")
                        .requires("partition")
                        ,
                )
                .arg(
                    Arg::new("codon")
                        .long("codon")
                        .requires("concat")
                        .help("Specify a codon model partition format")
                        ,
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        
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
                        
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        
                        .multiple(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .value_parser([
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
                        
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Use a custom output filename")
                        
                        .default_value("id")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("map")
                        .long("map")
                        .help("Map ID across all alignments")
                        ,
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        
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
                    
                    .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                            .short('i')
                            .long("input")
                            .help("Input path (include wildcard support)")
                            
                            .multiple(true)
                            .required_unless_present("dir")
                            .conflicts_with("dir")
                            .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        
                        // .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .value_parser([
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
                        
                        // .required(true)
                        .default_value("SEGUL-Remove")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
                        
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
                        .value_parser([
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
                        
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("re")
                        .long("re")
                        .help("Remove sequence that match regular expression")
                        .conflicts_with("id")
                        
                        .require_equals(true)
                        .value_name("REGEX")       
                )
                .arg(
                    Arg::new("id")
                        .long("id")
                        .help("Input sequence IDs using terminal commands (STDIN)")
                        .required_unless_present("re")
                        
                        .multiple(true)
                        .value_name("STRING")       
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        
                ),
    )
        .subcommand(
            Command::new("rename")
            .about("Batch renaming sequence IDs in multiple alignments")
                .arg(
                    Arg::new("dir")
                    .short('d')
                    .long("dir")
                    .help("Input a directory path")
                    
                    .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                            .short('i')
                            .long("input")
                            .help("Input path (include wildcard support)")
                            .multiple(true)
                            .required_unless_present("dir")
                            .conflicts_with("dir")
                            .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .value_parser([
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
                        
                        // .required(true)
                        .default_value("SEGUL-Rename")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
                        
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
                        .value_parser([
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
                        
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("replace-id")
                        .long("replace")
                        .help("Rename using input IDs in a file")
                        
                        .required_unless_present_any(["remove", "remove-re", "replace-from", "replace-from-re", "remove-re-all"])
                        .conflicts_with_all(["remove", "remove-re", "remove-re-all", "replace-from", "replace-to"])
                        .value_name("PATH"),

                )
                .arg(
                    Arg::new("remove")
                        .long("remove")
                        .help("Remove matching input string")
                        
                        .conflicts_with_all(["remove-re", "remove-re-all"])
                        .require_equals(true)
                        .value_name("STRING"),

                )
                .arg(
                    Arg::new("remove-re")
                        .long("remove-re")
                        .help("Remove first found matching regular expression")
                        
                        .conflicts_with_all(["remove-re-all"])
                        .require_equals(true)
                        .value_name("REGEX"),

                )
                .arg(
                    Arg::new("remove-re-all")
                        .long("remove-re-all")
                        .help("Remove all matching regular expression")
                        
                        .conflicts_with_all(["replace-id"])
                        .require_equals(true)
                        .value_name("REGEX"),

                )
                .arg(
                    Arg::new("replace-from")
                        .long("replace-from")
                        .help("Replace matching input string with the output string")
                        
                        .conflicts_with_all(["remove", "remove-re", "remove-re-all"])
                        .require_equals(true)
                        .requires("replace-to")
                        .value_name("INPUT-STRING"),

                )
                .arg(
                    Arg::new("replace-from-re")
                        .long("replace-from-re")
                        .help("Replace matching input regex with the output string/regex")
                        
                        .conflicts_with_all(["remove", "remove-re", "remove-re-all"])
                        .require_equals(true)
                        .requires("replace-to")
                        .value_name("INPUT-REGEX"),

                )
                .arg(
                    Arg::new("replace-to")
                        .long("replace-to")
                        .help("Replace matching input string with the output string")
                        
                        .conflicts_with_all(["remove", "remove-re", "remove-re-all"])
                        .require_equals(true)
                        .value_name("OUTPUT-STRING"),

                )
                .arg(
                    Arg::new("dry-run")
                        .long("dry-run")
                        .help("Check if the program can parse the input ids correctly")
                        ,
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        
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
                        
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        
                        .multiple(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        
                        // .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .value_parser([
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
                        
                        // .required(true)
                        .default_value("SEGUL-Summary")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("prefix")
                        .long("prefix")
                        .help("Use a costum output filename")
                        
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("percent-interval")
                        .long("interval")
                        .help("Specify a custom percentage interval value for counting data matrix completeness")
                        
                        .value_name("INTEGER")
                        .default_value("5")
                        .value_parser(["1", "2", "5", "10"]),
                )
                .arg(
                    Arg::new("per-locus")
                        .long("per-locus")
                        .help("Generate summary statistic for each locus")
                        .conflicts_with("percent-interval")
                        ,   
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        
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
                        
                        .multiple(true)
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("partition")
                        .short('p')
                        .long("input-part")
                        .help("Specify input partition format")
                        
                        .value_name("PART-FORMAT")
                        .value_parser(["nexus", "raxml", "charset"]),
                )
                .arg(
                    Arg::new("output-partition")
                        .short('P')
                        .long("output-part")
                        .help("Specify output partition format")
                        
                        .value_name("PART-FORMAT")
                        // .required(true)
                        .default_value("nexus")
                        .value_parser(["nexus", "raxml"]),
                )
                .arg(
                    Arg::new("datatype")
                        .long("datatype")
                        .help("Specify data type")
                        
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("codon")
                        .long("codon")
                        .help("Specify codon model partition format")
                        ,
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        
                )
                .arg(
                    Arg::new("uncheck")
                        .long("uncheck")
                        .help("Skip checking partition formats")
                        
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
                        
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        
                        .multiple(true)
                        .required_unless_present("dir")
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        
                        // .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .value_parser([
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
                        .conflicts_with_all(["id", "file"])
                        
                        .require_equals(true)
                        .value_name("REGEX")       
                )
                .arg(
                    Arg::new("file")
                        .long("file")
                        .help("Input sequence IDs in a file")
                        .conflicts_with_all(["id", "re"])
                        
                        .value_name("A-TEXT-FILE-PATH")       
                )
                .arg(
                    Arg::new("id")
                        .long("id")
                        .help("Input sequence IDs using terminal commands (STDIN)")
                        .conflicts_with_all(["re", "file"])
                        .required_unless_present_any(["re", "file"])
                        
                        .multiple(true)
                        .value_name("STRING")       
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Specify an output directory name")
                        
                        // .required(true)
                        .default_value("SEGUL-Extract")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
                        
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
                        .value_parser([
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
                        
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        
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
                        
                        .conflicts_with("input")
                        .value_name("PATH"),
                )
                .arg(
                    Arg::new("input")
                        .short('i')
                        .long("input")
                        .help("Input path (include wildcard support)")
                        
                        .multiple(true)
                        .required_unless_present_any(["dir", "show-tables"])
                        .conflicts_with("dir")
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        
                        // .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .value_parser([
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
                        
                        // .required(true)
                        .default_value("SEGUL-Translation")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("output-format")
                        .short('F')
                        .long("output-format")
                        .help("Specify an output sequence format")
                        
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
                        .value_parser([
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
                        
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("table")
                        .long("table")
                        .help("Specify the NCBI translation table")
                        
                        .default_value("1")
                        .value_name("INTEGER")
                        .value_parser([
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
                        
                        .default_value("1")
                        .value_name("INTEGER")
                        .value_parser([
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
                        
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrites the output file/directory if exist")
                        
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
                        
                        .multiple(false)
                        .value_name("INPUT-PATH"),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .help("Specify an output directory")
                        
                        .default_value("SEGUL-Split")
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("input-format")
                        .short('f')
                        .long("input-format")
                        .help("Specify an input format")
                        
                        // .required(true)
                        .value_name("SEQ-FORMAT")
                        .default_value("auto")
                        .value_parser([
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
                        
                        .default_value("nexus")
                        .value_name("SEQ-FORMAT")
                        .value_parser([
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
                        
                        // .required(true)
                        .value_name("DATATYPE")
                        .default_value("dna")
                        .value_parser(["dna", "aa", "ignore"]),
                )
                .arg(
                    Arg::new("input-partition")
                        .short('I')
                        .long("input-part")
                        .help("Input sequence partitions")
                        .required(true)
                        
                        .value_name("PART-PATH"),
                )
                .arg(
                    Arg::new("partition")
                        .short('p')
                        .long("part")
                        .help("Specify partition format")
                        
                        .value_name("PART-FORMAT")
                        .value_parser(["nexus", "raxml"]),
                )
                .arg(
                    Arg::new("prefix")
                        .long("prefix")
                        .help("Add prefix to output file names")
                        
                        .value_name("STRING"),
                )
                .arg(
                    Arg::new("overwrite")
                        .long("overwrite")
                        .help("Overwrite existing output file(s)/directory")
                        
                )
                .arg(
                    Arg::new("uncheck")
                        .long("uncheck")
                        .help("Skip checking partition formats")
                        
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

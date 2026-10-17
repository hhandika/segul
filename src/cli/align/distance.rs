use std::path::PathBuf;

use crate::core::align::distance::group::SampleGrouping;
use crate::core::align::distance::{AlignmentDistance, MatrixDiagonal, MatrixTriangle};
use crate::helper::finder::IDs;
use crate::helper::logger::AlignSeqLogger;
use crate::helper::types::{DataType, InputFmt};
use crate::stats::distance::DistanceModel;

use crate::cli::args::align::AlignDistanceArgs;
use crate::cli::{AlignSeqInput, InputCli, OutputCli, collect_paths};

impl InputCli for DistanceParser<'_> {}
impl OutputCli for DistanceParser<'_> {}
impl AlignSeqInput for DistanceParser<'_> {}

pub(in crate::cli) struct DistanceParser<'a> {
    args: &'a AlignDistanceArgs,
    input_dir: Option<PathBuf>,
}

impl<'a> DistanceParser<'a> {
    pub(in crate::cli) fn new(args: &'a AlignDistanceArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn compute(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let model = self.parse_model(&self.args.model);
        let triangle = self.parse_triangle(&self.args.triangle);
        let diagonal = self.parse_diagonal(&self.args.diagonal);
        if !model.is_valid_for(&datatype) {
            log::error!(
                "The {} model is not valid for {} data. \
                Use p, jc69, k2p, or tn93 for DNA, and p, jc69, or poisson for amino acid.",
                model,
                datatype
            );
            std::process::exit(1);
        }
        let task = "Pairwise Genetic Distance";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log(task);
        let grouping = self.parse_grouping(&files, &input_fmt, &datatype);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        let mut handle = AlignmentDistance::new(&input_fmt, &datatype, &model, &self.args.output)
            .with_triangle(triangle)
            .with_diagonal(diagonal)
            .with_long_format(self.args.long);
        if let Some(grouping) = grouping {
            handle = handle.with_grouping(grouping);
        }
        handle.compute(&files, self.args.prefix.as_deref());
    }

    fn parse_grouping(
        &self,
        files: &[PathBuf],
        input_fmt: &InputFmt,
        datatype: &DataType,
    ) -> Option<SampleGrouping> {
        let grouping = if let (Some(delimiter), Some(fields)) =
            (self.args.group_delimiter, self.args.group_fields)
        {
            SampleGrouping::from_delimiter(delimiter, fields)
        } else if let Some(pattern) = &self.args.group_regex {
            SampleGrouping::from_regex(pattern)
        } else if let Some(path) = &self.args.group_file {
            SampleGrouping::from_map_file(path)
        } else if self.args.group_auto {
            let ids: Vec<String> = IDs::new(files, input_fmt, datatype)
                .id_unique()
                .into_iter()
                .collect();
            match SampleGrouping::detect(&ids) {
                Some(suggestion) => {
                    log::info!(
                        "Detected groups: split on '{}', keep {} field(s), {} groups\n",
                        suggestion.delimiter,
                        suggestion.fields,
                        suggestion.groups
                    );
                    Ok(suggestion.to_grouping())
                }
                None => {
                    log::error!(
                        "Could not detect groups from sample IDs. \
                        Use --group-delimiter, --group-regex, or --group-file."
                    );
                    std::process::exit(1);
                }
            }
        } else {
            return None;
        };
        match grouping {
            Ok(grouping) => Some(grouping),
            Err(e) => {
                log::error!("{:#}", e);
                std::process::exit(1);
            }
        }
    }

    fn parse_model(&self, model: &str) -> DistanceModel {
        model
            .parse::<DistanceModel>()
            .expect("Failed to parse distance model")
    }

    fn parse_triangle(&self, triangle: &str) -> MatrixTriangle {
        triangle
            .parse::<MatrixTriangle>()
            .expect("Failed to parse matrix triangle")
    }

    fn parse_diagonal(&self, diagonal: &str) -> MatrixDiagonal {
        diagonal
            .parse::<MatrixDiagonal>()
            .expect("Failed to parse matrix diagonal")
    }
}

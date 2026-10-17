use std::path::PathBuf;

use crate::core::align::distance::{AlignmentDistance, MatrixTriangle};
use crate::helper::logger::AlignSeqLogger;
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
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        AlignmentDistance::new(&input_fmt, &datatype, &model, &self.args.output)
            .with_triangle(triangle)
            .compute(&files, self.args.prefix.as_deref());
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
}

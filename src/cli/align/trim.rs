use std::path::PathBuf;

use crate::{
    cli::{args::align::AlignTrimArgs, collect_paths, AlignSeqInput, InputCli, OutputCli},
    core::align::trim::{AlignmentTrimming, TrimmingParameters},
    helper::logger::AlignSeqLogger,
};

pub(in crate::cli) struct AlignTrimParser<'a> {
    pub args: &'a AlignTrimArgs,
    pub input_dir: Option<PathBuf>,
    pub missing_data: Option<f64>,
    pub pars_inf: Option<usize>,
    params: TrimmingParameters,
    param_counter: usize,
}

impl InputCli for AlignTrimParser<'_> {}
impl OutputCli for AlignTrimParser<'_> {}
impl AlignSeqInput for AlignTrimParser<'_> {}

impl<'a> AlignTrimParser<'a> {
    pub(in crate::cli) fn new(args: &'a AlignTrimArgs) -> Self {
        Self {
            args,
            input_dir: None,
            missing_data: None,
            pars_inf: None,
            param_counter: 0,
            params: TrimmingParameters::None,
        }
    }

    pub(in crate::cli) fn trim(&mut self) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        let output_fmt = self.parse_output_fmt(&self.args.out_fmt.output_fmt);
        let task = "Trim alignment";
        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);
        self.missing_data = self.args.missing;
        self.pars_inf = self.args.pars_inf;
        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log(task);
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        self.parse_params();
        self.check_multiple_params();
        let trim = AlignmentTrimming::new(
            &files,
            &input_fmt,
            &datatype,
            &self.args.output,
            &output_fmt,
            &self.params,
        );
        trim.trim();
    }

    fn parse_params(&mut self) {
        if let Some(missing_data) = self.missing_data {
            self.params = TrimmingParameters::MissingData(missing_data);
            self.param_counter += 1;
        }

        if let Some(pars_inf) = self.pars_inf {
            self.params = TrimmingParameters::ParsInf(pars_inf);
            self.param_counter += 1;
        }
    }

    fn check_multiple_params(&self) {
        if self.param_counter > 1 {
            panic!("Only one parameter is allowed");
        }
    }
}

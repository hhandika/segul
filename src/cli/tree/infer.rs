use std::path::PathBuf;

use segul_core::core::tree::estimate::{
    ParsimonyOptions, RootingMethod, TreeEstimateOptions, TreeEstimation,
};
use segul_core::helper::logger::AlignSeqLogger;
use segul_core::helper::types::DataType;
use segul_core::report::ReportContext;
use segul_core::writer::tree::TreeWriteOptions;

use crate::cli::args::tree::TreeInferArgs;
use crate::cli::{AlignSeqInput, InputCli, OutputCli, collect_paths};

pub(in crate::cli) struct TreeInferParser<'a> {
    args: &'a TreeInferArgs,
    input_dir: Option<PathBuf>,
}

impl InputCli for TreeInferParser<'_> {}
impl OutputCli for TreeInferParser<'_> {}
impl AlignSeqInput for TreeInferParser<'_> {}

impl<'a> TreeInferParser<'a> {
    pub(in crate::cli) fn new(args: &'a TreeInferArgs) -> Self {
        Self {
            args,
            input_dir: None,
        }
    }

    pub(in crate::cli) fn infer(&mut self, ctx: &ReportContext) {
        let input_fmt = self.parse_input_fmt(&self.args.in_fmt.input_fmt);
        let datatype = self.parse_datatype(&self.args.in_fmt.datatype);
        if matches!(datatype, DataType::Ignore) {
            log::error!(
                "Tree inference requires DNA or amino acid data. \
                Specify --datatype dna or --datatype aa."
            );
            std::process::exit(1);
        }
        if self.args.method.is_distance() && !self.args.model.is_valid_for(&datatype) {
            log::error!(
                "The {} model is not valid for {} data. \
                Use p, jc69, k2p, or tn93 for DNA, and p, jc69, or poisson for amino acid.",
                self.args.model,
                datatype
            );
            std::process::exit(1);
        }

        let dir = &self.args.io.dir;
        let files = collect_paths!(self, dir, input_fmt);

        AlignSeqLogger::new(
            self.input_dir.as_deref(),
            &input_fmt,
            &datatype,
            files.len(),
        )
        .log("Phylogenetic tree inference");

        self.check_output_dir_exist(&self.args.output, self.args.io.force);

        let rooting = if self.args.midpoint {
            RootingMethod::Midpoint
        } else if let Some(ref outgroup) = self.args.outgroup {
            RootingMethod::Outgroup(outgroup.clone())
        } else {
            RootingMethod::None
        };

        let options = TreeEstimateOptions {
            method: self.args.method,
            model: self.args.model,
            parsimony: ParsimonyOptions {
                search: self.args.search,
                swap: self.args.swap,
                replicates: self.args.replicates,
                gap_mode: self.args.gap_mode,
                branch_lengths: self.args.mp_branch_lengths,
                seed: self.args.seed,
            },
            bootstrap: self.args.bootstrap,
            seed: self.args.seed,
            rooting,
            output_fmt: self.args.output_fmt,
            single_file: self.args.single_file,
            write_opts: TreeWriteOptions {
                use_translate: !self.args.no_translate,
                gzip: self.args.gzip,
                ..Default::default()
            },
        };

        let estimation = TreeEstimation::new(&input_fmt, &datatype, options)
            .with_report_context(ctx.clone().with_input_dir(self.input_dir.clone()));

        if let Err(err) =
            estimation.estimate_with_report(&files, &self.args.output, self.args.prefix.as_deref())
        {
            log::error!("{:#}", err);
            std::process::exit(1);
        }
    }
}

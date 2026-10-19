use crate::{
    core::maf::{
        BlockFilter,
        convert::{MafConverter, MafGrouping},
    },
    helper::utils,
};

use super::MafInputCli;
use crate::cli::{InputCli, OutputCli, args::genomics::MafConvertArgs};

pub(in crate::cli) struct MafConvertParser<'a> {
    args: &'a MafConvertArgs,
}

impl InputCli for MafConvertParser<'_> {}
impl OutputCli for MafConvertParser<'_> {}
impl MafInputCli for MafConvertParser<'_> {}

impl<'a> MafConvertParser<'a> {
    pub(in crate::cli) fn new(args: &'a MafConvertArgs) -> Self {
        Self { args }
    }

    pub(in crate::cli) fn convert(&mut self) {
        let output_fmt = self.parse_output_fmt(&self.args.output_fmt);
        let files = self.collect_maf_files(&self.args.io);
        log::info!("{:18}: {}", "File counts", utils::fmt_num(&files.len()));
        log::info!("{:18}: {}", "Input format", "MAF");
        if let Some(bed) = &self.args.bed {
            log::info!("{:18}: {}", "BED file", bed.display());
        }
        log::info!("{:18}: {}\n", "Task", "MAF conversion");
        self.check_output_dir_exist(&self.args.output, self.args.io.force);
        let grouping: MafGrouping = self
            .args
            .group_by
            .parse()
            .expect("Group-by values are checked by the argument parser");
        let filter = BlockFilter {
            min_ref_len: self.args.min_block_len,
            min_species: self.args.min_species,
            min_score: self.args.min_score,
            species: self.parse_species(&self.args.common),
        };
        let mut converter = MafConverter::new(&files, &self.args.output, &output_fmt)
            .with_grouping(grouping)
            .with_clip(!self.args.no_clip)
            .with_strand_aware(self.args.strand_aware)
            .with_naming(self.parse_naming(&self.args.common))
            .with_ref_species(self.args.common.ref_species.clone())
            .with_filter(filter)
            .with_unmatched(self.args.write_unmatched)
            .with_prefix(self.args.prefix.clone());
        if let Some(bed) = &self.args.bed {
            converter = converter.with_bed(bed);
        }
        if let Err(e) = converter.convert() {
            log::error!("{:#}", e);
            std::process::exit(1);
        }
    }
}

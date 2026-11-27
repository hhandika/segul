use crate::cli::align::convert::ConvertParser;
use crate::cli::align::distance::DistanceParser;
use crate::cli::align::filter::FilterParser;
use crate::cli::align::partition::PartParser;
use crate::cli::align::split::SplitParser;
use crate::cli::align::summarize::SummaryParser;
use crate::cli::read::summarize::ReadSummaryCliParser;
use crate::cli::sequence::extract::ExtractParser;
use crate::cli::sequence::id::IdParser;
use crate::cli::sequence::remove::RemoveParser;
use crate::cli::sequence::rename::RenameParser;
use crate::cli::sequence::translate::TranslateParser;
use crate::cli::sequence::transpose::TransposeParser;

use segul_core::report::ReportContext;

use super::align::concat::ConcatParser;
use super::align::trim::AlignTrimParser;
use super::align::unalign::UnalignParser;
use super::args::MainSubcommand;
use super::args::align::{AlignmentSubcommand, PartitionSubcommand};
use super::args::genomics::{ContigSubcommand, MafSubcommand, SeqReadSubcommand};
use super::args::sequence::SequenceSubcommand;
use super::args::view::ViewSubcommand;
use super::contig::summarize::ContigCliParser;
use super::maf::convert::MafConvertParser;
use super::maf::summary::MafSummaryParser;
use super::report::{ReportCliParser, RunCliParser};
use super::sequence::addition::AdditionParser;
use super::sequence::filter::SequenceFilterParser;
use super::view::align::AlignViewParser;
use super::view::summary::SummaryViewParser;

pub(crate) fn match_cli_subcommand(subcommand: &MainSubcommand, ctx: &ReportContext) {
    match subcommand {
        MainSubcommand::RawRead(subcommand) => match_raw_read_subcommand(subcommand, ctx),
        MainSubcommand::Contig(subcommand) => match_contig_subcommand(subcommand, ctx),
        MainSubcommand::Maf(subcommand) => match_maf_subcommand(subcommand, ctx),
        MainSubcommand::Alignment(subcommand) => match_alignment_subcommand(subcommand, ctx),
        MainSubcommand::Partition(subcommand) => match_partition_subcommand(subcommand, ctx),
        MainSubcommand::Sequence(subcommand) => match_sequence_subcommand(subcommand, ctx),
        MainSubcommand::Run(args) => RunCliParser::new(args).run(ctx),
        MainSubcommand::Report(sub) => match sub {
            super::args::ReportSubcommand::Verify(args) => ReportCliParser::verify(args),
            super::args::ReportSubcommand::Methods(args) => ReportCliParser::methods(args),
            super::args::ReportSubcommand::Workflow(args) => ReportCliParser::workflow(args),
        },
        MainSubcommand::View(subcommand) => match_view_subcommand(subcommand),
    }
}

fn match_contig_subcommand(subcommand: &ContigSubcommand, ctx: &ReportContext) {
    match subcommand {
        ContigSubcommand::ContigSummary(args) => ContigCliParser::new(args).summarize(ctx),
    };
}

fn match_raw_read_subcommand(subcommand: &SeqReadSubcommand, ctx: &ReportContext) {
    match subcommand {
        SeqReadSubcommand::RawSummary(raw_args) => {
            ReadSummaryCliParser::new(raw_args).summarize(ctx)
        }
    };
}

fn match_maf_subcommand(subcommand: &MafSubcommand, ctx: &ReportContext) {
    match subcommand {
        MafSubcommand::Convert(args) => MafConvertParser::new(args).convert(ctx),
        MafSubcommand::Summary(args) => MafSummaryParser::new(args).summarize(ctx),
    };
}

fn match_partition_subcommand(subcommand: &PartitionSubcommand, ctx: &ReportContext) {
    match subcommand {
        PartitionSubcommand::Convert(part_args) => PartParser::new(part_args).convert(ctx),
    };
}

fn match_alignment_subcommand(subcommand: &AlignmentSubcommand, ctx: &ReportContext) {
    match subcommand {
        AlignmentSubcommand::Concat(concat_args) => ConcatParser::new(concat_args).concat(ctx),
        AlignmentSubcommand::Convert(convert_args) => ConvertParser::new(convert_args).convert(ctx),
        AlignmentSubcommand::Distance(distance_args) => {
            DistanceParser::new(distance_args).compute(ctx)
        }
        AlignmentSubcommand::Filter(filter_args) => FilterParser::new(filter_args).filter(ctx),
        AlignmentSubcommand::Split(split_args) => SplitParser::new(split_args).split(ctx),
        AlignmentSubcommand::Summary(summary_args) => {
            SummaryParser::new(summary_args).summarize(ctx)
        }
        AlignmentSubcommand::Trim(trim_args) => AlignTrimParser::new(trim_args).trim(ctx),
        AlignmentSubcommand::Unalign(unalign_args) => UnalignParser::new(unalign_args).unalign(ctx),
    };
}

fn match_sequence_subcommand(subcommand: &SequenceSubcommand, ctx: &ReportContext) {
    match subcommand {
        SequenceSubcommand::Extract(extract_args) => ExtractParser::new(extract_args).extract(ctx),
        SequenceSubcommand::Filter(filter_args) => {
            SequenceFilterParser::new(filter_args).filter(ctx)
        }
        SequenceSubcommand::Id(id_args) => IdParser::new(id_args).extract(ctx),
        SequenceSubcommand::Remove(remove_args) => RemoveParser::new(remove_args).remove(ctx),
        SequenceSubcommand::Rename(rename_args) => RenameParser::new(rename_args).rename(ctx),
        SequenceSubcommand::Translate(trans_args) => {
            TranslateParser::new(trans_args).translate(ctx)
        }
        SequenceSubcommand::Transpose(transpose_args) => {
            TransposeParser::new(transpose_args).transpose(ctx)
        }
        SequenceSubcommand::Add(add_args) => AdditionParser::new(add_args).add(ctx),
    };
}

pub(crate) fn match_view_subcommand(subcommand: &ViewSubcommand) {
    match subcommand {
        ViewSubcommand::Align(args) => AlignViewParser::new(args).view(),
        ViewSubcommand::Summary(args) => SummaryViewParser::new(args).view(),
    };
}

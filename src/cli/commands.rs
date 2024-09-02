use crate::cli::align::convert::ConvertParser;
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

use super::align::concat::ConcatParser;
use super::align::trim::AlignTrimParser;
use super::align::unalign::UnalignParser;
use super::args::align::{AlignmentSubcommand, PartitionSubcommand};
use super::args::genomics::{ContigSubcommand, MafSubcommand, SeqReadSubcommand};
use super::args::sequence::SequenceSubcommand;
use super::args::MainSubcommand;
use super::contig::summarize::ContigCliParser;
use super::maf::convert::MafConvertParser;
use super::sequence::addition::AdditionParser;
use super::sequence::filter::SequenceFilterParser;

pub(crate) fn match_cli_subcommand(subcommand: &MainSubcommand) {
    match subcommand {
        MainSubcommand::RawRead(subcommand) => match_raw_read_subcommand(subcommand),
        MainSubcommand::Contig(subcommand) => match_contig_subcommand(subcommand),
        MainSubcommand::Maf(subcommand) => match_maf_subcommand(subcommand),
        MainSubcommand::Alignment(subcommand) => match_alignment_subcommand(subcommand),
        MainSubcommand::Partition(subcommand) => match_partition_subcommand(subcommand),
        MainSubcommand::Sequence(subcommand) => match_sequence_subcommand(subcommand),
    }
}

fn match_contig_subcommand(subcommand: &ContigSubcommand) {
    match subcommand {
        ContigSubcommand::ContigSummary(args) => ContigCliParser::new(args).summarize(),
    };
}

fn match_raw_read_subcommand(subcommand: &SeqReadSubcommand) {
    match subcommand {
        SeqReadSubcommand::RawSummary(raw_args) => ReadSummaryCliParser::new(raw_args).summarize(),
    };
}

fn match_maf_subcommand(subcommand: &MafSubcommand) {
    match subcommand {
        MafSubcommand::MafConvert(maf_args) => {
            MafConvertParser::new(maf_args).convert();
        }
    };
}

fn match_partition_subcommand(subcommand: &PartitionSubcommand) {
    match subcommand {
        PartitionSubcommand::Convert(part_args) => PartParser::new(part_args).convert(),
    };
}

fn match_alignment_subcommand(subcommand: &AlignmentSubcommand) {
    match subcommand {
        AlignmentSubcommand::Concat(concat_args) => ConcatParser::new(concat_args).concat(),
        AlignmentSubcommand::Convert(convert_args) => ConvertParser::new(convert_args).convert(),
        AlignmentSubcommand::Filter(filter_args) => FilterParser::new(filter_args).filter(),
        AlignmentSubcommand::Split(split_args) => SplitParser::new(split_args).split(),
        AlignmentSubcommand::Summary(summary_args) => SummaryParser::new(summary_args).summarize(),
        AlignmentSubcommand::Trim(trim_args) => AlignTrimParser::new(trim_args).trim(),
        AlignmentSubcommand::Unalign(unalign_args) => UnalignParser::new(unalign_args).unalign(),
    };
}

fn match_sequence_subcommand(subcommand: &SequenceSubcommand) {
    match subcommand {
        SequenceSubcommand::Extract(extract_args) => ExtractParser::new(extract_args).extract(),
        SequenceSubcommand::Filter(filter_args) => SequenceFilterParser::new(filter_args).filter(),
        SequenceSubcommand::Id(id_args) => IdParser::new(id_args).extract(),
        SequenceSubcommand::Remove(remove_args) => RemoveParser::new(remove_args).remove(),
        SequenceSubcommand::Rename(rename_args) => RenameParser::new(rename_args).rename(),
        SequenceSubcommand::Translate(trans_args) => TranslateParser::new(trans_args).translate(),
        SequenceSubcommand::Add(add_args) => AdditionParser::new(add_args).add(),
    };
}

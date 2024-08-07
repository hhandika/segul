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
use super::args::align::{AlignmentSubcommand, PartitionSubcommand};
use super::args::genomics::{ContigSubcommand, SeqReadSubcommand};
use super::args::sequence::SequenceSubcommand;
use super::args::MainSubcommand;
use super::contig::summarize::ContigCliParser;
use super::sequence::filter::SequenceFilterParser;

pub(crate) fn match_cli_subcommand(subcommand: &MainSubcommand) {
    match subcommand {
        MainSubcommand::RawRead(subcommand) => match_raw_read_subcommand(subcommand),
        MainSubcommand::Contig(subcommand) => match_contig_subcommand(subcommand),
        MainSubcommand::Alignment(subcommand) => match_alignment_subcommand(subcommand),
        MainSubcommand::Partition(subcommand) => match_partition_subcommand(subcommand),
        MainSubcommand::Sequence(subcommand) => match_sequence_subcommand(subcommand),
    }
}

fn match_contig_subcommand(subcommand: &ContigSubcommand) {
    match subcommand {
        ContigSubcommand::ContigSummary(summary_args) => {
            ContigCliParser::new(summary_args).summarize();
        }
    };
}

fn match_raw_read_subcommand(subcommand: &SeqReadSubcommand) {
    match subcommand {
        SeqReadSubcommand::RawSummary(raw_args) => {
            ReadSummaryCliParser::new(raw_args).summarize();
        }
    };
}

fn match_partition_subcommand(subcommand: &PartitionSubcommand) {
    match subcommand {
        PartitionSubcommand::Convert(part_args) => {
            PartParser::new(part_args).convert();
        }
    };
}

fn match_alignment_subcommand(subcommand: &AlignmentSubcommand) {
    match subcommand {
        AlignmentSubcommand::Concat(concat_args) => {
            ConcatParser::new(concat_args).concat();
        }
        AlignmentSubcommand::Convert(convert_args) => {
            ConvertParser::new(convert_args).convert();
        }
        AlignmentSubcommand::Filter(filter_args) => {
            FilterParser::new(filter_args).filter();
        }
        AlignmentSubcommand::Split(split_args) => SplitParser::new(split_args).split(),
        AlignmentSubcommand::AlignSummary(summary_args) => {
            SummaryParser::new(summary_args).summarize();
        }
    };
}

fn match_sequence_subcommand(subcommand: &SequenceSubcommand) {
    match subcommand {
        SequenceSubcommand::Extract(extract_args) => {
            ExtractParser::new(extract_args).extract();
        }
        SequenceSubcommand::Filter(filter_args) => {
            SequenceFilterParser::new(filter_args).filter();
        }
        SequenceSubcommand::Id(id_args) => {
            IdParser::new(id_args).extract();
        }
        SequenceSubcommand::Remove(remove_args) => {
            RemoveParser::new(remove_args).remove();
        }
        SequenceSubcommand::Rename(rename_args) => {
            RenameParser::new(rename_args).rename();
        }
        SequenceSubcommand::Translate(translate_args) => {
            TranslateParser::new(translate_args).translate();
        }
    };
}

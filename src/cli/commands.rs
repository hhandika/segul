use crate::cli::args::MainSubcommand;
use crate::cli::concat::ConcatParser;
use crate::cli::convert::ConvertParser;
use crate::cli::extract::ExtractParser;
use crate::cli::filter::FilterParser;
use crate::cli::id::IdParser;
use crate::cli::partition::PartParser;
use crate::cli::read::ReadSummaryCliParser;
use crate::cli::remove::RemoveParser;
use crate::cli::rename::RenameParser;
use crate::cli::split::SplitParser;
use crate::cli::summarize::SummaryParser;
use crate::cli::translate::TranslateParser;

use super::{
    args::{
        AlignmentSubcommand, ContigSubcommand, PartitionSubcommand, SeqReadSubcommand,
        SequenceSubcommand,
    },
    contigs::ContigCliParser,
};

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
        SequenceSubcommand::Id(id_args) => {
            IdParser::new(id_args).find();
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

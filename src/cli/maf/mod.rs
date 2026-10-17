pub(crate) mod convert;
pub(crate) mod summary;

use std::path::{Path, PathBuf};

use crate::core::maf::{BlockFilter, SampleNaming};
use crate::helper::finder::MafFileFinder;

use crate::cli::InputCli;
use crate::cli::args::IOArgs;
use crate::cli::args::genomics::MafCommonArgs;

/// Shared input handling for the MAF subcommands.
trait MafInputCli: InputCli {
    fn collect_maf_files(&self, io: &IOArgs) -> Vec<PathBuf> {
        match &io.dir {
            Some(dir) => {
                log::info!("{:18}: {}", "Input dir", dir);
                MafFileFinder::new(Path::new(dir)).find()
            }
            None => {
                log::info!("{:18}: {}", "Input path", "STDIN");
                self.collect_paths(&io.input)
            }
        }
    }

    fn parse_naming(&self, common: &MafCommonArgs) -> SampleNaming {
        if common.full_source_name {
            SampleNaming::FullSource
        } else {
            SampleNaming::AssemblyPrefix
        }
    }

    /// Species from `--species` or `--species-file`. Exits on an unreadable file.
    fn parse_species(&self, common: &MafCommonArgs) -> Option<Vec<String>> {
        if let Some(list) = &common.species {
            return Some(BlockFilter::parse_species_list(list));
        }
        common.species_file.as_ref().map(|path| {
            BlockFilter::read_species_file(path).unwrap_or_else(|e| {
                log::error!("{:#}", e);
                std::process::exit(1);
            })
        })
    }
}

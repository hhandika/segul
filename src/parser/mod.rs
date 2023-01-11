//! All the parsers for the supported formats.
//!
//! Includes the following parsers:
//! 1. `delimited`: Delimited text format parser: csv and tsv.
//! 2. `fasta`: Fasta format parser.
//! 3. `nexus`: Nexus format parser.
//! 4. `partition`: Partition format parser.
//! 5. `phylip`: Phylip format parser.
//! 6. `txt`: Plain text format parser.
pub mod delimited;
pub mod fasta;
pub mod nexus;
pub mod partition;
pub mod phylip;
pub mod txt;

macro_rules! insert_matrix {
    ($self: ident, $id: ident, $seq: ident ) => {
        match $self.matrix.get($id) {
            Some(original_seq) => panic!(
                "Found duplicate IDs for file {}. \
                First duplicate ID found: {}. \
                Both sequences are the same: {}.",
                $self.input.display(),
                $id,
                original_seq == $seq
            ),
            None => {
                $self.matrix.insert($id.to_string(), $seq.to_string());
            }
        }
    };
}

macro_rules! warn_duplicate_ids {
    ($self: ident, $ids: ident) => {
        if $ids.len() != $self.header.ntax {
            log::warn!(
                "{} Found problematic alignment: {}.\n\
                Number of taxa in the matrix ({} taxa) \
                does not match the number of taxa in the header ({} taxa).\n\
                If this is `id` subcommand, try to use the same file as an input \
                for other subcommands.\n\
                If the cause is duplicate IDs, \
                the other subcommands will show the first found duplicate.\n\n",
                "WARNING!".red(),
                $self.input.display(),
                $ids.len(),
                $self.header.ntax
            );
        }
    };
}

pub(crate) use insert_matrix;
pub(crate) use warn_duplicate_ids;

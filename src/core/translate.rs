use std::fs;
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::helper::sequence::{SeqCheck, Sequence};
use crate::helper::types::{DataType, Header, InputFmt, OutputFmt, SeqMatrix};

pub enum TransTables {
    StandardCode,
    MtDna,
}

pub struct Translate<'a> {
    input_fmt: &'a InputFmt,
    trans_table: &'a TransTables,
    datatype: &'a DataType,
}

impl<'a> Translate<'a> {
    pub fn new(
        trans_table: &'a TransTables,
        input_fmt: &'a InputFmt,
        datatype: &'a DataType,
    ) -> Self {
        Self {
            trans_table,
            input_fmt,
            datatype,
        }
    }

    pub fn translate_sequences(&self, files: &[PathBuf], output: &Path, output_fmt: &OutputFmt) {
        fs::create_dir_all(output).expect("Failed creating an output directory");
        files.par_iter().for_each(|file| {
            let (mut seq, _) = Sequence::new(file, self.datatype).get(self.input_fmt);
            let _translate = self.translate(&mut seq);
        });

        match output_fmt {
            OutputFmt::Fasta => println!("Fasta"),
            _ => println!("Not fasta"),
        }
    }

    fn translate(&self, matrix: &mut SeqMatrix) {
        matrix.iter().for_each(|(id, seq)| {
            println!("{}: {}", id, seq);
        });

        matrix.clear();
    }
}

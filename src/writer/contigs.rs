//! Contig summary writer
use std::{
    io::{Result, Write},
    path::Path,
};

use crate::stats::contigs::ContigSummary;

use super::FileWriter;

const OUTPUT_FILENAME: &str = "contig-summary.csv";

pub struct ContigSummaryWriter<'a> {
    summary: &'a [ContigSummary],
    output: &'a Path,
}

impl FileWriter for ContigSummaryWriter<'_> {}

impl<'a> ContigSummaryWriter<'a> {
    pub fn new(summary: &'a [ContigSummary], output: &'a Path) -> Self {
        Self { summary, output }
    }

    pub fn write(&self) -> Result<()> {
        let output_path = self.output.join(OUTPUT_FILENAME);
        let mut writer = self
            .create_output_file(&output_path)
            .expect("Failed writing to file");
        self.write_records(&mut writer)?;
        writer.flush()?;
        Ok(())
    }

    fn write_records<W: Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(
            writer,
            "file_path,contig_name,\
            contig_count,base_count,nucleotide,\
            GC_content,AT_content,\
            sum,min,max,mean,median,\
            N50,N75,N90,\
            contig_750,contig_1000,contig_1500,\
            G_Count,C_Count,A_Count,T_Count"
        )?;
        for summary in self.summary {
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                summary.file_path,
                summary.contig_name,
                summary.contig_count,
                summary.base_count,
                summary.nucleotide,
                summary.gc_content,
                summary.at_content,
                summary.stats.sum,
                summary.stats.min,
                summary.stats.max,
                summary.stats.mean,
                summary.stats.median,
                summary.nstats.n50,
                summary.nstats.n75,
                summary.nstats.n90,
                summary.contig750,
                summary.contig1000,
                summary.contig1500,
                summary.g_count,
                summary.c_count,
                summary.a_count,
                summary.t_count,
            )?;
        }

        Ok(())
    }
}

//! Writer for read summary statistics.
use std::collections::BTreeMap;
use std::io::{Result, Write};
use std::path::{Path, PathBuf};

use super::FileWriter;

use crate::stats::fastq::{FastqSummary, FastqSummaryMin};
use crate::stats::qscores::ReadQScore;
use crate::stats::read::ReadRecord;

const DEFAULT_READ_SUFFIX: &str = "default-read-summary";
const MINIMAL_READ_SUFFIX: &str = "minimal-read-summary";
const DEFAULT_EXTENSION: &str = "csv";

pub struct ReadSummaryWriter<'a> {
    output: &'a Path,
    prefix: Option<&'a str>,
}

impl FileWriter for ReadSummaryWriter<'_> {}

impl<'a> ReadSummaryWriter<'a> {
    pub fn new(output: &'a Path, prefix: Option<&'a str>) -> Self {
        Self { output, prefix }
    }

    pub fn write(&self, records: &[FastqSummary]) -> Result<()> {
        let final_path = self.create_final_output_path(DEFAULT_READ_SUFFIX);
        let mut writer = self
            .create_output_file(&final_path)
            .expect("Failed writing to file");
        self.write_records(&mut writer, records)?;
        writer.flush()?;

        Ok(())
    }

    pub fn write_read_count_only(&self, records: &[FastqSummaryMin]) -> Result<()> {
        let output_path = self.create_final_output_path(MINIMAL_READ_SUFFIX);
        let mut writer = self
            .create_output_file(&output_path)
            .expect("Failed writing to file");
        writeln!(writer, "file_path,file_name,read_count")?;
        for rec in records {
            writeln!(
                writer,
                "{},{},{}",
                rec.path.display(),
                rec.file_name,
                rec.read_count
            )?;
        }
        writer.flush()?;
        Ok(())
    }

    fn create_final_output_path(&self, suffix: &str) -> PathBuf {
        match self.prefix {
            Some(prefix) => {
                let file_name = format!("{}_{}", prefix, suffix);
                self.output
                    .join(file_name)
                    .with_extension(DEFAULT_EXTENSION)
            }
            None => self.output.join(suffix).with_extension(DEFAULT_EXTENSION),
        }
    }

    /// Write the summary records to a file.
    fn write_records<W: Write>(&self, writer: &mut W, records: &[FastqSummary]) -> Result<()> {
        writeln!(
            writer,
            "file_path,file_name,read_count,base_count,\
                    mean_read_length,min_read_length,max_read_length,\
                    GC_count,GC_content,AT_count,AT_content,\
                    N_content,\
                    G,C,A,T,N,\
                    low_QScore,mean_QScore,min_QScore,max_QScore\
                    "
        )?;

        for rec in records {
            writeln!(
                writer,
                "{},{},{},{},\
                {},{},{},\
                {},{},{},{},\
                {},{},{},{},{},{},\
                {},{},{},{}",
                rec.path.display(),
                rec.file_name,
                rec.reads.stats.count,
                rec.reads.len,
                rec.reads.stats.mean,
                rec.reads.stats.min.unwrap_or(0),
                rec.reads.stats.max.unwrap_or(0),
                rec.read_summary.gc_count,
                rec.read_summary.gc_content,
                rec.read_summary.at_count,
                rec.read_summary.at_content,
                rec.read_summary.n_content,
                rec.reads.g_count,
                rec.reads.c_count,
                rec.reads.a_count,
                rec.reads.t_count,
                rec.reads.n_count,
                rec.qscores.low_q,
                rec.qscores.stats.mean,
                rec.qscores.stats.min.unwrap_or(0),
                rec.qscores.stats.max.unwrap_or(0),
            )?;
        }

        Ok(())
    }

    pub fn write_per_read_records(
        &self,
        fpath: &Path,
        reads: &BTreeMap<i32, ReadRecord>,
        qscores: &BTreeMap<i32, ReadQScore>,
    ) {
        let fname = format!(
            "{}_{}",
            fpath
                .file_stem()
                .expect("Failed getting file name")
                .to_str()
                .expect("Failed converting file name to string"),
            "read_summary.csv"
        );
        let output_dir = self.output.join("reads").join(fname);
        let mut writer = self
            .create_output_file(&output_dir)
            .expect("Failed writing to file");
        writeln!(
            writer,
            "index,G,C,A,T\
        ,proportion_G,proportion_C,proportion_A,proportion_T\
        ,mean_QScore,min_QScore,max_QScore",
        )
        .expect("Failed writing to file");
        reads.iter().for_each(|(i, r)| {
            let scores = if let Some(q) = qscores.get(i) {
                q
            } else {
                panic!("Failed getting quality scores for index {}", i);
            };
            let sum = r.g_count + r.c_count + r.a_count + r.t_count;
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{},{},{},{}",
                i,
                r.g_count,
                r.c_count,
                r.a_count,
                r.t_count,
                r.g_count as f64 / sum as f64,
                r.c_count as f64 / sum as f64,
                r.a_count as f64 / sum as f64,
                r.t_count as f64 / sum as f64,
                scores.stats.mean,
                scores.stats.min.unwrap_or(0),
                scores.stats.max.unwrap_or(0)
            )
            .expect("Failed writing to file");
        });
    }
}

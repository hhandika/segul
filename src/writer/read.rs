//! Writer for read summary statistics.
use std::collections::BTreeMap;
use std::io::{Result, Write};
use std::path::{Path, PathBuf};
use zip::write::FileOptions;

use super::FileWriter;

use crate::stats::fastq::{FastqMappedRead, FastqSummary, FastqSummaryMin};
use crate::stats::qscores::ReadQScore;
use crate::stats::read::ReadRecord;

const DEFAULT_READ_SUFFIX: &str = "default-read-summary";
const MINIMAL_READ_SUFFIX: &str = "minimal-read-summary";
const PER_READ_SUFFIX: &str = "read-pos-summary";

const DEFAULT_EXTENSION: &str = "csv";
const ZIP_EXTENSION: &str = "zip";

const READ_HEADER: &str = "file_path,file_name,read_count,base_count,\
mean_read_length,min_read_length,max_read_length,\
GC_count,GC_content,AT_count,AT_content,\
N_content,\
G,C,A,T,N,\
low_QScore,mean_QScore,min_QScore,max_QScore";

const PER_READ_HEADER: &str = "index,G,C,A,T\
,proportion_G,proportion_C,proportion_A,proportion_T\
,mean_QScore,min_QScore,max_QScore";

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
        writeln!(writer, "{}", READ_HEADER)?;

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
}

pub struct ReadPosSummaryWriter<'a> {
    output: &'a Path,
    prefix: Option<&'a str>,
}

impl FileWriter for ReadPosSummaryWriter<'_> {}

impl<'a> ReadPosSummaryWriter<'a> {
    pub fn new(output: &'a Path, prefix: Option<&'a str>) -> Self {
        Self { output, prefix }
    }

    pub fn write(&self, reads: &[FastqMappedRead]) -> Result<()> {
        let zip_path = self.create_final_output_path();
        let mut zip = self
            .create_zip_writer(&zip_path)
            .expect("Failed writing to file");
        let options = FileOptions::default().compression_method(zip::CompressionMethod::Deflated);
        reads.iter().for_each(|r| {
            let file_name = format!(
                "{}.{}",
                r.file_path
                    .file_stem()
                    .expect("Failed getting file name")
                    .to_str()
                    .expect("Failed converting file name to string"),
                DEFAULT_EXTENSION
            );

            zip.start_file(file_name, options)
                .expect("Failed writing to file");
            let content = self.parse_content(&r.reads, &r.qscores);
            zip.write(&content).expect("Failed writing to file");
        });

        zip.finish().expect("Failed writing to file");
        Ok(())
    }

    fn parse_content(
        &self,
        reads: &BTreeMap<i32, ReadRecord>,
        qscores: &BTreeMap<i32, ReadQScore>,
    ) -> Vec<u8> {
        let mut all_content: Vec<u8> = PER_READ_HEADER.as_bytes().to_vec();
        // Add new line after header
        all_content.push(b'\n');
        reads.iter().for_each(|(i, r)| {
            let scores = if let Some(q) = qscores.get(i) {
                q
            } else {
                panic!("Failed getting quality scores for index {}", i);
            };

            let content = self.format_content(i, r, scores);
            all_content.extend_from_slice(&content);
            // Add new line after each record
            all_content.push(b'\n');
        });
        all_content
    }

    fn format_content(&self, index: &i32, read: &ReadRecord, qscores: &ReadQScore) -> Vec<u8> {
        let sum = read.g_count + read.c_count + read.a_count + read.t_count;
        let data = format!(
            "{},{},{},{},{},{},{},{},{},{},{},{}",
            index,
            read.g_count,
            read.c_count,
            read.a_count,
            read.t_count,
            read.g_count as f64 / sum as f64,
            read.c_count as f64 / sum as f64,
            read.a_count as f64 / sum as f64,
            read.t_count as f64 / sum as f64,
            qscores.stats.mean,
            qscores.stats.min.unwrap_or(0),
            qscores.stats.max.unwrap_or(0)
        );
        data.as_bytes().to_vec()
    }

    fn create_final_output_path(&self) -> PathBuf {
        match self.prefix {
            Some(prefix) => {
                let file_name = format!("{}_{}", prefix, PER_READ_SUFFIX);
                self.output.join(file_name).with_extension(ZIP_EXTENSION)
            }
            None => self
                .output
                .join(PER_READ_SUFFIX)
                .with_extension(ZIP_EXTENSION),
        }
    }
}

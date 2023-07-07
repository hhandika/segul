use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Result, Write};
use std::path::Path;

use super::FileWriter;

use crate::stats::fastq::FastqSummary;
use crate::stats::read::{QScoreStream, ReadRecord};

const DEFAULT_OUTPUT: &str = "read-summary.csv";

pub struct ReadSummaryWriter<'a> {
    output: &'a Path,
}

impl FileWriter for ReadSummaryWriter<'_> {}

impl<'a> ReadSummaryWriter<'a> {
    pub fn new(output: &'a Path) -> Self {
        Self { output }
    }

    pub fn write(&self, records: &[FastqSummary]) -> Result<()> {
        let output_path = self.output.join(DEFAULT_OUTPUT);
        let mut writer = self
            .create_output_file(&output_path)
            .expect("Failed writing to file");
        self.write_records(&mut writer, records)?;
        writer.flush()?;

        Ok(())
    }

    /// Return a buffered writer for appending to the output file.
    pub fn write_append(&self) -> BufWriter<File> {
        let output_path = self.output.join(DEFAULT_OUTPUT);
        self.append_output_file(&output_path)
            .expect("Failed writing to file")
    }

    pub fn write_read_count_only(&self, records: &BTreeMap<String, usize>) -> Result<()> {
        let output_path = self.output.join(DEFAULT_OUTPUT);
        let mut writer = self
            .create_output_file(&output_path)
            .expect("Failed writing to file");
        writeln!(writer, "File,NumReads")?;
        for (path, count) in records {
            writeln!(writer, "{},{}", path, count)?;
        }
        writer.flush()?;
        Ok(())
    }

    /// Write the summary records to a file.
    pub fn write_records<W: Write>(&self, writer: &mut W, records: &[FastqSummary]) -> Result<()> {
        writeln!(
            writer,
            "File,NumReads,NumBases,\
                    MeanReadLen,MinReadLen,MaxReadLen,\
                    GCcount,GCcontent,ATcount,ATContent,\
                    Ncontent,\
                    G,C,A,T,N,\
                    LowQ,Mean,Min,Max\
                    "
        )?;

        for rec in records {
            writeln!(
                writer,
                "{},{},{},\
                {},{},{},\
                {},{},{},{},\
                {},{},{},{},{},{},\
                {},{},{},{}",
                rec.path.display(),
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
                rec.qscores.sum,
                rec.qscores.min,
                rec.qscores.max,
            )?;
        }

        Ok(())
    }

    pub fn write_per_read_records(
        &self,
        fpath: &Path,
        reads: &BTreeMap<i32, ReadRecord>,
        qscores: &BTreeMap<i32, QScoreStream>,
    ) {
        let fname = format!(
            "{}_{}",
            fpath
                .file_stem()
                .expect("Failed getting file name")
                .to_str()
                .expect("Failed converting file name to string"),
            "read_summary.tsv"
        );
        let output_dir = self.output.join("reads").join(fname);
        let mut writer = self
            .create_output_file(&output_dir)
            .expect("Failed writing to file");
        writeln!(
            writer,
            "index,G,C,A,T\
        ,ProportionG,ProportionC,ProportionA,ProportionT\
        ,MeanQ,MinQ,MaxQ",
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
                scores.mean,
                scores.min.unwrap_or(0),
                scores.max.unwrap_or(0)
            )
            .expect("Failed writing to file");
        });
    }
}

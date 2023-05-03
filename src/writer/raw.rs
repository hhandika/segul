use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Result, Write};
use std::path::Path;

use super::FileWriter;

use crate::stats::raw::{FastqRecords, QScoreRecords, QScoreStream, ReadRecord};

impl FileWriter for RawSummaryWriter<'_> {}

const DEFAULT_OUTPUT: &str = "summary.csv";

pub struct RawSummaryWriter<'a> {
    output: &'a Path,
}

impl<'a> RawSummaryWriter<'a> {
    pub fn new(output: &'a Path) -> Self {
        Self { output }
    }

    pub fn write(&self, records: &[(FastqRecords, QScoreRecords)]) -> Result<()> {
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
    pub fn write_records<W: Write>(
        &self,
        writer: &mut W,
        records: &[(FastqRecords, QScoreRecords)],
    ) -> Result<()> {
        writeln!(
            writer,
            "File,NumReads,NumBases,\
                    MinReadLen,MeanReadLen,MaxReadLen,\
                    GCcount,GCcontent,ATcount,ATContent,\
                    Ncount,Ncontent,\
                    LowQ,Mean,Min,Max\
                    "
        )?;
        for (seq, q) in records {
            writeln!(
                writer,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                seq.path.display(),
                seq.num_reads,
                seq.num_bases,
                seq.min_read_len,
                seq.mean_read_len,
                seq.max_read_len,
                seq.gc_count,
                seq.gc_content,
                seq.at_count,
                seq.at_content,
                seq.n_count,
                seq.n_content,
                q.low_q,
                q.mean,
                q.min,
                q.max,
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

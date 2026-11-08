//! Writer for read summary statistics.
use std::io::{Result, Write};
use std::path::{Path, PathBuf};

use zip::ZipWriter;
use zip::write::SimpleFileOptions;

use super::FileWriter;

use crate::stats::fastq::report::{ReadReport, SampleReadReport};
use crate::stats::fastq::{FastqSummary, FastqSummaryMin};

const DEFAULT_READ_SUFFIX: &str = "default-read-summary";
const MINIMAL_READ_SUFFIX: &str = "minimal-read-summary";
const PER_READ_SUFFIX: &str = "read-pos-summary";
const DISTRIBUTION_SUFFIX: &str = "read-dist-summary";
const REPORT_SUFFIX: &str = "read-report";

const DEFAULT_EXTENSION: &str = "csv";
const ZIP_EXTENSION: &str = "zip";
const JSON_EXTENSION: &str = "json";

// New columns are appended so readers that use column names keep working.
const READ_HEADER: &str = "file_path,file_name,read_count,base_count,\
mean_read_length,min_read_length,max_read_length,\
GC_count,GC_content,AT_count,AT_content,\
N_content,\
G,C,A,T,N,\
low_QScore,mean_QScore,min_QScore,max_QScore,\
median_read_length,stdev_read_length,N50_read_length,\
Q20_bases,Q20_pct,Q30_bases,Q30_pct,\
mean_read_QScore,low_quality_reads,phred_offset";

const COMPLETE_READ_HEADER: &str =
    ",duplicate_pct_est,adapter_reads_pct,polyG_reads_pct,overrepresented_count";

const PER_READ_HEADER: &str = "index,G,C,A,T\
,proportion_G,proportion_C,proportion_A,proportion_T\
,mean_QScore,min_QScore,max_QScore\
,N,proportion_N\
,P10_QScore,Q1_QScore,median_QScore,Q3_QScore,P90_QScore,end_index";

/// Output file path with an optional prefix.
fn output_path(output: &Path, prefix: Option<&str>, suffix: &str, extension: &str) -> PathBuf {
    let file_name = match prefix {
        Some(prefix) => format!("{}_{}", prefix, suffix),
        None => suffix.to_string(),
    };
    output.join(file_name).with_extension(extension)
}

pub struct ReadSummaryWriter<'a> {
    output: &'a Path,
    prefix: Option<&'a str>,
}

impl FileWriter for ReadSummaryWriter<'_> {}

impl<'a> ReadSummaryWriter<'a> {
    pub fn new(output: &'a Path, prefix: Option<&'a str>) -> Self {
        Self { output, prefix }
    }

    /// Write default statistics for each file.
    pub fn write(&self, records: &[FastqSummary]) -> Result<PathBuf> {
        self.write_summary(records, false)
    }

    /// Write default statistics, plus the Complete mode columns when `complete` is true.
    pub fn write_summary(&self, records: &[FastqSummary], complete: bool) -> Result<PathBuf> {
        let final_path = output_path(
            self.output,
            self.prefix,
            DEFAULT_READ_SUFFIX,
            DEFAULT_EXTENSION,
        );
        let mut writer = self
            .create_output_file(&final_path)
            .expect("Failed writing to file");
        self.write_records(&mut writer, records, complete)?;
        writer.flush()?;

        Ok(final_path)
    }

    pub fn write_read_count_only(&self, records: &[FastqSummaryMin]) -> Result<PathBuf> {
        let output_path = output_path(
            self.output,
            self.prefix,
            MINIMAL_READ_SUFFIX,
            DEFAULT_EXTENSION,
        );
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
        Ok(output_path)
    }

    /// Write the summary records to a file.
    fn write_records<W: Write>(
        &self,
        writer: &mut W,
        records: &[FastqSummary],
        complete: bool,
    ) -> Result<()> {
        if complete {
            writeln!(writer, "{}{}", READ_HEADER, COMPLETE_READ_HEADER)?;
        } else {
            writeln!(writer, "{}", READ_HEADER)?;
        }

        for rec in records {
            let m = rec.metrics();
            write!(
                writer,
                "{},{},{},{},\
                {},{},{},\
                {},{},{},{},\
                {},{},{},{},{},{},\
                {},{},{},{},\
                {},{},{},\
                {},{},{},{},\
                {},{},{}",
                rec.path.display(),
                rec.file_name,
                m.read_count,
                m.base_count,
                m.mean_read_length,
                m.min_read_length,
                m.max_read_length,
                rec.read_summary.gc_count,
                rec.read_summary.gc_content,
                rec.read_summary.at_count,
                rec.read_summary.at_content,
                rec.read_summary.n_content,
                m.g_count,
                m.c_count,
                m.a_count,
                m.t_count,
                m.n_count,
                m.low_qscore_bases,
                m.mean_qscore,
                m.min_qscore,
                m.max_qscore,
                m.median_read_length,
                m.stdev_read_length,
                m.n50_read_length,
                m.q20_bases,
                m.q20_pct,
                m.q30_bases,
                m.q30_pct,
                m.mean_read_qscore,
                m.low_quality_reads,
                m.phred_offset,
            )?;
            if complete {
                write!(
                    writer,
                    ",{},{},{},{}",
                    m.duplicate_pct_est.unwrap_or(0.0),
                    m.adapter_reads_pct.unwrap_or(0.0),
                    m.poly_g_reads_pct.unwrap_or(0.0),
                    m.overrepresented_count.unwrap_or(0),
                )?;
            }
            writeln!(writer)?;
        }

        Ok(())
    }
}

/// Writes per-position statistics, one CSV per sample, to `read-pos-summary.zip`.
pub struct ReadPosSummaryWriter<'a> {
    output: &'a Path,
    prefix: Option<&'a str>,
}

impl FileWriter for ReadPosSummaryWriter<'_> {}

impl<'a> ReadPosSummaryWriter<'a> {
    pub fn new(output: &'a Path, prefix: Option<&'a str>) -> Self {
        Self { output, prefix }
    }

    pub fn write(&self, reports: &[SampleReadReport]) -> Result<PathBuf> {
        let zip_path = output_path(self.output, self.prefix, PER_READ_SUFFIX, ZIP_EXTENSION);
        let mut zip = self
            .create_zip_writer(&zip_path)
            .expect("Failed writing to file");
        for report in reports {
            start_entry(
                &mut zip,
                &format!("{}.{}", report.sample, DEFAULT_EXTENSION),
            )?;
            self.write_content(&mut zip, report)?;
        }
        zip.finish()?;
        Ok(zip_path)
    }

    fn write_content<W: Write>(&self, zip: &mut W, report: &SampleReadReport) -> Result<()> {
        writeln!(zip, "{}", PER_READ_HEADER)?;
        for p in &report.positions {
            // Proportions of G, C, A, and T leave out ambiguous bases, as before.
            let acgt = (p.g_count + p.c_count + p.a_count + p.t_count) as f64;
            writeln!(
                zip,
                "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                p.start,
                p.g_count,
                p.c_count,
                p.a_count,
                p.t_count,
                p.g_count as f64 / acgt,
                p.c_count as f64 / acgt,
                p.a_count as f64 / acgt,
                p.t_count as f64 / acgt,
                p.mean_qscore,
                p.min_qscore,
                p.max_qscore,
                p.n_count,
                p.n_count as f64 / p.base_count() as f64,
                p.p10_qscore,
                p.q1_qscore,
                p.median_qscore,
                p.q3_qscore,
                p.p90_qscore,
                p.end,
            )?;
        }
        Ok(())
    }
}

/// Writes read distributions to `read-dist-summary.zip`, in one folder per sample.
pub struct ReadDistSummaryWriter<'a> {
    output: &'a Path,
    prefix: Option<&'a str>,
}

impl FileWriter for ReadDistSummaryWriter<'_> {}

impl<'a> ReadDistSummaryWriter<'a> {
    pub fn new(output: &'a Path, prefix: Option<&'a str>) -> Self {
        Self { output, prefix }
    }

    pub fn write(&self, reports: &[SampleReadReport]) -> Result<PathBuf> {
        let zip_path = output_path(self.output, self.prefix, DISTRIBUTION_SUFFIX, ZIP_EXTENSION);
        let mut zip = self
            .create_zip_writer(&zip_path)
            .expect("Failed writing to file");
        for report in reports {
            self.write_sample(&mut zip, report)?;
        }
        zip.finish()?;
        Ok(zip_path)
    }

    fn write_sample<W: Write + std::io::Seek>(
        &self,
        zip: &mut ZipWriter<W>,
        report: &SampleReadReport,
    ) -> Result<()> {
        let entry = |name: &str| format!("{}/{}.{}", report.sample, name, DEFAULT_EXTENSION);

        start_entry(zip, &entry("length"))?;
        writeln!(zip, "length,reads")?;
        for l in &report.length_distribution {
            writeln!(zip, "{},{}", l.length, l.reads)?;
        }

        start_entry(zip, &entry("read-qscore"))?;
        writeln!(zip, "mean_QScore,reads")?;
        for (q, reads) in report.read_quality_distribution.iter().enumerate() {
            writeln!(zip, "{},{}", q, reads)?;
        }

        start_entry(zip, &entry("gc"))?;
        writeln!(zip, "GC_percent,reads,theoretical_reads")?;
        for (gc, (reads, theoretical)) in report
            .gc_distribution
            .iter()
            .zip(&report.gc_theoretical)
            .enumerate()
        {
            writeln!(zip, "{},{},{}", gc, reads, theoretical)?;
        }

        if let Some(duplication) = &report.duplication {
            start_entry(zip, &entry("duplication"))?;
            writeln!(zip, "level,distinct_pct,reads_pct")?;
            for level in &duplication.levels {
                writeln!(
                    zip,
                    "{},{},{}",
                    level.label, level.distinct_pct, level.reads_pct
                )?;
            }
        }

        start_entry(zip, &entry("overrepresented"))?;
        writeln!(zip, "sequence,count,percent,source")?;
        for seq in &report.overrepresented {
            writeln!(
                zip,
                "{},{},{},{}",
                seq.sequence, seq.count, seq.percent, seq.source
            )?;
        }

        if !report.adapters.is_empty() {
            start_entry(zip, &entry("adapters"))?;
            let names: Vec<&str> = report.adapters.iter().map(|a| a.name.as_str()).collect();
            writeln!(zip, "index,end_index,{}", names.join(","))?;
            for (i, p) in report.positions.iter().enumerate() {
                let values: Vec<String> = report
                    .adapters
                    .iter()
                    .map(|a| a.cumulative_pct[i].to_string())
                    .collect();
                writeln!(zip, "{},{},{}", p.start, p.end, values.join(","))?;
            }
        }
        Ok(())
    }
}

/// Writes all statistics to `read-report.json`.
pub struct ReadReportWriter<'a> {
    output: &'a Path,
    prefix: Option<&'a str>,
}

impl FileWriter for ReadReportWriter<'_> {}

impl<'a> ReadReportWriter<'a> {
    pub fn new(output: &'a Path, prefix: Option<&'a str>) -> Self {
        Self { output, prefix }
    }

    pub fn write(&self, report: &ReadReport) -> anyhow::Result<PathBuf> {
        let path = output_path(self.output, self.prefix, REPORT_SUFFIX, JSON_EXTENSION);
        let mut writer = self.create_output_file(&path)?;
        serde_json::to_writer(&mut writer, report)?;
        writer.flush()?;
        Ok(path)
    }
}

fn start_entry<W: Write + std::io::Seek>(zip: &mut ZipWriter<W>, name: &str) -> Result<()> {
    let options = SimpleFileOptions::default().compression_method(zip::CompressionMethod::Deflated);
    zip.start_file(name, options)?;
    Ok(())
}

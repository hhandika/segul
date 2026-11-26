//! Assembly decompression, report parsing, and ENA contig filtering.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::Result;
use flate2::read::GzDecoder;

use super::client::CancelSignal;
use super::fasta::open_maybe_gzip;
use super::types::{DownloadError, base_accession};

const BUF_SIZE: usize = 256 * 1024;

/// Decompresses a gzipped FASTA `source` file into `dest`.
pub fn decompress_fasta(source: &Path, dest: &Path, cancel: &CancelSignal) -> Result<()> {
    let file = File::open(source)?;
    let mut decoder = GzDecoder::new(BufReader::with_capacity(BUF_SIZE, file));
    let mut writer = BufWriter::with_capacity(BUF_SIZE, File::create(dest)?);
    let mut buf = vec![0u8; BUF_SIZE];
    loop {
        cancel.throw_if_cancelled()?;
        match decoder.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => writer.write_all(&buf[..n])?,
            Err(err)
                if matches!(
                    err.kind(),
                    std::io::ErrorKind::InvalidData
                        | std::io::ErrorKind::InvalidInput
                        | std::io::ErrorKind::UnexpectedEof
                        | std::io::ErrorKind::Other
                ) =>
            {
                return Err(
                    DownloadError("The downloaded assembly is damaged.".to_string()).into(),
                );
            }
            Err(err) => return Err(err.into()),
        }
    }
    writer.flush()?;
    Ok(())
}

/// Parses accessions from an ENA sequence report, without version suffixes.
pub fn parse_sequence_report(report: &str) -> Result<HashSet<String>> {
    let mut lines = report.lines();
    let Some(first) = lines.next() else {
        return Err(DownloadError("ENA sent an unreadable assembly report.".to_string()).into());
    };
    if !first.starts_with("accession\t") {
        return Err(DownloadError("ENA sent an unreadable assembly report.".to_string()).into());
    }
    let mut set = HashSet::new();
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(acc) = trimmed.split('\t').next() {
            set.insert(base_accession(acc));
        }
    }
    Ok(set)
}

/// Filters an ENA assembly FASTA (plain or gzipped) to `members` and checks `expected_bases`.
pub fn write_ena_assembly(
    source: &Path,
    dest: &Path,
    members: Option<&HashSet<String>>,
    expected_bases: u64,
    cancel: &CancelSignal,
) -> Result<()> {
    let mut reader = open_maybe_gzip(source)?;
    let mut writer = BufWriter::with_capacity(BUF_SIZE, File::create(dest)?);
    let mut seen = HashSet::new();
    let mut bases: u64 = 0;
    let mut is_kept = false;
    let mut line_buf = Vec::with_capacity(4096);

    loop {
        line_buf.clear();
        let read_bytes = reader.read_until(b'\n', &mut line_buf)?;
        if read_bytes == 0 {
            break;
        }
        cancel.throw_if_cancelled()?;
        let line_lossy = String::from_utf8_lossy(&line_buf);
        let line = line_lossy.trim_end_matches(['\r', '\n']);

        if let Some(header) = line.strip_prefix('>') {
            let first = header.split_whitespace().next().unwrap_or("");
            let id = base_accession(first.split('|').next_back().unwrap_or(first));
            let in_members = members.is_none_or(|m| m.contains(&id));
            is_kept = in_members && seen.insert(id);
        } else if is_kept {
            bases = bases.saturating_add(line.trim().len() as u64);
        }

        if is_kept {
            writeln!(writer, "{line}")?;
        }
    }
    writer.flush()?;

    let has_all = members.is_none_or(|m| m.iter().all(|id| seen.contains(id)));
    if !has_all || bases != expected_bases {
        return Err(DownloadError(
            "ENA sent an incomplete or different assembly. Try GenBank.".to_string(),
        )
        .into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn test_write_ena_assembly_filters_gzipped_fasta() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("ena.fa.gz");
        let dst = dir.path().join("out.fna");
        let fasta = ">ENA|CM000001|CM000001.1 chr\nACGTACGT\n\
                     >ENA|JACO02000001|JACO02000001.1 contig\nACGTACGT\n\
                     >ENA|JACO02000002|JACO02000002.1 unplaced\nGGCC\n";
        {
            let mut enc = GzEncoder::new(File::create(&src).unwrap(), Compression::default());
            enc.write_all(fasta.as_bytes()).unwrap();
            enc.finish().unwrap();
        }
        let report = "accession\tsequence-name\nCM000001.1\tI\nJACO02000002.1\tunplaced\n";
        let members = parse_sequence_report(report).unwrap();
        write_ena_assembly(&src, &dst, Some(&members), 12, &CancelSignal::new()).unwrap();
        let out = fs::read_to_string(&dst).unwrap();
        assert!(!out.contains("JACO02000001"));
        assert!(out.contains("JACO02000002"));
    }

    #[test]
    fn test_decompress_fasta_rejects_corrupt() {
        let dir = TempDir::new().unwrap();
        let src = dir.path().join("bad.gz");
        let dst = dir.path().join("out.fna");
        fs::write(&src, b"not gzip").unwrap();
        let err = decompress_fasta(&src, &dst, &CancelSignal::new()).unwrap_err();
        assert!(err.to_string().contains("damaged"));
    }
}

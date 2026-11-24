//! Streams archive FASTA into SEGUL-safe FASTA with unique headers.

use std::collections::{HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::LazyLock;

use anyhow::Result;
use flate2::read::GzDecoder;
use regex::Regex;

use super::client::CancelSignal;
use super::types::{
    DatabaseSource, DownloadContent, DownloadError, DownloadOptions, DownloadRecord,
    accession_matches, base_accession, safe_name,
};

static GENE_KEY_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[\s_\-]+").expect("Valid gene key regex"));
static NOT_SEQUENCE_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[^A-Za-z*.~?\-]").expect("Valid sequence regex"));
static NOT_DNA_RE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"[^?\-ACGTRYSWKMBDHVNacgtryswkmbdhvn.]").expect("Valid DNA regex")
});
static NOT_AA_RE: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"[^?\-ARNDCQEGHILKMFPSTWYVXBZJU*.~]").expect("Valid AA regex"));

const IO_BUF_SIZE: usize = 256 * 1024;

/// Parsed origin metadata from an archive FASTA header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaOrigin {
    pub accession: String,
    pub gene: String,
    pub product: String,
    pub codon_start: Option<u32>,
}

/// A sequence written by [`prepare_batch`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreparedSequence {
    pub record: DownloadRecord,
    pub id: String,
    pub header: String,
    pub trimmed: u32,
}

/// Result of preparing a single downloaded FASTA batch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreparedBatch {
    pub sequences: Vec<PreparedSequence>,
    pub has_unknown_residues: bool,
}

/// Extracts accession and coding metadata from a GenBank or ENA FASTA header.
pub fn read_header(source: DatabaseSource, header: &str) -> FastaOrigin {
    match source {
        DatabaseSource::GenBank => read_genbank_header(header),
        DatabaseSource::Ena => read_ena_header(header),
    }
}

/// Opens `path` for buffered reading, decoding gzip transparently if magic bytes match.
pub fn open_maybe_gzip(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;
    file.seek(SeekFrom::Start(0))?;
    if n == 2 && magic == [0x1f, 0x8b] {
        let decoder = GzDecoder::new(BufReader::with_capacity(IO_BUF_SIZE, file));
        Ok(Box::new(BufReader::with_capacity(IO_BUF_SIZE, decoder)))
    } else {
        Ok(Box::new(BufReader::with_capacity(IO_BUF_SIZE, file)))
    }
}

/// Streams archive FASTA from `raw_path`, filters/trims records, and writes to `part_paths`.
pub fn prepare_batch(
    raw_path: &Path,
    source: DatabaseSource,
    records: &[DownloadRecord],
    options: &DownloadOptions,
    trims_to_frame: bool,
    part_paths: &HashMap<String, PathBuf>,
    cancel: &CancelSignal,
) -> Result<PreparedBatch> {
    let mut used_ids = HashSet::new();
    let wanted: HashMap<String, &DownloadRecord> = records
        .iter()
        .map(|r| (base_accession(&r.accession), r))
        .collect();

    let mut reader = open_maybe_gzip(raw_path)?;
    let mut sequences = Vec::new();
    let mut has_unknown_residues = false;
    let mut has_header = false;
    let mut sink: Option<BufWriter<File>> = None;
    let mut sink_path: Option<PathBuf> = None;

    struct PendingSeq<'a> {
        record: &'a DownloadRecord,
        base_id: String,
        header: String,
        trimmed: u32,
    }

    let mut current: Option<PendingSeq<'_>> = None;
    let mut is_written = false;
    let mut to_trim: u32 = 0;
    let mut line_buf = Vec::with_capacity(4096);

    loop {
        line_buf.clear();
        let read_bytes = reader.read_until(b'\n', &mut line_buf)?;
        if read_bytes == 0 {
            break;
        }
        cancel.throw_if_cancelled()?;
        let line_lossy = String::from_utf8_lossy(&line_buf);
        let line = line_lossy.trim();
        if line.is_empty() {
            continue;
        }

        if let Some(rest) = line.strip_prefix('>') {
            has_header = true;
            current = None;
            is_written = false;
            let header = rest.trim().to_string();
            let origin = read_header(source, &header);
            let base = base_accession(&origin.accession);
            let Some(&record) = wanted.get(&base) else {
                continue;
            };
            if !accession_matches(&record.accession, &origin.accession)
                || !matches_gene(record, &origin, &options.gene)
            {
                continue;
            }
            to_trim = if trims_to_frame && record.content == DownloadContent::Coding {
                origin
                    .codon_start
                    .unwrap_or(record.codon_start)
                    .saturating_sub(1)
                    .min(2)
            } else {
                0
            };
            current = Some(PendingSeq {
                record,
                base_id: base_id(record, &origin, options),
                header,
                trimmed: to_trim,
            });
            continue;
        }

        if !has_header || NOT_SEQUENCE_RE.is_match(line) {
            return Err(DownloadError(format!(
                "{} sent an error instead of sequences. Try again.",
                source.label()
            ))
            .into());
        }

        let Some(seq) = current.as_ref() else {
            continue;
        };

        let residues: &str = if to_trim > 0 {
            let cut = (to_trim as usize).min(line.len());
            to_trim -= cut as u32;
            let rem = &line[cut..];
            if rem.is_empty() {
                continue;
            }
            rem
        } else {
            line
        };

        if !is_written {
            is_written = true;
            let key = seq.record.key();
            let output = part_paths
                .get(&key)
                .ok_or_else(|| io::Error::new(io::ErrorKind::NotFound, "Missing part path"))?;
            if sink_path.as_ref() != Some(output) {
                if let Some(mut prev) = sink.take() {
                    prev.flush()?;
                }
                if let Some(parent) = output.parent() {
                    fs::create_dir_all(parent)?;
                }
                let file = OpenOptions::new().create(true).append(true).open(output)?;
                sink = Some(BufWriter::with_capacity(IO_BUF_SIZE, file));
                sink_path = Some(output.clone());
            }
            let id = unique_id(&seq.base_id, &mut used_ids);
            sequences.push(PreparedSequence {
                record: seq.record.clone(),
                id: id.clone(),
                header: seq.header.clone(),
                trimmed: seq.trimmed,
            });
            if let Some(writer) = sink.as_mut() {
                writeln!(writer, ">{id}")?;
            }
        }

        let invalid = if seq.record.content == DownloadContent::Protein {
            &*NOT_AA_RE
        } else {
            &*NOT_DNA_RE
        };
        if !has_unknown_residues && invalid.is_match(residues) {
            has_unknown_residues = true;
        }
        if let Some(writer) = sink.as_mut() {
            writeln!(writer, "{residues}")?;
        }
    }

    if let Some(mut writer) = sink {
        writer.flush()?;
    }

    Ok(PreparedBatch {
        sequences,
        has_unknown_residues,
    })
}

fn read_genbank_header(header: &str) -> FastaOrigin {
    let mut id = header.split_whitespace().next().unwrap_or("");
    let owned_id;
    if let Some(rest) = id.strip_prefix("lcl|") {
        if let Some(pos) = rest.find("_cds") {
            id = if pos > 0 { &rest[..pos] } else { rest };
        } else {
            id = rest;
        }
    } else if id.contains('|') {
        let parts: Vec<&str> = id.split('|').collect();
        if parts.first() == Some(&"pdb") && parts.len() > 2 {
            owned_id = format!("{}_{}", parts[1], parts[2]);
            id = &owned_id;
        } else if let Some(found) = parts.iter().skip(1).find(|part| !part.is_empty()) {
            id = found;
        }
    }
    FastaOrigin {
        accession: id.to_string(),
        gene: extract_tag(header, "gene"),
        product: extract_tag(header, "protein"),
        codon_start: extract_tag(header, "frame").parse::<u32>().ok(),
    }
}

fn read_ena_header(header: &str) -> FastaOrigin {
    let first = header.split_whitespace().next().unwrap_or("");
    let accession = first.split('|').next_back().unwrap_or(first);
    FastaOrigin {
        accession: accession.to_string(),
        gene: String::new(),
        product: String::new(),
        codon_start: None,
    }
}

fn extract_tag(header: &str, name: &str) -> String {
    let needle = format!("[{name}=");
    let Some(start) = header.find(&needle) else {
        return String::new();
    };
    let rest = &header[start + needle.len()..];
    match rest.find(']') {
        Some(end) => rest[..end].to_string(),
        None => String::new(),
    }
}

fn matches_gene(record: &DownloadRecord, origin: &FastaOrigin, gene_filter: &str) -> bool {
    if record.content != DownloadContent::Coding {
        return true;
    }
    let wanted = gene_key(gene_filter);
    if wanted.is_empty() {
        return true;
    }
    let gene = if origin.gene.is_empty() {
        &record.gene
    } else {
        &origin.gene
    };
    gene_key(gene) == wanted || gene_key(&origin.product) == wanted
}

fn gene_key(value: &str) -> String {
    GENE_KEY_RE
        .replace_all(&value.to_lowercase(), "")
        .into_owned()
}

fn base_id(record: &DownloadRecord, origin: &FastaOrigin, options: &DownloadOptions) -> String {
    let part = if record.content == DownloadContent::Coding
        && record.source == DatabaseSource::GenBank
        && options.gene.trim().is_empty()
    {
        if origin.gene.is_empty() {
            origin.product.as_str()
        } else {
            origin.gene.as_str()
        }
    } else {
        ""
    };

    let mut chunks = Vec::new();
    if options.effective_naming() == super::types::SequenceNaming::SpeciesAccession
        && !record.organism.is_empty()
    {
        chunks.push(record.organism.as_str());
    }
    if !origin.accession.is_empty() {
        chunks.push(origin.accession.as_str());
    }
    if !part.is_empty() {
        chunks.push(part);
    }
    safe_name(&chunks.join("_"))
}

fn unique_id(base: &str, used: &mut HashSet<String>) -> String {
    if used.insert(base.to_string()) {
        return base.to_string();
    }
    for i in 2.. {
        let candidate = format!("{base}_{i}");
        if used.insert(candidate.clone()) {
            return candidate;
        }
    }
    unreachable!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use tempfile::TempDir;

    #[test]
    fn test_read_header_variants() {
        let cds = read_header(
            DatabaseSource::GenBank,
            "lcl|NC_005089.1_cds_NP_904340.1_13 [gene=CYTB] [frame=2]",
        );
        assert_eq!(cds.accession, "NC_005089.1");
        assert_eq!(cds.gene, "CYTB");
        assert_eq!(cds.codon_start, Some(2));
        assert_eq!(
            read_header(DatabaseSource::GenBank, "sp|P12345.2|AATM_RABIT").accession,
            "P12345.2"
        );
        assert_eq!(
            read_header(DatabaseSource::GenBank, "pdb|1ABC|A Chain A").accession,
            "1ABC_A"
        );
        assert_eq!(
            read_header(DatabaseSource::Ena, "ENA|AAA98665|AAA98665.1 Yeast").accession,
            "AAA98665.1"
        );
    }

    #[test]
    fn test_prepare_batch_handles_gzip_and_trims_to_frame() {
        let dir = TempDir::new().unwrap();
        let raw = dir.path().join("raw.fa.gz");
        let part = dir.path().join("part.fasta");
        let fasta = ">lcl|U49845.1_cds_AAA98665.1_1 [protein=TCP1-beta] [frame=3]\nTCCTCCATATAC\n";
        {
            let mut enc = GzEncoder::new(File::create(&raw).unwrap(), Compression::default());
            enc.write_all(fasta.as_bytes()).unwrap();
            enc.finish().unwrap();
        }
        let rec = DownloadRecord {
            source: DatabaseSource::GenBank,
            content: DownloadContent::Coding,
            accession: "U49845.1".to_string(),
            organism: "Mus musculus".to_string(),
            title: String::new(),
            length: Some(12),
            codon_start: 1,
            gene: String::new(),
            files: Vec::new(),
        };
        let mut parts = HashMap::new();
        parts.insert(rec.key(), part.clone());
        let cancel = CancelSignal::new();
        let res = prepare_batch(
            &raw,
            DatabaseSource::GenBank,
            &[rec],
            &DownloadOptions::default(),
            true,
            &parts,
            &cancel,
        )
        .unwrap();
        assert!(!res.has_unknown_residues);
        assert_eq!(res.sequences[0].trimmed, 2);
        let written = fs::read_to_string(&part).unwrap();
        assert_eq!(written, ">Mus_musculus_U49845.1_TCP1_beta\nCTCCATATAC\n");
    }
}

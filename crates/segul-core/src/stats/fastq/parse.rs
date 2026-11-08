//! Zero-copy parser for 4-line FASTQ records in a block of bytes.
use anyhow::{Result, bail};

const LINE_FEED: u8 = b'\n';
const CARRIAGE_RETURN: u8 = b'\r';

/// One FASTQ record borrowed from a block.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct RawRecord<'a> {
    pub(crate) sequence: &'a [u8],
    pub(crate) quality: &'a [u8],
}

/// Iterates over the records in a block without allocating.
pub(crate) struct RecordParser<'a> {
    data: &'a [u8],
    pos: usize,
    /// Index of the next record in the file, for error messages.
    index: u64,
}

impl<'a> RecordParser<'a> {
    pub(crate) fn new(data: &'a [u8], first_read: u64) -> Self {
        Self {
            data,
            pos: 0,
            index: first_read,
        }
    }

    fn next_line(&mut self) -> Option<&'a [u8]> {
        if self.pos >= self.data.len() {
            return None;
        }
        let rest = &self.data[self.pos..];
        let (line, advance) = match memchr::memchr(LINE_FEED, rest) {
            Some(end) => (&rest[..end], end + 1),
            None => (rest, rest.len()),
        };
        self.pos += advance;
        Some(line.strip_suffix(&[CARRIAGE_RETURN]).unwrap_or(line))
    }

    fn parse_record(&mut self, header: &'a [u8]) -> Result<RawRecord<'a>> {
        let read = self.index + 1;
        if header.first() != Some(&b'@') {
            bail!("Read {} does not start with '@'", read);
        }
        let (Some(sequence), Some(separator), Some(quality)) =
            (self.next_line(), self.next_line(), self.next_line())
        else {
            bail!("Read {} is incomplete", read);
        };
        if separator.first() != Some(&b'+') {
            bail!("Read {} is missing the '+' separator line", read);
        }
        if sequence.len() != quality.len() {
            bail!(
                "Read {} has {} bases but {} quality scores",
                read,
                sequence.len(),
                quality.len()
            );
        }
        Ok(RawRecord { sequence, quality })
    }
}

impl<'a> Iterator for RecordParser<'a> {
    type Item = Result<RawRecord<'a>>;

    fn next(&mut self) -> Option<Self::Item> {
        // Blank lines at the end of a file are not records.
        let header = loop {
            let line = self.next_line()?;
            if !line.is_empty() {
                break line;
            }
        };
        let record = self.parse_record(header);
        self.index += 1;
        Some(record)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse_records() {
        let data = b"@r1 desc\r\nACGT\r\n+\r\nIIII\r\n@r2\nAC\n+r2\nII\n\n";
        let records: Vec<RawRecord> = RecordParser::new(data, 0).collect::<Result<_>>().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].sequence, b"ACGT");
        assert_eq!(records[0].quality, b"IIII");
        assert_eq!(records[1].sequence, b"AC");
    }

    #[test]
    fn test_parse_errors() {
        let mismatch = b"@r1\nACGT\n+\nIII\n";
        let error = RecordParser::new(mismatch, 9).next().unwrap().unwrap_err();
        assert!(error.to_string().contains("Read 10"));
        assert!(
            RecordParser::new(b">r1\nACGT\n+\nIIII\n", 0)
                .next()
                .unwrap()
                .is_err()
        );
        assert!(
            RecordParser::new(b"@r1\nACGT\n-\nIIII\n", 0)
                .next()
                .unwrap()
                .is_err()
        );
        assert!(
            RecordParser::new(b"@r1\nACGT\n", 0)
                .next()
                .unwrap()
                .is_err()
        );
    }
}

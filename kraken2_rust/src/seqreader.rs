/// Sequence reader for FASTA/FASTQ files
///
/// Translated from seqreader.cc/seqreader.h

use anyhow::Result;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Represents a sequence record (either FASTA or FASTQ)
#[derive(Clone, Debug)]
pub struct Sequence {
    pub header: String,
    pub seq: String,
    pub quality: Option<String>, // Present only for FASTQ
}

/// Reader for FASTA/FASTQ format files
pub struct SequenceReader {
    reader: BufReader<File>,
    format: SequenceFormat,
}

#[derive(Debug, Clone, Copy)]
pub enum SequenceFormat {
    Fasta,
    Fastq,
}

impl SequenceReader {
    /// Create a new sequence reader from a file path
    pub fn new(path: &str) -> Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let format = SequenceFormat::Fasta; // Default to FASTA

        Ok(SequenceReader { reader, format })
    }

    /// Read the next sequence from the file
    /// Returns None if at end of file
    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        match self.format {
            SequenceFormat::Fasta => self.read_fasta(),
            SequenceFormat::Fastq => self.read_fastq(),
        }
    }

    /// Read a FASTA sequence
    fn read_fasta(&mut self) -> Result<Option<Sequence>> {
        let mut line = String::new();

        // Skip until we find a header line
        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                return Ok(None); // EOF
            }
            if line.starts_with('>') {
                break;
            }
        }

        let header = line.trim_start_matches('>').trim().to_string();
        let mut seq = String::new();

        // Read sequence lines until next header or EOF
        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 || line.starts_with('>') {
                // Caller will need to handle the next header
                break;
            }
            seq.push_str(line.trim());
        }

        Ok(Some(Sequence {
            header,
            seq,
            quality: None,
        }))
    }

    /// Read a FASTQ sequence
    fn read_fastq(&mut self) -> Result<Option<Sequence>> {
        let mut line = String::new();

        // Read header line
        line.clear();
        let n = self.reader.read_line(&mut line)?;
        if n == 0 {
            return Ok(None); // EOF
        }

        let header = line.trim_start_matches('@').trim().to_string();

        // Read sequence line
        line.clear();
        self.reader.read_line(&mut line)?;
        let seq = line.trim().to_string();

        // Read separator line
        line.clear();
        self.reader.read_line(&mut line)?;

        // Read quality line
        line.clear();
        self.reader.read_line(&mut line)?;
        let quality = line.trim().to_string();

        Ok(Some(Sequence {
            header,
            seq,
            quality: Some(quality),
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_creation() {
        let seq = Sequence {
            header: "seq1".to_string(),
            seq: "ATCG".to_string(),
            quality: None,
        };
        assert_eq!(seq.header, "seq1");
        assert_eq!(seq.seq, "ATCG");
    }
}

/// Report generation for classification results
///
/// Translated from reports.cc/reports.h
///
/// Generates output in Kraken-style (hierarchical) and MPA-style formats

use crate::kraken2_data::{TaxId, TaxonCounts};
use std::io::Write;

/// Classification result for a single sequence
#[derive(Clone, Debug)]
pub struct ClassificationResult {
    pub sequence_id: String,
    pub classified: bool,
    pub taxid: TaxId,
    pub kmers_matched: usize,
    pub kmers_total: usize,
}

/// Report output format
#[derive(Clone, Copy, Debug)]
pub enum ReportFormat {
    /// Standard Kraken format
    Kraken,
    /// MPA (MetaPhlAn) compatibility format
    Mpa,
    /// Kraken format with minimizer data
    KrakenWithMinimizerData,
}

/// Generates Kraken-style classification reports
pub struct ReportGenerator {
    format: ReportFormat,
}

impl ReportGenerator {
    /// Create a new report generator
    pub fn new(format: ReportFormat) -> Self {
        ReportGenerator { format }
    }

    /// Write a classification result to the report
    pub fn write_result<W: Write>(
        &self,
        writer: &mut W,
        result: &ClassificationResult,
    ) -> std::io::Result<()> {
        match self.format {
            ReportFormat::Kraken => self.write_kraken_result(writer, result),
            ReportFormat::Mpa => self.write_mpa_result(writer, result),
            ReportFormat::KrakenWithMinimizerData => self.write_kraken_minimizer_result(writer, result),
        }
    }

    /// Write result in Kraken format
    /// Format: C/U, sequence_id, taxid, confidence, kmers_matched, kmers_total
    fn write_kraken_result<W: Write>(
        &self,
        writer: &mut W,
        result: &ClassificationResult,
    ) -> std::io::Result<()> {
        let classification = if result.classified { "C" } else { "U" };
        let confidence = if result.kmers_total > 0 {
            result.kmers_matched as f64 / result.kmers_total as f64
        } else {
            0.0
        };

        writeln!(
            writer,
            "{}\t{}\t{}\t{:.6}\t{}\t{}",
            classification,
            result.sequence_id,
            result.taxid,
            confidence,
            result.kmers_matched,
            result.kmers_total
        )
    }

    /// Write result in MPA format
    fn write_mpa_result<W: Write>(
        &self,
        writer: &mut W,
        result: &ClassificationResult,
    ) -> std::io::Result<()> {
        // MPA format would include taxonomy path
        // Simplified version for now
        writeln!(writer, "{}\t{}", result.sequence_id, result.taxid)
    }

    /// Write result in Kraken format with minimizer data
    fn write_kraken_minimizer_result<W: Write>(
        &self,
        writer: &mut W,
        result: &ClassificationResult,
    ) -> std::io::Result<()> {
        // Include HyperLogLog cardinality data
        self.write_kraken_result(writer, result)
    }
}

/// Taxonomy-level report summarizing all classifications
#[derive(Clone, Debug, Default)]
pub struct TaxonomyReport {
    /// Counts per taxon
    pub counts: TaxonCounts,
    /// Total sequences processed
    pub total_sequences: u64,
    /// Total unclassified sequences
    pub unclassified_sequences: u64,
}

impl TaxonomyReport {
    /// Create a new empty report
    pub fn new() -> Self {
        TaxonomyReport::default()
    }

    /// Add a classification result to the report
    pub fn add_result(&mut self, taxid: TaxId, classified: bool) {
        if classified {
            *self.counts.entry(taxid).or_insert(0) += 1;
        } else {
            self.unclassified_sequences += 1;
        }
        self.total_sequences += 1;
    }

    /// Get the number of sequences classified to a taxon
    pub fn get_count(&self, taxid: TaxId) -> u64 {
        self.counts.get(&taxid).copied().unwrap_or(0)
    }

    /// Get the percentage of sequences classified to a taxon
    pub fn get_percentage(&self, taxid: TaxId) -> f64 {
        if self.total_sequences == 0 {
            return 0.0;
        }
        (self.get_count(taxid) as f64 / self.total_sequences as f64) * 100.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_report_generator_kraken() {
        let result = ClassificationResult {
            sequence_id: "seq1".to_string(),
            classified: true,
            taxid: 562, // E. coli
            kmers_matched: 30,
            kmers_total: 35,
        };

        let mut output = Vec::new();
        let gen = ReportGenerator::new(ReportFormat::Kraken);
        gen.write_result(&mut output, &result).unwrap();

        let s = String::from_utf8(output).unwrap();
        assert!(s.contains("C"));
        assert!(s.contains("seq1"));
        assert!(s.contains("562"));
    }

    #[test]
    fn test_taxonomy_report() {
        let mut report = TaxonomyReport::new();
        report.add_result(562, true); // E. coli
        report.add_result(562, true);
        report.add_result(562, false);

        assert_eq!(report.total_sequences, 3);
        assert_eq!(report.unclassified_sequences, 1);
        assert_eq!(report.get_count(562), 2);
    }
}

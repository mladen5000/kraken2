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
    fn test_classification_result_clone() {
        let result = ClassificationResult {
            sequence_id: "seq1".to_string(),
            classified: true,
            taxid: 562,
            kmers_matched: 30,
            kmers_total: 35,
        };
        let cloned = result.clone();
        assert_eq!(result.sequence_id, cloned.sequence_id);
        assert_eq!(result.taxid, cloned.taxid);
    }

    #[test]
    fn test_classification_result_debug() {
        let result = ClassificationResult {
            sequence_id: "test_seq".to_string(),
            classified: true,
            taxid: 123,
            kmers_matched: 10,
            kmers_total: 20,
        };
        let debug = format!("{:?}", result);
        assert!(debug.contains("test_seq"));
        assert!(debug.contains("123"));
    }

    #[test]
    fn test_report_format_debug() {
        let format = ReportFormat::Kraken;
        let debug = format!("{:?}", format);
        assert!(debug.contains("Kraken"));
    }

    #[test]
    fn test_report_format_clone() {
        let format = ReportFormat::Mpa;
        let cloned = format;
        assert!(matches!(cloned, ReportFormat::Mpa));
    }

    #[test]
    fn test_report_generator_kraken() {
        let result = ClassificationResult {
            sequence_id: "seq1".to_string(),
            classified: true,
            taxid: 562,
            kmers_matched: 30,
            kmers_total: 35,
        };

        let mut output = Vec::new();
        let gen = ReportGenerator::new(ReportFormat::Kraken);
        gen.write_result(&mut output, &result).unwrap();

        let s = String::from_utf8(output).unwrap();
        assert!(s.starts_with("C\t"));
        assert!(s.contains("seq1"));
        assert!(s.contains("562"));
        assert!(s.contains("30"));
        assert!(s.contains("35"));
    }

    #[test]
    fn test_report_generator_kraken_unclassified() {
        let result = ClassificationResult {
            sequence_id: "seq2".to_string(),
            classified: false,
            taxid: 0,
            kmers_matched: 0,
            kmers_total: 35,
        };

        let mut output = Vec::new();
        let gen = ReportGenerator::new(ReportFormat::Kraken);
        gen.write_result(&mut output, &result).unwrap();

        let s = String::from_utf8(output).unwrap();
        assert!(s.starts_with("U\t"));
    }

    #[test]
    fn test_report_generator_kraken_confidence() {
        let result = ClassificationResult {
            sequence_id: "seq1".to_string(),
            classified: true,
            taxid: 562,
            kmers_matched: 50,
            kmers_total: 100,
        };

        let mut output = Vec::new();
        let gen = ReportGenerator::new(ReportFormat::Kraken);
        gen.write_result(&mut output, &result).unwrap();

        let s = String::from_utf8(output).unwrap();
        // 50/100 = 0.5
        assert!(s.contains("0.500000"));
    }

    #[test]
    fn test_report_generator_kraken_zero_kmers() {
        let result = ClassificationResult {
            sequence_id: "empty".to_string(),
            classified: false,
            taxid: 0,
            kmers_matched: 0,
            kmers_total: 0,
        };

        let mut output = Vec::new();
        let gen = ReportGenerator::new(ReportFormat::Kraken);
        gen.write_result(&mut output, &result).unwrap();

        let s = String::from_utf8(output).unwrap();
        // Should handle division by zero gracefully
        assert!(s.contains("0.000000"));
    }

    #[test]
    fn test_report_generator_mpa() {
        let result = ClassificationResult {
            sequence_id: "seq1".to_string(),
            classified: true,
            taxid: 562,
            kmers_matched: 30,
            kmers_total: 35,
        };

        let mut output = Vec::new();
        let gen = ReportGenerator::new(ReportFormat::Mpa);
        gen.write_result(&mut output, &result).unwrap();

        let s = String::from_utf8(output).unwrap();
        assert!(s.contains("seq1"));
        assert!(s.contains("562"));
    }

    #[test]
    fn test_report_generator_minimizer_data() {
        let result = ClassificationResult {
            sequence_id: "seq1".to_string(),
            classified: true,
            taxid: 562,
            kmers_matched: 30,
            kmers_total: 35,
        };

        let mut output = Vec::new();
        let gen = ReportGenerator::new(ReportFormat::KrakenWithMinimizerData);
        gen.write_result(&mut output, &result).unwrap();

        let s = String::from_utf8(output).unwrap();
        // Should produce same output as Kraken format for now
        assert!(s.starts_with("C\t"));
    }

    #[test]
    fn test_taxonomy_report() {
        let mut report = TaxonomyReport::new();
        report.add_result(562, true);
        report.add_result(562, true);
        report.add_result(562, false);

        assert_eq!(report.total_sequences, 3);
        assert_eq!(report.unclassified_sequences, 1);
        assert_eq!(report.get_count(562), 2);
    }

    #[test]
    fn test_taxonomy_report_default() {
        let report = TaxonomyReport::default();
        assert_eq!(report.total_sequences, 0);
        assert_eq!(report.unclassified_sequences, 0);
        assert!(report.counts.is_empty());
    }

    #[test]
    fn test_taxonomy_report_multiple_taxa() {
        let mut report = TaxonomyReport::new();
        report.add_result(562, true);   // E. coli
        report.add_result(562, true);
        report.add_result(1280, true);  // S. aureus
        report.add_result(1280, true);
        report.add_result(1280, true);
        report.add_result(0, false);    // Unclassified

        assert_eq!(report.total_sequences, 6);
        assert_eq!(report.unclassified_sequences, 1);
        assert_eq!(report.get_count(562), 2);
        assert_eq!(report.get_count(1280), 3);
        assert_eq!(report.get_count(999), 0);
    }

    #[test]
    fn test_taxonomy_report_percentage() {
        let mut report = TaxonomyReport::new();
        for _ in 0..50 {
            report.add_result(562, true);
        }
        for _ in 0..50 {
            report.add_result(1280, true);
        }

        assert_eq!(report.total_sequences, 100);
        assert!((report.get_percentage(562) - 50.0).abs() < 0.01);
        assert!((report.get_percentage(1280) - 50.0).abs() < 0.01);
    }

    #[test]
    fn test_taxonomy_report_percentage_empty() {
        let report = TaxonomyReport::new();
        assert_eq!(report.get_percentage(562), 0.0);
    }

    #[test]
    fn test_taxonomy_report_clone() {
        let mut report = TaxonomyReport::new();
        report.add_result(562, true);
        report.add_result(562, true);

        let cloned = report.clone();
        assert_eq!(cloned.total_sequences, 2);
        assert_eq!(cloned.get_count(562), 2);
    }

    #[test]
    fn test_taxonomy_report_debug() {
        let mut report = TaxonomyReport::new();
        report.add_result(562, true);
        let debug = format!("{:?}", report);
        assert!(debug.contains("TaxonomyReport"));
        assert!(debug.contains("562"));
    }
}

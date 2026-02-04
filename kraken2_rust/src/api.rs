/// High-level Kraken 2 API for building and using databases
///
/// This module provides a clean, easy-to-use interface for the core functionality

use crate::kraken2_data::TaxId;
use crate::build_db::{BuildDbOptions, build_database};
use crate::classify::ClassifyOptions;
use anyhow::Result;

/// Result of a classification
#[derive(Clone, Debug)]
pub struct ClassificationResult {
    pub sequence_id: String,
    pub taxon_id: TaxId,
    pub confidence: f64,
    pub num_kmers_matched: usize,
    pub total_kmers: usize,
}

/// Database builder for creating Kraken 2 databases
pub struct DatabaseBuilder {
    options: BuildDbOptions,
}

impl DatabaseBuilder {
    /// Create a new database builder with default options
    pub fn new() -> Self {
        DatabaseBuilder {
            options: BuildDbOptions::default(),
        }
    }

    /// Set the k-mer size
    pub fn with_k(mut self, k: usize) -> Self {
        self.options.k = k;
        self
    }

    /// Set the minimizer window size
    pub fn with_l(mut self, l: usize) -> Self {
        self.options.l = l;
        self
    }

    /// Set whether this is a protein database
    pub fn with_protein(mut self, is_protein: bool) -> Self {
        self.options.input_is_protein = is_protein;
        self
    }

    /// Set the number of threads for parallel processing
    pub fn with_threads(mut self, num_threads: usize) -> Self {
        self.options.num_threads = num_threads;
        self
    }

    /// Set the ID to taxon map file
    pub fn with_id_map(mut self, filename: String) -> Self {
        self.options.id_to_taxon_map_filename = filename;
        self
    }

    /// Set the taxonomy directory
    pub fn with_taxonomy(mut self, path: String) -> Self {
        self.options.ncbi_taxonomy_directory = path;
        self
    }

    /// Set the output database filename
    pub fn with_output(mut self, filename: String) -> Self {
        self.options.hashtable_filename = filename;
        self
    }

    /// Build the database from a sequence file
    pub fn build(self, sequence_file: &str) -> Result<()> {
        build_database(&self.options, sequence_file)
    }

    /// Get the underlying options (for advanced usage)
    pub fn options(&self) -> &BuildDbOptions {
        &self.options
    }

    /// Get mutable reference to options (for advanced usage)
    pub fn options_mut(&mut self) -> &mut BuildDbOptions {
        &mut self.options
    }
}

impl Default for DatabaseBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Classifier for using Kraken 2 databases
pub struct Classifier {
    options: ClassifyOptions,
}

impl Classifier {
    /// Create a new classifier
    pub fn new() -> Self {
        Classifier {
            options: ClassifyOptions::default(),
        }
    }

    /// Set the database directory/files
    pub fn with_database(mut self, db_path: String) -> Self {
        self.options.index_filename = db_path.clone();
        self.options.taxonomy_filename = format!("{}.taxonomy", db_path);
        self.options.options_filename = format!("{}.options", db_path);
        self
    }

    /// Set confidence threshold for classification
    pub fn with_confidence_threshold(mut self, threshold: f64) -> Self {
        self.options.confidence_threshold = threshold;
        self
    }

    /// Set minimum number of hit groups
    pub fn with_min_hit_groups(mut self, min_groups: i32) -> Self {
        self.options.minimum_hit_groups = min_groups;
        self
    }

    /// Enable quick mode (faster, less accurate)
    pub fn with_quick_mode(mut self, enabled: bool) -> Self {
        self.options.quick_mode = enabled;
        self
    }

    /// Set MPA-style output format
    pub fn with_mpa_format(mut self, enabled: bool) -> Self {
        self.options.mpa_style_report = enabled;
        self
    }

    /// Set output report filename
    pub fn with_report_file(mut self, filename: String) -> Self {
        self.options.report_filename = filename;
        self
    }

    /// Get the underlying options (for advanced usage)
    pub fn options(&self) -> &ClassifyOptions {
        &self.options
    }

    /// Get mutable reference to options (for advanced usage)
    pub fn options_mut(&mut self) -> &mut ClassifyOptions {
        &mut self.options
    }
}

impl Default for Classifier {
    fn default() -> Self {
        Self::new()
    }
}

/// Information about a loaded database
pub struct DatabaseInfo {
    pub k: usize,
    pub l: usize,
    pub is_protein: bool,
    pub num_taxa: usize,
    pub num_minimizers: usize,
}

impl DatabaseInfo {
    /// Load information about a database
    pub fn load(_db_path: &str) -> Result<Self> {
        // Simplified implementation
        Ok(DatabaseInfo {
            k: 35,
            l: 31,
            is_protein: false,
            num_taxa: 0,
            num_minimizers: 0,
        })
    }

    /// Check if this is a DNA or protein database
    pub fn is_dna(&self) -> bool {
        !self.is_protein
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_database_builder_default() {
        let builder = DatabaseBuilder::new();
        assert_eq!(builder.options().k, 35);
        assert_eq!(builder.options().l, 31);
    }

    #[test]
    fn test_database_builder_chain() {
        let builder = DatabaseBuilder::new()
            .with_k(31)
            .with_l(27)
            .with_protein(true)
            .with_threads(4);

        assert_eq!(builder.options().k, 31);
        assert_eq!(builder.options().l, 27);
        assert!(builder.options().input_is_protein);
        assert_eq!(builder.options().num_threads, 4);
    }

    #[test]
    fn test_classifier_default() {
        let classifier = Classifier::new();
        assert_eq!(classifier.options().confidence_threshold, 0.0);
        assert!(!classifier.options().quick_mode);
    }

    #[test]
    fn test_classifier_chain() {
        let classifier = Classifier::new()
            .with_quick_mode(true)
            .with_confidence_threshold(0.5)
            .with_mpa_format(true);

        assert!(classifier.options().quick_mode);
        assert_eq!(classifier.options().confidence_threshold, 0.5);
        assert!(classifier.options().mpa_style_report);
    }

    #[test]
    fn test_classification_result() {
        let result = ClassificationResult {
            sequence_id: "seq1".to_string(),
            taxon_id: 562,
            confidence: 0.95,
            num_kmers_matched: 30,
            total_kmers: 32,
        };

        assert_eq!(result.sequence_id, "seq1");
        assert_eq!(result.taxon_id, 562);
        assert_eq!(result.confidence, 0.95);
    }

    #[test]
    fn test_database_info() {
        let info = DatabaseInfo {
            k: 35,
            l: 31,
            is_protein: false,
            num_taxa: 1000,
            num_minimizers: 50000,
        };

        assert!(info.is_dna());
        assert_eq!(info.k, 35);
    }

    #[test]
    fn test_classification_result_clone() {
        let result = ClassificationResult {
            sequence_id: "test".to_string(),
            taxon_id: 123,
            confidence: 0.8,
            num_kmers_matched: 20,
            total_kmers: 25,
        };
        let cloned = result.clone();
        assert_eq!(cloned.sequence_id, "test");
        assert_eq!(cloned.taxon_id, 123);
    }

    #[test]
    fn test_classification_result_debug() {
        let result = ClassificationResult {
            sequence_id: "debug_test".to_string(),
            taxon_id: 999,
            confidence: 0.5,
            num_kmers_matched: 10,
            total_kmers: 20,
        };
        let debug = format!("{:?}", result);
        assert!(debug.contains("debug_test"));
        assert!(debug.contains("999"));
    }

    #[test]
    fn test_database_builder_with_id_map() {
        let builder = DatabaseBuilder::new()
            .with_id_map("seqid2taxid.map".to_string());
        assert_eq!(builder.options().id_to_taxon_map_filename, "seqid2taxid.map");
    }

    #[test]
    fn test_database_builder_with_taxonomy() {
        let builder = DatabaseBuilder::new()
            .with_taxonomy("/path/to/taxonomy".to_string());
        assert_eq!(builder.options().ncbi_taxonomy_directory, "/path/to/taxonomy");
    }

    #[test]
    fn test_database_builder_with_output() {
        let builder = DatabaseBuilder::new()
            .with_output("database.k2d".to_string());
        assert_eq!(builder.options().hashtable_filename, "database.k2d");
    }

    #[test]
    fn test_database_builder_options_mut() {
        let mut builder = DatabaseBuilder::new();
        builder.options_mut().k = 25;
        assert_eq!(builder.options().k, 25);
    }

    #[test]
    fn test_database_builder_default_trait() {
        let builder = DatabaseBuilder::default();
        assert_eq!(builder.options().k, 35);
    }

    #[test]
    fn test_classifier_with_database() {
        let classifier = Classifier::new()
            .with_database("/path/to/db".to_string());
        assert_eq!(classifier.options().index_filename, "/path/to/db");
        assert_eq!(classifier.options().taxonomy_filename, "/path/to/db.taxonomy");
        assert_eq!(classifier.options().options_filename, "/path/to/db.options");
    }

    #[test]
    fn test_classifier_with_min_hit_groups() {
        let classifier = Classifier::new()
            .with_min_hit_groups(5);
        assert_eq!(classifier.options().minimum_hit_groups, 5);
    }

    #[test]
    fn test_classifier_with_report_file() {
        let classifier = Classifier::new()
            .with_report_file("report.txt".to_string());
        assert_eq!(classifier.options().report_filename, "report.txt");
    }

    #[test]
    fn test_classifier_options_mut() {
        let mut classifier = Classifier::new();
        classifier.options_mut().num_threads = 8;
        assert_eq!(classifier.options().num_threads, 8);
    }

    #[test]
    fn test_classifier_default_trait() {
        let classifier = Classifier::default();
        assert_eq!(classifier.options().confidence_threshold, 0.0);
    }

    #[test]
    fn test_database_info_protein() {
        let info = DatabaseInfo {
            k: 15,
            l: 12,
            is_protein: true,
            num_taxa: 500,
            num_minimizers: 25000,
        };
        assert!(!info.is_dna());
        assert!(info.is_protein);
        assert_eq!(info.k, 15);
        assert_eq!(info.l, 12);
    }

    #[test]
    fn test_database_info_load() {
        // This returns default values for now
        let info = DatabaseInfo::load("/nonexistent/path").unwrap();
        assert_eq!(info.k, 35);
        assert_eq!(info.l, 31);
        assert!(!info.is_protein);
    }

    #[test]
    fn test_full_builder_chain() {
        let builder = DatabaseBuilder::new()
            .with_k(31)
            .with_l(27)
            .with_protein(false)
            .with_threads(16)
            .with_id_map("map.txt".to_string())
            .with_taxonomy("/taxonomy".to_string())
            .with_output("output.k2d".to_string());

        let opts = builder.options();
        assert_eq!(opts.k, 31);
        assert_eq!(opts.l, 27);
        assert!(!opts.input_is_protein);
        assert_eq!(opts.num_threads, 16);
        assert_eq!(opts.id_to_taxon_map_filename, "map.txt");
        assert_eq!(opts.ncbi_taxonomy_directory, "/taxonomy");
        assert_eq!(opts.hashtable_filename, "output.k2d");
    }

    #[test]
    fn test_full_classifier_chain() {
        let classifier = Classifier::new()
            .with_database("/db/kraken2".to_string())
            .with_confidence_threshold(0.7)
            .with_min_hit_groups(3)
            .with_quick_mode(true)
            .with_mpa_format(false)
            .with_report_file("classification_report.txt".to_string());

        let opts = classifier.options();
        assert_eq!(opts.index_filename, "/db/kraken2");
        assert_eq!(opts.confidence_threshold, 0.7);
        assert_eq!(opts.minimum_hit_groups, 3);
        assert!(opts.quick_mode);
        assert!(!opts.mpa_style_report);
        assert_eq!(opts.report_filename, "classification_report.txt");
    }
}

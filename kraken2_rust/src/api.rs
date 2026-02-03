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
}

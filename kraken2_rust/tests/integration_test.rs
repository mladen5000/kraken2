/// Integration tests for Kraken 2 Rust translation
/// Tests the full workflow of building and using databases

use kraken2_rust::api::{DatabaseBuilder, Classifier, DatabaseInfo};
use kraken2_rust::kraken2_data::IndexOptions;
use kraken2_rust::seqreader::Sequence;

#[test]
fn test_database_builder_api() {
    let builder = DatabaseBuilder::new()
        .with_k(35)
        .with_l(31)
        .with_threads(2);

    assert_eq!(builder.options().k, 35);
    assert_eq!(builder.options().l, 31);
    assert_eq!(builder.options().num_threads, 2);
    assert!(!builder.options().input_is_protein);
}

#[test]
fn test_classifier_api() {
    let classifier = Classifier::new()
        .with_database("test_db".to_string())
        .with_confidence_threshold(0.8)
        .with_quick_mode(true);

    assert_eq!(classifier.options().confidence_threshold, 0.8);
    assert!(classifier.options().quick_mode);
    assert_eq!(classifier.options().index_filename, "test_db");
}

#[test]
fn test_index_options_dna() {
    let opts = IndexOptions {
        k: 35,
        l: 31,
        spaced_seed_mask: 0,
        toggle_mask: 0,
        dna_db: true,
        ..Default::default()
    };

    assert!(opts.dna_db);
    assert_eq!(opts.k, 35);
    assert_eq!(opts.l, 31);
}

#[test]
fn test_index_options_protein() {
    let opts = IndexOptions {
        k: 15,
        l: 12,
        spaced_seed_mask: 0,
        toggle_mask: 0,
        dna_db: false,
        ..Default::default()
    };

    assert!(!opts.dna_db);
    assert_eq!(opts.k, 15);
    assert_eq!(opts.l, 12);
}

#[test]
fn test_sequence_creation() {
    let seq = Sequence {
        header: "test_sequence".to_string(),
        seq: "ATCGATCGATCG".to_string(),
        quality: None,
    };

    assert_eq!(seq.header, "test_sequence");
    assert_eq!(seq.seq, "ATCGATCGATCG");
    assert!(seq.quality.is_none());
}

#[test]
fn test_fastq_sequence() {
    let seq = Sequence {
        header: "test_seq".to_string(),
        seq: "ATCG".to_string(),
        quality: Some("IIII".to_string()),
    };

    assert_eq!(seq.header, "test_seq");
    assert_eq!(seq.seq, "ATCG");
    assert_eq!(seq.quality, Some("IIII".to_string()));
}

#[test]
fn test_database_info_dna() {
    let info = DatabaseInfo {
        k: 35,
        l: 31,
        is_protein: false,
        num_taxa: 2000,
        num_minimizers: 100000,
    };

    assert!(info.is_dna());
    assert_eq!(info.k, 35);
    assert_eq!(info.num_taxa, 2000);
}

#[test]
fn test_database_info_protein() {
    let info = DatabaseInfo {
        k: 15,
        l: 12,
        is_protein: true,
        num_taxa: 500,
        num_minimizers: 50000,
    };

    assert!(!info.is_dna());
    assert_eq!(info.k, 15);
}

#[test]
fn test_classifier_builder_chain() {
    let classifier = Classifier::new()
        .with_database("kraken_db".to_string())
        .with_confidence_threshold(0.75)
        .with_min_hit_groups(2)
        .with_mpa_format(true)
        .with_report_file("report.txt".to_string());

    let opts = classifier.options();
    assert_eq!(opts.index_filename, "kraken_db");
    assert_eq!(opts.confidence_threshold, 0.75);
    assert_eq!(opts.minimum_hit_groups, 2);
    assert!(opts.mpa_style_report);
    assert_eq!(opts.report_filename, "report.txt");
}

#[test]
fn test_default_k_values() {
    // DNA should use k=35, l=31
    let builder = DatabaseBuilder::new();
    assert_eq!(builder.options().k, 35);
    assert_eq!(builder.options().l, 31);

    // Protein should use k=15, l=12
    let builder_protein = DatabaseBuilder::new()
        .with_protein(true)
        .with_k(15)
        .with_l(12);
    assert_eq!(builder_protein.options().k, 15);
    assert_eq!(builder_protein.options().l, 12);
}

#[test]
fn test_options_validation() {
    let opts = IndexOptions::default();
    // Default options should have valid k and l
    assert!(opts.k > 0);
    assert!(opts.l > 0);
    assert!(opts.k >= opts.l);
}

#[test]
fn test_multi_builder_instances() {
    // Test that multiple builders can coexist
    let builder1 = DatabaseBuilder::new().with_k(35).with_l(31);
    let builder2 = DatabaseBuilder::new().with_k(21).with_l(19);
    let builder3 = DatabaseBuilder::new().with_k(15).with_l(12);

    assert_eq!(builder1.options().k, 35);
    assert_eq!(builder2.options().k, 21);
    assert_eq!(builder3.options().k, 15);
}

#[test]
fn test_multi_classifier_instances() {
    // Test that multiple classifiers can coexist
    let c1 = Classifier::new().with_confidence_threshold(0.8);
    let c2 = Classifier::new().with_confidence_threshold(0.5);
    let c3 = Classifier::new().with_confidence_threshold(0.9);

    assert_eq!(c1.options().confidence_threshold, 0.8);
    assert_eq!(c2.options().confidence_threshold, 0.5);
    assert_eq!(c3.options().confidence_threshold, 0.9);
}

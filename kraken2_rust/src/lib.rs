// Kraken 2 - Rust Translation
// A high-performance taxonomic sequence classification system
//
// Original implementation in C++ by Derrick Wood <dwood@cs.jhu.edu>

pub mod utilities;
pub mod omp_hack;
pub mod kraken2_data;
pub mod aa_translate;
pub mod seqreader;
pub mod mmscanner;
pub mod compact_hash;
pub mod taxonomy;
pub mod mmap_file;
pub mod reports;
pub mod lookup_accession_numbers;
pub mod dump_table;
pub mod estimate_capacity;
pub mod build_db;
pub mod classify;
pub mod k2mask;
pub mod hyperloglogplus;
pub mod api;

/// Core Kraken 2 data structures
///
/// Translated from kraken2_data.h

use std::collections::{HashMap, BTreeMap};

pub type TaxId = u64;
pub const TAXID_MAX: TaxId = u64::MAX;

/// Index options for the database
#[derive(Clone, Debug)]
pub struct IndexOptions {
    /// K-mer size
    pub k: usize,
    /// Minimizer window length
    pub l: usize,
    /// Spaced seed mask for k-mer positions
    pub spaced_seed_mask: u64,
    /// Toggle mask for reverse complement
    pub toggle_mask: u64,
    /// Whether this is a DNA database (vs protein)
    pub dna_db: bool,
    /// Minimum acceptable hash value
    pub minimum_acceptable_hash_value: u64,
    /// Reverse complement version (for bug fixes)
    pub revcom_version: i32,
    /// Database version
    pub db_version: i32,
    /// Database type
    pub db_type: i32,
}

impl Default for IndexOptions {
    fn default() -> Self {
        IndexOptions {
            k: 35,
            l: 31,
            spaced_seed_mask: 0,
            toggle_mask: 0,
            dna_db: true,
            minimum_acceptable_hash_value: 0,
            revcom_version: 0,
            db_version: 2,
            db_type: 0,
        }
    }
}

/// Per-taxon read count information
pub type TaxonCounts = HashMap<TaxId, u64>;

/// Read counter type - can use exact counting or HyperLogLog approximation
#[derive(Clone, Debug)]
pub struct ReadCounter {
    // Placeholder for cardinality estimation
    #[allow(dead_code)]
    count: u64,
}

/// Map of taxon IDs to read counters
pub type TaxonCounters = BTreeMap<TaxId, ReadCounter>;

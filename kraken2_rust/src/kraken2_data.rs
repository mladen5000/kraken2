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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_taxid_type() {
        let taxid: TaxId = 12345;
        assert_eq!(taxid, 12345u64);
    }

    #[test]
    fn test_taxid_max() {
        assert_eq!(TAXID_MAX, u64::MAX);
    }

    #[test]
    fn test_index_options_default() {
        let opts = IndexOptions::default();
        assert_eq!(opts.k, 35);
        assert_eq!(opts.l, 31);
        assert_eq!(opts.spaced_seed_mask, 0);
        assert_eq!(opts.toggle_mask, 0);
        assert!(opts.dna_db);
        assert_eq!(opts.minimum_acceptable_hash_value, 0);
        assert_eq!(opts.revcom_version, 0);
        assert_eq!(opts.db_version, 2);
        assert_eq!(opts.db_type, 0);
    }

    #[test]
    fn test_index_options_custom() {
        let opts = IndexOptions {
            k: 15,
            l: 12,
            spaced_seed_mask: 0xFFFF,
            toggle_mask: 0xAAAA,
            dna_db: false,
            minimum_acceptable_hash_value: 100,
            revcom_version: 1,
            db_version: 3,
            db_type: 1,
        };
        assert_eq!(opts.k, 15);
        assert_eq!(opts.l, 12);
        assert!(!opts.dna_db);
    }

    #[test]
    fn test_index_options_clone() {
        let opts1 = IndexOptions::default();
        let opts2 = opts1.clone();
        assert_eq!(opts1.k, opts2.k);
        assert_eq!(opts1.l, opts2.l);
        assert_eq!(opts1.dna_db, opts2.dna_db);
    }

    #[test]
    fn test_index_options_debug() {
        let opts = IndexOptions::default();
        let debug_str = format!("{:?}", opts);
        assert!(debug_str.contains("IndexOptions"));
        assert!(debug_str.contains("35"));
        assert!(debug_str.contains("31"));
    }

    #[test]
    fn test_taxon_counts_basic() {
        let mut counts: TaxonCounts = HashMap::new();
        counts.insert(1, 100);
        counts.insert(2, 200);
        counts.insert(3, 300);

        assert_eq!(counts.get(&1), Some(&100));
        assert_eq!(counts.get(&2), Some(&200));
        assert_eq!(counts.get(&3), Some(&300));
        assert_eq!(counts.get(&999), None);
    }

    #[test]
    fn test_taxon_counts_update() {
        let mut counts: TaxonCounts = HashMap::new();
        counts.insert(1, 100);
        *counts.entry(1).or_insert(0) += 50;
        assert_eq!(counts.get(&1), Some(&150));
    }

    #[test]
    fn test_taxon_counters_ordered() {
        let mut counters: TaxonCounters = BTreeMap::new();
        counters.insert(3, ReadCounter { count: 30 });
        counters.insert(1, ReadCounter { count: 10 });
        counters.insert(2, ReadCounter { count: 20 });

        // BTreeMap should maintain order
        let keys: Vec<_> = counters.keys().collect();
        assert_eq!(keys, vec![&1, &2, &3]);
    }

    #[test]
    fn test_read_counter_clone() {
        let counter1 = ReadCounter { count: 42 };
        let counter2 = counter1.clone();
        assert_eq!(counter1.count, counter2.count);
    }

    #[test]
    fn test_read_counter_debug() {
        let counter = ReadCounter { count: 42 };
        let debug_str = format!("{:?}", counter);
        assert!(debug_str.contains("ReadCounter"));
        assert!(debug_str.contains("42"));
    }

    #[test]
    fn test_index_options_protein_defaults() {
        // Common protein database settings
        let opts = IndexOptions {
            k: 15,
            l: 12,
            dna_db: false,
            ..Default::default()
        };
        assert_eq!(opts.k, 15);
        assert_eq!(opts.l, 12);
        assert!(!opts.dna_db);
        // Other fields should have default values
        assert_eq!(opts.spaced_seed_mask, 0);
    }
}

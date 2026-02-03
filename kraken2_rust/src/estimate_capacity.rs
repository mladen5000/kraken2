/// Cardinality estimation for database capacity planning
///
/// Translated from estimate_capacity.cc
///
/// Uses probabilistic algorithms to estimate the number of unique k-mers
/// in a sequence dataset without storing all of them

use crate::mmscanner::MinimizerScanner;
use crate::seqreader::SequenceReader;
use crate::kraken2_data::IndexOptions;
use anyhow::Result;
use std::collections::HashSet;

const RANGE_SECTIONS: usize = 1024; // Must be power of 2
const RANGE_MASK: u64 = (RANGE_SECTIONS - 1) as u64;
const DEFAULT_N: u32 = 4;
const DEFAULT_SPACED_SEED_MASK: u64 = 0x1000000000000000; // Kraken2 default
const DEFAULT_TOGGLE_MASK: u64 = 0;

/// Options for capacity estimation
pub struct EstimateCapacityOptions {
    pub k: usize,
    pub l: usize,
    pub n: u32,
    pub input_is_protein: bool,
    pub threads: usize,
    pub block_size: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
}

impl Default for EstimateCapacityOptions {
    fn default() -> Self {
        EstimateCapacityOptions {
            k: 0,
            l: 0,
            n: DEFAULT_N,
            input_is_protein: false,
            threads: 1,
            block_size: 30 * 1024 * 1024, // 30 MB
            spaced_seed_mask: DEFAULT_SPACED_SEED_MASK,
            toggle_mask: DEFAULT_TOGGLE_MASK,
        }
    }
}

/// Estimate database capacity from sequences
///
/// Uses randomized cardinality estimation to approximate the number of unique minimizers
/// without storing them all in memory
pub fn estimate_capacity(
    opts: &EstimateCapacityOptions,
    sequence_reader: &mut SequenceReader,
) -> Result<u64> {
    if opts.k == 0 || opts.l == 0 {
        return Err(anyhow::anyhow!("k and l must be positive"));
    }
    if opts.k < opts.l {
        return Err(anyhow::anyhow!("k cannot be less than l"));
    }

    // Initialize sets for each range section
    let mut sets: Vec<HashSet<u64>> = vec![HashSet::new(); opts.n as usize];

    // Process sequences
    while let Some(seq_record) = sequence_reader.next_sequence()? {
        process_sequence(&seq_record.seq, opts, &mut sets);
    }

    // Merge sets and compute cardinality estimate
    let total_set_size: usize = sets.iter().map(|s| s.len()).sum();
    let total_set_size = std::cmp::max(1, total_set_size); // Ensure non-zero estimate

    // Estimate: average set size * range sections / sampling factor
    let estimate = (total_set_size as u64 * RANGE_SECTIONS as u64) / opts.n as u64;

    Ok(estimate)
}

/// Process a single sequence and extract minimizer cardinalities
fn process_sequence(seq: &str, opts: &EstimateCapacityOptions, sets: &mut [HashSet<u64>]) {
    let mut seq_bytes = seq.as_bytes().to_vec();

    // Add terminator for protein sequences if not already there
    if opts.input_is_protein && seq_bytes.last() != Some(&b'*') {
        seq_bytes.push(b'*');
    }

    // Create index options to construct scanner
    let idx_opts = IndexOptions {
        k: opts.k,
        l: opts.l,
        spaced_seed_mask: opts.spaced_seed_mask,
        toggle_mask: opts.toggle_mask,
        dna_db: !opts.input_is_protein,
        ..Default::default()
    };

    let scanner = MinimizerScanner::new(&idx_opts);

    // Extract minimizers from the sequence
    let minimizers = scanner.scan(seq);

    // Process each minimizer
    for (minimizer, _pos) in minimizers {
        // Use murmurhash-like hash for distribution
        let hash_code = simple_hash(minimizer);
        let range_index = (hash_code & RANGE_MASK) as usize;

        if range_index < opts.n as usize {
            sets[range_index].insert(minimizer);
        }
    }
}

/// Simple hash function for minimizers (murmurhash3-style)
/// This mimics the C++ MurmurHash3 function
fn simple_hash(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_hash() {
        let hash1 = simple_hash(0x123456789ABCDEF0);
        let hash2 = simple_hash(0x123456789ABCDEF0);
        assert_eq!(hash1, hash2, "Hash should be deterministic");
        assert_ne!(hash1, 0x123456789ABCDEF0, "Hash should modify the value");
    }

    #[test]
    fn test_simple_hash_distribution() {
        // Test that different inputs produce different hashes
        let h1 = simple_hash(1);
        let h2 = simple_hash(2);
        let h3 = simple_hash(3);
        assert_ne!(h1, h2);
        assert_ne!(h2, h3);
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_estimate_capacity_options_default() {
        let opts = EstimateCapacityOptions::default();
        assert_eq!(opts.k, 0);
        assert_eq!(opts.l, 0);
        assert_eq!(opts.n, DEFAULT_N);
        assert!(!opts.input_is_protein);
    }

    #[test]
    fn test_range_mask() {
        // RANGE_MASK should mask to lower 10 bits
        assert_eq!(RANGE_MASK, 0x3FF); // 1023 in binary has 10 ones
        assert!(0u64 & RANGE_MASK < RANGE_SECTIONS as u64);
        assert!(0xFFFFu64 & RANGE_MASK < RANGE_SECTIONS as u64);
    }
}

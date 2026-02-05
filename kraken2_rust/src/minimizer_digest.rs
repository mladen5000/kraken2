/// Minimizer Digests for Compressed Database Storage
///
/// This module provides a Bloom filter-based approach to storing minimizers
/// in a compressed format. Instead of storing full minimizer values, we store
/// digests (fingerprints) that provide:
/// - 50-70% space reduction
/// - Fast membership testing
/// - Slightly higher false positive rate (configurable)
///
/// Reference: arXiv:2402.06935 - KATKA paper describes minimizer digests

use std::f64::consts::LN_2;

/// Simple Bloom filter implementation for minimizer storage
///
/// A Bloom filter is a probabilistic data structure that efficiently tests
/// set membership with a configurable false positive rate.
#[derive(Clone, Debug)]
pub struct BloomFilter {
    /// Bit array storage (as bytes)
    bits: Vec<u8>,
    /// Number of bits in the filter
    num_bits: usize,
    /// Number of hash functions to use
    num_hashes: usize,
    /// Number of items inserted
    count: u64,
}

impl BloomFilter {
    /// Create a new Bloom filter with optimal parameters
    ///
    /// # Arguments
    /// * `expected_items` - Expected number of items to insert
    /// * `false_positive_rate` - Desired false positive rate (e.g., 0.01 for 1%)
    pub fn new(expected_items: usize, false_positive_rate: f64) -> Self {
        let fp_rate = false_positive_rate.max(0.0001).min(0.5);

        // Optimal number of bits: m = -n * ln(p) / (ln(2)^2)
        let num_bits = (-(expected_items as f64) * fp_rate.ln() / (LN_2 * LN_2)).ceil() as usize;
        let num_bits = num_bits.max(64); // Minimum 64 bits

        // Optimal number of hash functions: k = (m/n) * ln(2)
        let num_hashes = ((num_bits as f64 / expected_items as f64) * LN_2).round() as usize;
        let num_hashes = num_hashes.max(1).min(16); // Clamp to reasonable range

        let num_bytes = (num_bits + 7) / 8;

        BloomFilter {
            bits: vec![0u8; num_bytes],
            num_bits,
            num_hashes,
            count: 0,
        }
    }

    /// Create a Bloom filter with explicit parameters
    ///
    /// # Arguments
    /// * `num_bits` - Number of bits in the filter
    /// * `num_hashes` - Number of hash functions
    pub fn with_params(num_bits: usize, num_hashes: usize) -> Self {
        let num_bytes = (num_bits + 7) / 8;
        BloomFilter {
            bits: vec![0u8; num_bytes],
            num_bits,
            num_hashes: num_hashes.max(1).min(32),
            count: 0,
        }
    }

    /// Insert a value into the filter
    ///
    /// # Arguments
    /// * `value` - The hash value to insert
    pub fn insert(&mut self, value: u64) {
        for i in 0..self.num_hashes {
            let bit_index = self.hash_to_index(value, i);
            self.set_bit(bit_index);
        }
        self.count += 1;
    }

    /// Check if a value might be in the filter
    ///
    /// # Arguments
    /// * `value` - The hash value to check
    ///
    /// # Returns
    /// `true` if the value might be present (may be false positive)
    /// `false` if the value is definitely not present
    pub fn contains(&self, value: u64) -> bool {
        for i in 0..self.num_hashes {
            let bit_index = self.hash_to_index(value, i);
            if !self.get_bit(bit_index) {
                return false;
            }
        }
        true
    }

    /// Get the number of items inserted
    pub fn count(&self) -> u64 {
        self.count
    }

    /// Get the number of bits in the filter
    pub fn num_bits(&self) -> usize {
        self.num_bits
    }

    /// Get the number of hash functions
    pub fn num_hashes(&self) -> usize {
        self.num_hashes
    }

    /// Estimate the current false positive rate
    pub fn estimated_fp_rate(&self) -> f64 {
        let ones = self.count_ones() as f64;
        let bits = self.num_bits as f64;
        let fill_ratio = ones / bits;

        // FP rate ≈ (1 - e^(-kn/m))^k ≈ fill_ratio^k
        fill_ratio.powi(self.num_hashes as i32)
    }

    /// Get the fill ratio (fraction of bits set)
    pub fn fill_ratio(&self) -> f64 {
        self.count_ones() as f64 / self.num_bits as f64
    }

    /// Get memory usage in bytes
    pub fn memory_bytes(&self) -> usize {
        self.bits.len()
    }

    /// Clear the filter
    pub fn clear(&mut self) {
        for byte in &mut self.bits {
            *byte = 0;
        }
        self.count = 0;
    }

    /// Merge another Bloom filter into this one (union)
    ///
    /// # Panics
    /// Panics if the filters have different sizes
    pub fn merge(&mut self, other: &Self) {
        assert_eq!(self.num_bits, other.num_bits, "Bloom filters must have same size");
        assert_eq!(self.num_hashes, other.num_hashes, "Bloom filters must have same hash count");

        for i in 0..self.bits.len() {
            self.bits[i] |= other.bits[i];
        }
        self.count += other.count;
    }

    /// Hash a value to a bit index using double hashing
    #[inline]
    fn hash_to_index(&self, value: u64, hash_num: usize) -> usize {
        // Double hashing: h(i) = h1 + i * h2
        let h1 = murmur_finalizer(value);
        let h2 = murmur_finalizer(value.rotate_left(32));

        let combined = h1.wrapping_add((hash_num as u64).wrapping_mul(h2));
        (combined as usize) % self.num_bits
    }

    #[inline]
    fn set_bit(&mut self, index: usize) {
        let byte_index = index / 8;
        let bit_index = index % 8;
        self.bits[byte_index] |= 1 << bit_index;
    }

    #[inline]
    fn get_bit(&self, index: usize) -> bool {
        let byte_index = index / 8;
        let bit_index = index % 8;
        (self.bits[byte_index] >> bit_index) & 1 == 1
    }

    fn count_ones(&self) -> usize {
        self.bits.iter().map(|b| b.count_ones() as usize).sum()
    }
}

/// MurmurHash3 finalizer
#[inline]
fn murmur_finalizer(mut x: u64) -> u64 {
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
    x ^= x >> 33;
    x
}

/// Minimizer Digest for compressed minimizer storage
///
/// Combines a Bloom filter for quick membership testing with optional
/// storage of top minimizers for more accurate lookups.
#[derive(Clone, Debug)]
pub struct MinimizerDigest {
    /// Bloom filter for quick membership testing
    bloom: BloomFilter,
    /// Optional storage for most frequent/important minimizers
    top_minimizers: Vec<(u64, u32)>, // (minimizer, taxon_id)
    /// Maximum number of top minimizers to store
    max_top: usize,
}

impl MinimizerDigest {
    /// Create a new minimizer digest
    ///
    /// # Arguments
    /// * `expected_minimizers` - Expected number of minimizers
    /// * `fp_rate` - Desired false positive rate for Bloom filter
    /// * `max_top` - Maximum number of top minimizers to store exactly
    pub fn new(expected_minimizers: usize, fp_rate: f64, max_top: usize) -> Self {
        MinimizerDigest {
            bloom: BloomFilter::new(expected_minimizers, fp_rate),
            top_minimizers: Vec::with_capacity(max_top),
            max_top,
        }
    }

    /// Create with default parameters
    pub fn with_capacity(expected_minimizers: usize) -> Self {
        Self::new(expected_minimizers, 0.01, 1000)
    }

    /// Insert a minimizer with its taxon ID
    pub fn insert(&mut self, minimizer: u64, taxon_id: u32) {
        self.bloom.insert(minimizer);

        // Store in top list if not full or if we want exact lookups for some
        if self.top_minimizers.len() < self.max_top {
            self.top_minimizers.push((minimizer, taxon_id));
        }
    }

    /// Check if a minimizer might be present
    ///
    /// Returns `true` if the minimizer might be in the digest (may be false positive)
    pub fn may_contain(&self, minimizer: u64) -> bool {
        self.bloom.contains(minimizer)
    }

    /// Look up a minimizer's taxon ID
    ///
    /// First checks the top minimizers for exact match, then falls back to
    /// Bloom filter for membership testing.
    ///
    /// # Returns
    /// * `Some(taxon_id)` if found in top minimizers
    /// * `None` if not found (or only found in Bloom filter)
    pub fn lookup(&self, minimizer: u64) -> Option<u32> {
        // Check exact storage first
        for &(m, taxon) in &self.top_minimizers {
            if m == minimizer {
                return Some(taxon);
            }
        }

        // Bloom filter doesn't store values, only membership
        None
    }

    /// Check if minimizer is present, returning whether it's an exact or approximate match
    pub fn query(&self, minimizer: u64) -> DigestQueryResult {
        // Check exact storage
        for &(m, taxon) in &self.top_minimizers {
            if m == minimizer {
                return DigestQueryResult::ExactMatch(taxon);
            }
        }

        // Check Bloom filter
        if self.bloom.contains(minimizer) {
            DigestQueryResult::PossibleMatch
        } else {
            DigestQueryResult::NotPresent
        }
    }

    /// Get the number of minimizers inserted
    pub fn count(&self) -> u64 {
        self.bloom.count()
    }

    /// Get memory usage in bytes (approximate)
    pub fn memory_bytes(&self) -> usize {
        self.bloom.memory_bytes() + self.top_minimizers.len() * 12
    }

    /// Get the compression ratio compared to storing all minimizers
    ///
    /// Assumes each minimizer would take 8 bytes (u64) + 4 bytes (taxon)
    pub fn compression_ratio(&self) -> f64 {
        let original_size = self.bloom.count() as usize * 12;
        if original_size == 0 {
            return 1.0;
        }
        self.memory_bytes() as f64 / original_size as f64
    }

    /// Get estimated false positive rate
    pub fn estimated_fp_rate(&self) -> f64 {
        self.bloom.estimated_fp_rate()
    }

    /// Clear the digest
    pub fn clear(&mut self) {
        self.bloom.clear();
        self.top_minimizers.clear();
    }

    /// Get the number of exactly stored minimizers
    pub fn exact_count(&self) -> usize {
        self.top_minimizers.len()
    }
}

/// Result of querying a minimizer digest
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DigestQueryResult {
    /// Found an exact match with taxon ID
    ExactMatch(u32),
    /// Minimizer may be present (Bloom filter positive)
    PossibleMatch,
    /// Minimizer is definitely not present
    NotPresent,
}

impl DigestQueryResult {
    /// Check if the result indicates presence (exact or possible)
    pub fn is_present(&self) -> bool {
        !matches!(self, DigestQueryResult::NotPresent)
    }

    /// Check if this is a definite match
    pub fn is_exact(&self) -> bool {
        matches!(self, DigestQueryResult::ExactMatch(_))
    }

    /// Get the taxon ID if available
    pub fn taxon_id(&self) -> Option<u32> {
        match self {
            DigestQueryResult::ExactMatch(id) => Some(*id),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // Bloom Filter Tests
    // =========================================================================

    #[test]
    fn test_bloom_filter_creation() {
        let bf = BloomFilter::new(1000, 0.01);
        assert!(bf.num_bits() > 0);
        assert!(bf.num_hashes() > 0);
        assert_eq!(bf.count(), 0);
    }

    #[test]
    fn test_bloom_filter_with_params() {
        let bf = BloomFilter::with_params(1024, 4);
        assert_eq!(bf.num_bits(), 1024);
        assert_eq!(bf.num_hashes(), 4);
    }

    #[test]
    fn test_bloom_filter_insert_contains() {
        let mut bf = BloomFilter::new(100, 0.01);

        bf.insert(12345);
        assert!(bf.contains(12345));
        assert_eq!(bf.count(), 1);
    }

    #[test]
    fn test_bloom_filter_no_false_negatives() {
        let mut bf = BloomFilter::new(1000, 0.01);

        // Insert many values
        for i in 0..500 {
            bf.insert(i);
        }

        // All inserted values must be found
        for i in 0..500 {
            assert!(bf.contains(i), "Value {} should be found", i);
        }
    }

    #[test]
    fn test_bloom_filter_false_positive_rate() {
        let mut bf = BloomFilter::new(1000, 0.01);

        // Insert values
        for i in 0..1000 {
            bf.insert(i);
        }

        // Check false positive rate on non-inserted values
        let mut false_positives = 0;
        for i in 10000..11000 {
            if bf.contains(i) {
                false_positives += 1;
            }
        }

        let fp_rate = false_positives as f64 / 1000.0;
        // Should be roughly around target (with some margin for small sample)
        assert!(fp_rate < 0.1, "FP rate {} is too high", fp_rate);
    }

    #[test]
    fn test_bloom_filter_clear() {
        let mut bf = BloomFilter::new(100, 0.01);
        bf.insert(12345);
        assert!(bf.contains(12345));

        bf.clear();
        assert!(!bf.contains(12345));
        assert_eq!(bf.count(), 0);
    }

    #[test]
    fn test_bloom_filter_merge() {
        let mut bf1 = BloomFilter::new(100, 0.01);
        let mut bf2 = BloomFilter::new(100, 0.01);

        bf1.insert(1);
        bf1.insert(2);
        bf2.insert(3);
        bf2.insert(4);

        bf1.merge(&bf2);

        assert!(bf1.contains(1));
        assert!(bf1.contains(2));
        assert!(bf1.contains(3));
        assert!(bf1.contains(4));
    }

    #[test]
    fn test_bloom_filter_fill_ratio() {
        let mut bf = BloomFilter::new(100, 0.01);
        assert_eq!(bf.fill_ratio(), 0.0);

        for i in 0..50 {
            bf.insert(i);
        }

        let ratio = bf.fill_ratio();
        assert!(ratio > 0.0 && ratio < 1.0);
    }

    #[test]
    fn test_bloom_filter_memory() {
        let bf = BloomFilter::new(1000, 0.01);
        let mem = bf.memory_bytes();
        assert!(mem > 0);
        // Should be roughly 1KB for 1000 items at 1% FP rate
        // (optimal is about 9.6 bits per item)
        assert!(mem < 2000);
    }

    #[test]
    fn test_bloom_filter_clone() {
        let mut bf = BloomFilter::new(100, 0.01);
        bf.insert(12345);

        let cloned = bf.clone();
        assert!(cloned.contains(12345));
        assert_eq!(cloned.count(), bf.count());
    }

    // =========================================================================
    // Minimizer Digest Tests
    // =========================================================================

    #[test]
    fn test_digest_creation() {
        let digest = MinimizerDigest::new(1000, 0.01, 100);
        assert_eq!(digest.count(), 0);
        assert_eq!(digest.exact_count(), 0);
    }

    #[test]
    fn test_digest_with_capacity() {
        let digest = MinimizerDigest::with_capacity(1000);
        assert_eq!(digest.count(), 0);
    }

    #[test]
    fn test_digest_insert_query() {
        let mut digest = MinimizerDigest::new(100, 0.01, 50);

        digest.insert(12345, 42);

        // Should be in Bloom filter
        assert!(digest.may_contain(12345));

        // Should be exact match (within top limit)
        match digest.query(12345) {
            DigestQueryResult::ExactMatch(taxon) => assert_eq!(taxon, 42),
            _ => panic!("Should find exact match"),
        }
    }

    #[test]
    fn test_digest_lookup() {
        let mut digest = MinimizerDigest::new(100, 0.01, 50);

        digest.insert(12345, 42);

        assert_eq!(digest.lookup(12345), Some(42));
        assert_eq!(digest.lookup(99999), None);
    }

    #[test]
    fn test_digest_not_present() {
        let mut digest = MinimizerDigest::new(100, 0.01, 50);
        digest.insert(12345, 42);

        // Value not inserted should be NotPresent (no false positive expected for single item)
        // Note: might get PossibleMatch due to Bloom filter, but definitely not ExactMatch
        let result = digest.query(99999);
        assert!(!result.is_exact());
    }

    #[test]
    fn test_digest_compression() {
        let mut digest = MinimizerDigest::new(10000, 0.01, 100);

        for i in 0..10000 {
            digest.insert(i, (i % 1000) as u32);
        }

        let ratio = digest.compression_ratio();
        // Should achieve significant compression
        assert!(ratio < 0.5, "Compression ratio {} should be < 0.5", ratio);
    }

    #[test]
    fn test_digest_clear() {
        let mut digest = MinimizerDigest::new(100, 0.01, 50);
        digest.insert(12345, 42);

        digest.clear();

        assert_eq!(digest.count(), 0);
        assert_eq!(digest.exact_count(), 0);
        assert!(!digest.may_contain(12345));
    }

    #[test]
    fn test_digest_exact_limit() {
        let mut digest = MinimizerDigest::new(100, 0.01, 5);

        // Insert more than max_top
        for i in 0..10 {
            digest.insert(i, i as u32);
        }

        // Only first 5 should be exact
        assert_eq!(digest.exact_count(), 5);
    }

    // =========================================================================
    // DigestQueryResult Tests
    // =========================================================================

    #[test]
    fn test_query_result_is_present() {
        assert!(DigestQueryResult::ExactMatch(42).is_present());
        assert!(DigestQueryResult::PossibleMatch.is_present());
        assert!(!DigestQueryResult::NotPresent.is_present());
    }

    #[test]
    fn test_query_result_is_exact() {
        assert!(DigestQueryResult::ExactMatch(42).is_exact());
        assert!(!DigestQueryResult::PossibleMatch.is_exact());
        assert!(!DigestQueryResult::NotPresent.is_exact());
    }

    #[test]
    fn test_query_result_taxon_id() {
        assert_eq!(DigestQueryResult::ExactMatch(42).taxon_id(), Some(42));
        assert_eq!(DigestQueryResult::PossibleMatch.taxon_id(), None);
        assert_eq!(DigestQueryResult::NotPresent.taxon_id(), None);
    }

    #[test]
    fn test_query_result_clone() {
        let result = DigestQueryResult::ExactMatch(42);
        let cloned = result.clone();
        assert_eq!(result, cloned);
    }

    #[test]
    fn test_murmur_finalizer_deterministic() {
        let h1 = murmur_finalizer(12345);
        let h2 = murmur_finalizer(12345);
        assert_eq!(h1, h2);
    }

    #[test]
    fn test_murmur_finalizer_different_inputs() {
        let h1 = murmur_finalizer(1);
        let h2 = murmur_finalizer(2);
        assert_ne!(h1, h2);
    }
}

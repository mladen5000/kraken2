/// Minimizer scanner for k-mer extraction
///
/// Translated from mmscanner.cc/mmscanner.h
///
/// This module extracts minimizers (the lexicographically smallest k-mer in a window)
/// from sequences for efficient indexing.

use crate::kraken2_data::IndexOptions;

/// Minimizer scanner that extracts minimizers from sequences
///
/// Supports both batch-style (scan()) and iterator-style (load_sequence/next_minimizer) interfaces.
pub struct MinimizerScanner {
    k: usize,
    l: usize,
    #[allow(dead_code)]
    spaced_seed_mask: u64,
    #[allow(dead_code)]
    toggle_mask: u64,
    dna_db: bool,
    // Iterator state
    sequence: Vec<u8>,
    position: usize,
    last_ambig: bool,
    last_minimizer: u64,
}

impl MinimizerScanner {
    /// Create a new minimizer scanner from index options
    pub fn new(opts: &IndexOptions) -> Self {
        MinimizerScanner {
            k: opts.k,
            l: opts.l,
            spaced_seed_mask: opts.spaced_seed_mask,
            toggle_mask: opts.toggle_mask,
            dna_db: opts.dna_db,
            sequence: Vec::new(),
            position: 0,
            last_ambig: false,
            last_minimizer: 0,
        }
    }

    /// Load a sequence for iterator-style scanning
    pub fn load_sequence(&mut self, seq: &str) {
        self.sequence = seq.as_bytes().to_vec();
        self.position = 0;
        self.last_ambig = false;
        self.last_minimizer = 0;
    }

    /// Get the next minimizer from the loaded sequence
    /// Returns None when no more minimizers are available
    pub fn next_minimizer(&mut self) -> Option<u64> {
        if self.sequence.is_empty() {
            return None;
        }

        let window_size = if self.l > 0 { self.l } else { self.k };

        while self.position + window_size <= self.sequence.len() {
            // Check for ambiguous bases in this window
            let window_end = self.position + window_size;
            let window = &self.sequence[self.position..window_end];

            // Check for ambiguous (N) bases
            self.last_ambig = window.iter().any(|&b| {
                !matches!(b, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't')
            });

            if self.last_ambig {
                self.position += 1;
                continue;
            }

            // Compute minimizer for this window
            if let Some(min_hash) = self.get_minimizer_hash(&self.sequence, self.position, window_size) {
                self.last_minimizer = min_hash;
                self.position += 1;
                return Some(min_hash);
            }

            self.position += 1;
        }

        None
    }

    /// Check if the last position was ambiguous
    pub fn is_ambiguous(&self) -> bool {
        self.last_ambig
    }

    /// Get the last minimizer computed
    pub fn last_minimizer(&self) -> u64 {
        self.last_minimizer
    }

    /// Extract all minimizers from a sequence
    /// Returns a vector of (minimizer_value, position) tuples
    pub fn scan(&self, seq: &str) -> Vec<(u64, usize)> {
        let mut minimizers = Vec::new();

        if seq.len() < self.k {
            return minimizers;
        }

        // Window size for minimizer selection (l parameter)
        let window_size = if self.l > 0 { self.l } else { self.k };

        // Sliding window approach: for each position, find the minimizer in the window
        for i in 0..=(seq.len().saturating_sub(window_size)) {
            // Get the window of k-mers
            let window_end = std::cmp::min(i + window_size, seq.len());

            if let Some(min_hash) = self.get_minimizer_hash(seq.as_bytes(), i, window_end - i) {
                minimizers.push((min_hash, i));
            }
        }

        minimizers
    }

    /// Get the lexicographically smallest k-mer hash in a window
    #[allow(dead_code)]
    fn get_minimizer_hash(&self, seq: &[u8], start: usize, window_size: usize) -> Option<u64> {
        if start + window_size > seq.len() {
            return None;
        }

        let mut min_hash = u64::MAX;
        let mut _min_pos = 0;

        for i in 0..window_size {
            let kmer = &seq[start + i..start + i + self.k];
            if let Some(hash) = self.hash_kmer(kmer) {
                if hash < min_hash {
                    min_hash = hash;
                    _min_pos = i;
                }
            }
        }

        if min_hash == u64::MAX {
            None
        } else {
            Some(min_hash)
        }
    }

    /// Hash a k-mer to a 64-bit value
    #[allow(dead_code)]
    fn hash_kmer(&self, kmer: &[u8]) -> Option<u64> {
        if kmer.len() != self.k {
            return None;
        }

        let mut hash: u64 = 0;

        if self.dna_db {
            for &byte in kmer {
                hash <<= 2;
                match byte {
                    b'A' | b'a' => {},
                    b'C' | b'c' => hash |= 1,
                    b'G' | b'g' => hash |= 2,
                    b'T' | b't' => hash |= 3,
                    _ => return None, // Ambiguous base
                }
            }
        } else {
            // Protein sequences - simplified hashing
            for &byte in kmer {
                hash = hash.wrapping_mul(31).wrapping_add(byte as u64);
            }
        }

        Some(hash)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimizer_scanner_creation() {
        let opts = IndexOptions::default();
        let scanner = MinimizerScanner::new(&opts);
        assert_eq!(scanner.k, 35);
        assert_eq!(scanner.l, 31);
    }

    #[test]
    fn test_kmer_hash_dna() {
        let opts = IndexOptions::default();
        let scanner = MinimizerScanner::new(&opts);
        let kmer = b"ATCG";
        let hash = scanner.hash_kmer(kmer);
        // ATCG in 2-bit form: A=0, T=3, C=1, G=2 = 00111011 in binary
        assert!(hash.is_some());
    }
}

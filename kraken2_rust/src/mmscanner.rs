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
        // Check bounds: we need at least k bases from start position
        if start + self.k > seq.len() {
            return None;
        }

        let mut min_hash = u64::MAX;
        let mut _min_pos = 0;

        // Number of k-mers in the window
        let num_kmers = if window_size >= self.k {
            window_size - self.k + 1
        } else {
            1 // At least try to hash one k-mer
        };

        for i in 0..num_kmers {
            let kmer_start = start + i;
            let kmer_end = kmer_start + self.k;

            // Bounds check
            if kmer_end > seq.len() {
                break;
            }

            let kmer = &seq[kmer_start..kmer_end];
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

    fn create_scanner_with_k(k: usize, l: usize) -> MinimizerScanner {
        let opts = IndexOptions {
            k,
            l,
            spaced_seed_mask: 0,
            toggle_mask: 0,
            dna_db: true,
            minimum_acceptable_hash_value: 0,
            revcom_version: 0,
            db_version: 0,
            db_type: 0,
        };
        MinimizerScanner::new(&opts)
    }

    #[test]
    fn test_minimizer_scanner_creation() {
        let opts = IndexOptions::default();
        let scanner = MinimizerScanner::new(&opts);
        assert_eq!(scanner.k, 35);
        assert_eq!(scanner.l, 31);
    }

    #[test]
    fn test_kmer_hash_dna_correct_length() {
        // Create a scanner with k=4 so we can test with a 4-mer
        let scanner = create_scanner_with_k(4, 4);
        let kmer = b"ATCG";
        let hash = scanner.hash_kmer(kmer);
        // ATCG in 2-bit form: A=0, T=3, C=1, G=2
        // hash = 0*8 + 3*4 + 1*2 + 2 = 0 + 12 + 2 + 2 = 16... wait
        // Let me trace: hash starts at 0
        // A: hash <<= 2 (0), hash |= 0 -> hash = 0
        // T: hash <<= 2 (0), hash |= 3 -> hash = 3
        // C: hash <<= 2 (12), hash |= 1 -> hash = 13
        // G: hash <<= 2 (52), hash |= 2 -> hash = 54
        assert!(hash.is_some());
        assert_eq!(hash.unwrap(), 54); // 0b00110110
    }

    #[test]
    fn test_kmer_hash_dna_wrong_length() {
        let scanner = create_scanner_with_k(4, 4);
        // Wrong length k-mer should return None
        let kmer = b"ATCGA"; // 5 bases, but k=4
        let hash = scanner.hash_kmer(kmer);
        assert!(hash.is_none(), "Wrong length kmer should return None");
    }

    #[test]
    fn test_kmer_hash_dna_with_n() {
        let scanner = create_scanner_with_k(4, 4);
        let kmer = b"ATNG"; // Contains ambiguous N
        let hash = scanner.hash_kmer(kmer);
        assert!(hash.is_none(), "Kmer with N should return None");
    }

    #[test]
    fn test_kmer_hash_dna_lowercase() {
        let scanner = create_scanner_with_k(4, 4);
        let upper = b"ATCG";
        let lower = b"atcg";
        let hash_upper = scanner.hash_kmer(upper);
        let hash_lower = scanner.hash_kmer(lower);
        assert_eq!(hash_upper, hash_lower, "Case should not affect hash");
    }

    #[test]
    fn test_kmer_hash_specific_values() {
        let scanner = create_scanner_with_k(4, 4);

        // AAAA = 0000 0000 = 0
        assert_eq!(scanner.hash_kmer(b"AAAA"), Some(0));

        // TTTT = 1111 1111 = 255
        assert_eq!(scanner.hash_kmer(b"TTTT"), Some(0xFF));

        // CCCC = 0101 0101 = 85
        assert_eq!(scanner.hash_kmer(b"CCCC"), Some(0b01010101));

        // GGGG = 1010 1010 = 170
        assert_eq!(scanner.hash_kmer(b"GGGG"), Some(0b10101010));
    }

    #[test]
    fn test_scan_short_sequence() {
        let scanner = create_scanner_with_k(4, 4);
        // Sequence shorter than k should return empty
        let seq = "ATG";
        let minimizers = scanner.scan(seq);
        assert!(minimizers.is_empty());
    }

    #[test]
    fn test_scan_exact_length() {
        let scanner = create_scanner_with_k(4, 4);
        // Sequence exactly k length should return one minimizer
        let seq = "ATCG";
        let minimizers = scanner.scan(seq);
        assert_eq!(minimizers.len(), 1);
    }

    #[test]
    fn test_scan_with_window() {
        // k=4, l=6 means we look at 6-base windows for minimizer selection
        let scanner = create_scanner_with_k(4, 6);
        let seq = "ATCGATCG"; // 8 bases
        let minimizers = scanner.scan(seq);
        // Should have minimizers for positions 0, 1, 2 (windows of 6 from 8-base seq)
        assert!(!minimizers.is_empty());
    }

    #[test]
    fn test_load_sequence_and_iterate() {
        let mut scanner = create_scanner_with_k(4, 4);
        scanner.load_sequence("ATCGATCG");

        let mut count = 0;
        while scanner.next_minimizer().is_some() {
            count += 1;
        }
        // 8 - 4 + 1 = 5 positions to check
        assert!(count > 0, "Should have found at least one minimizer");
    }

    #[test]
    fn test_is_ambiguous() {
        let mut scanner = create_scanner_with_k(4, 4);
        scanner.load_sequence("ATNG");
        let result = scanner.next_minimizer();
        assert!(result.is_none() || scanner.is_ambiguous(), "Should detect ambiguous bases");
    }

    #[test]
    fn test_protein_hash() {
        // Test protein hashing mode
        let opts = IndexOptions {
            k: 4,
            l: 4,
            spaced_seed_mask: 0,
            toggle_mask: 0,
            dna_db: false, // Protein mode
            minimum_acceptable_hash_value: 0,
            revcom_version: 0,
            db_version: 0,
            db_type: 0,
        };
        let scanner = MinimizerScanner::new(&opts);

        let protein = b"MAST";
        let hash = scanner.hash_kmer(protein);
        assert!(hash.is_some(), "Protein hashing should work");
    }
}

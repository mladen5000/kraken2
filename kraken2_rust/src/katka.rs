/// KATKA Kernels - MEM-based Classification Enhancement
///
/// This module implements Maximal Exact Matches (MEMs) for improved strain-level
/// classification accuracy. Based on the KATKA paper (arXiv:2402.06935).
///
/// Key concepts:
/// - MEMs are maximal exact matches between query and reference sequences
/// - Longer MEMs provide stronger evidence for classification
/// - MEM-based scoring complements minimizer-based classification
///
/// Benefits:
/// - Better strain-level differentiation
/// - Improved accuracy for closely related species
/// - More robust to sequencing errors at match boundaries

use std::collections::HashMap;
use crate::kraken2_data::TaxId;

/// Suffix array for efficient MEM finding
///
/// Uses a simplified suffix array construction suitable for
/// moderate-sized reference sequences.
#[derive(Clone, Debug)]
pub struct SuffixArray {
    /// The text (concatenated reference sequences with separators)
    text: Vec<u8>,
    /// Suffix array (sorted indices into text)
    sa: Vec<usize>,
    /// Inverse suffix array (position -> rank)
    isa: Vec<usize>,
    /// LCP array (longest common prefix between adjacent suffixes)
    lcp: Vec<usize>,
    /// Taxon ID for each position in text
    taxon_map: Vec<TaxId>,
}

impl SuffixArray {
    /// Build a suffix array from reference sequences
    ///
    /// # Arguments
    /// * `sequences` - Vec of (sequence, taxon_id) pairs
    /// * `separator` - Byte used to separate sequences (typically '$' or 0)
    pub fn build(sequences: &[(Vec<u8>, TaxId)], separator: u8) -> Self {
        // Concatenate all sequences with separators
        let mut text = Vec::new();
        let mut taxon_map = Vec::new();

        for (seq, taxon) in sequences {
            for &b in seq {
                text.push(b);
                taxon_map.push(*taxon);
            }
            text.push(separator);
            taxon_map.push(0); // Separator has no taxon
        }

        let n = text.len();

        // Build suffix array using simple doubling algorithm
        // For production, consider using libdivsufsort or similar
        let sa = Self::build_suffix_array(&text);

        // Build inverse suffix array
        let mut isa = vec![0; n];
        for (rank, &pos) in sa.iter().enumerate() {
            isa[pos] = rank;
        }

        // Build LCP array using Kasai's algorithm
        let lcp = Self::build_lcp_array(&text, &sa, &isa);

        SuffixArray {
            text,
            sa,
            isa,
            lcp,
            taxon_map,
        }
    }

    /// Build suffix array using prefix doubling (O(n log^2 n))
    fn build_suffix_array(text: &[u8]) -> Vec<usize> {
        let n = text.len();
        if n == 0 {
            return Vec::new();
        }

        // Initialize ranks from characters
        let mut rank: Vec<i64> = text.iter().map(|&c| c as i64).collect();
        let mut sa: Vec<usize> = (0..n).collect();
        let mut tmp = vec![0i64; n];

        let mut k = 1usize;
        while k < n {
            // Sort by (rank[i], rank[i+k])
            sa.sort_by(|&a, &b| {
                let ra = rank[a];
                let rb = rank[b];
                if ra != rb {
                    return ra.cmp(&rb);
                }
                let ra2 = if a + k < n { rank[a + k] } else { -1 };
                let rb2 = if b + k < n { rank[b + k] } else { -1 };
                ra2.cmp(&rb2)
            });

            // Compute new ranks
            tmp[sa[0]] = 0;
            for i in 1..n {
                let prev = sa[i - 1];
                let curr = sa[i];

                let same = rank[prev] == rank[curr] && {
                    let r1 = if prev + k < n { rank[prev + k] } else { -1 };
                    let r2 = if curr + k < n { rank[curr + k] } else { -1 };
                    r1 == r2
                };

                tmp[curr] = tmp[prev] + if same { 0 } else { 1 };
            }

            std::mem::swap(&mut rank, &mut tmp);

            if rank[sa[n - 1]] == (n - 1) as i64 {
                break; // All ranks are unique
            }

            k *= 2;
        }

        sa
    }

    /// Build LCP array using Kasai's algorithm (O(n))
    fn build_lcp_array(text: &[u8], sa: &[usize], isa: &[usize]) -> Vec<usize> {
        let n = text.len();
        if n == 0 {
            return Vec::new();
        }

        let mut lcp = vec![0usize; n];
        let mut k = 0usize;

        for i in 0..n {
            if isa[i] == 0 {
                k = 0;
                continue;
            }

            let j = sa[isa[i] - 1];
            while i + k < n && j + k < n && text[i + k] == text[j + k] {
                k += 1;
            }

            lcp[isa[i]] = k;
            if k > 0 {
                k -= 1;
            }
        }

        lcp
    }

    /// Find all MEMs (Maximal Exact Matches) of a query against the reference
    ///
    /// # Arguments
    /// * `query` - The query sequence
    /// * `min_len` - Minimum MEM length to report
    ///
    /// # Returns
    /// Vector of MEMs as (query_start, ref_start, length, taxon_id)
    pub fn find_mems(&self, query: &[u8], min_len: usize) -> Vec<MEM> {
        let mut mems = Vec::new();

        if query.is_empty() || self.text.is_empty() {
            return mems;
        }

        // For each position in query, find maximal matches
        for q_start in 0..query.len() {
            // Binary search for the range of suffixes that match query[q_start..]
            let matches = self.find_matches_at_position(query, q_start, min_len);

            for (ref_pos, match_len, taxon) in matches {
                // Check if this is maximal (can't extend left)
                let is_maximal = q_start == 0
                    || ref_pos == 0
                    || query[q_start - 1] != self.text[ref_pos - 1];

                if is_maximal && match_len >= min_len {
                    mems.push(MEM {
                        query_start: q_start,
                        ref_start: ref_pos,
                        length: match_len,
                        taxon_id: taxon,
                    });
                }
            }
        }

        mems
    }

    /// Find matches at a specific query position using binary search
    fn find_matches_at_position(
        &self,
        query: &[u8],
        q_start: usize,
        min_len: usize,
    ) -> Vec<(usize, usize, TaxId)> {
        let mut results = Vec::new();
        let query_suffix = &query[q_start..];

        if query_suffix.len() < min_len {
            return results;
        }

        // Binary search for lower bound
        let mut lo = 0;
        let mut hi = self.sa.len();

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let ref_suffix = &self.text[self.sa[mid]..];

            if ref_suffix < query_suffix {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }

        let start = lo;

        // Binary search for upper bound (first suffix > query prefix of min_len)
        hi = self.sa.len();

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let ref_suffix = &self.text[self.sa[mid]..];

            // Compare with query_suffix prefix
            let cmp_len = std::cmp::min(query_suffix.len(), ref_suffix.len());
            let ref_prefix = &ref_suffix[..std::cmp::min(min_len, cmp_len)];
            let query_prefix = &query_suffix[..std::cmp::min(min_len, query_suffix.len())];

            if ref_prefix <= query_prefix {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }

        let end = lo;

        // Collect all matches in range
        for i in start..end {
            let ref_pos = self.sa[i];
            let ref_suffix = &self.text[ref_pos..];

            // Find actual match length
            let match_len = query_suffix
                .iter()
                .zip(ref_suffix.iter())
                .take_while(|(a, b)| a == b)
                .count();

            if match_len >= min_len {
                let taxon = self.taxon_map[ref_pos];
                if taxon != 0 {
                    results.push((ref_pos, match_len, taxon));
                }
            }
        }

        results
    }

    /// Get the length of the indexed text
    pub fn text_len(&self) -> usize {
        self.text.len()
    }

    /// Get the number of suffixes
    pub fn num_suffixes(&self) -> usize {
        self.sa.len()
    }
}

/// A Maximal Exact Match (MEM)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MEM {
    /// Start position in query
    pub query_start: usize,
    /// Start position in reference
    pub ref_start: usize,
    /// Length of the match
    pub length: usize,
    /// Taxon ID of the reference sequence
    pub taxon_id: TaxId,
}

impl MEM {
    /// Get the end position in query (exclusive)
    pub fn query_end(&self) -> usize {
        self.query_start + self.length
    }

    /// Get the end position in reference (exclusive)
    pub fn ref_end(&self) -> usize {
        self.ref_start + self.length
    }
}

/// MEM-based classification scorer
///
/// Computes classification scores based on MEM coverage and length distribution
#[derive(Clone, Debug)]
pub struct MEMScorer {
    /// Minimum MEM length to consider
    min_mem_length: usize,
    /// Weight factor for longer MEMs
    length_weight: f64,
    /// Minimum coverage fraction for confident classification
    min_coverage: f64,
}

impl Default for MEMScorer {
    fn default() -> Self {
        MEMScorer {
            min_mem_length: 31,  // Same as default minimizer window
            length_weight: 1.5,  // Longer MEMs weighted more heavily
            min_coverage: 0.5,   // At least 50% query coverage
        }
    }
}

impl MEMScorer {
    /// Create a new MEM scorer with custom parameters
    pub fn new(min_mem_length: usize, length_weight: f64, min_coverage: f64) -> Self {
        MEMScorer {
            min_mem_length,
            length_weight,
            min_coverage,
        }
    }

    /// Score MEMs for classification
    ///
    /// Returns a map of taxon_id -> score
    pub fn score_mems(&self, mems: &[MEM], query_len: usize) -> HashMap<TaxId, f64> {
        let mut scores: HashMap<TaxId, f64> = HashMap::new();

        if mems.is_empty() || query_len == 0 {
            return scores;
        }

        for mem in mems {
            if mem.length < self.min_mem_length {
                continue;
            }

            // Score = length^weight (longer MEMs count more)
            let score = (mem.length as f64).powf(self.length_weight);
            *scores.entry(mem.taxon_id).or_insert(0.0) += score;
        }

        // Normalize by query length
        let norm_factor = (query_len as f64).powf(self.length_weight);
        for score in scores.values_mut() {
            *score /= norm_factor;
        }

        scores
    }

    /// Compute coverage of query by MEMs for a specific taxon
    pub fn compute_coverage(&self, mems: &[MEM], query_len: usize, taxon: TaxId) -> f64 {
        if query_len == 0 {
            return 0.0;
        }

        // Collect intervals for this taxon
        let mut intervals: Vec<(usize, usize)> = mems
            .iter()
            .filter(|m| m.taxon_id == taxon && m.length >= self.min_mem_length)
            .map(|m| (m.query_start, m.query_end()))
            .collect();

        if intervals.is_empty() {
            return 0.0;
        }

        // Sort by start position
        intervals.sort_by_key(|&(start, _)| start);

        // Merge overlapping intervals
        let mut merged = Vec::new();
        let mut current = intervals[0];

        for &(start, end) in intervals.iter().skip(1) {
            if start <= current.1 {
                current.1 = std::cmp::max(current.1, end);
            } else {
                merged.push(current);
                current = (start, end);
            }
        }
        merged.push(current);

        // Sum covered bases
        let covered: usize = merged.iter().map(|(s, e)| e - s).sum();

        covered as f64 / query_len as f64
    }

    /// Get the best classification based on MEM scores
    ///
    /// Returns (taxon_id, score, coverage) or None if no confident classification
    pub fn classify(&self, mems: &[MEM], query_len: usize) -> Option<(TaxId, f64, f64)> {
        let scores = self.score_mems(mems, query_len);

        if scores.is_empty() {
            return None;
        }

        // Find taxon with highest score
        let (best_taxon, best_score) = scores
            .iter()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))?;

        let coverage = self.compute_coverage(mems, query_len, *best_taxon);

        if coverage >= self.min_coverage {
            Some((*best_taxon, *best_score, coverage))
        } else {
            None
        }
    }

    /// Get minimum MEM length
    pub fn min_mem_length(&self) -> usize {
        self.min_mem_length
    }

    /// Get minimum coverage threshold
    pub fn min_coverage(&self) -> f64 {
        self.min_coverage
    }
}

/// KATKA-enhanced classification result
///
/// Combines minimizer-based and MEM-based classification
#[derive(Clone, Debug)]
pub struct KatkaClassification {
    /// Taxon ID from combined scoring
    pub taxon_id: TaxId,
    /// Minimizer-based score (0-1)
    pub minimizer_score: f64,
    /// MEM-based score (0-1)
    pub mem_score: f64,
    /// Combined confidence score
    pub confidence: f64,
    /// Query coverage by MEMs
    pub mem_coverage: f64,
    /// Number of MEMs found
    pub num_mems: usize,
    /// Average MEM length
    pub avg_mem_length: f64,
}

impl KatkaClassification {
    /// Create a classification result
    pub fn new(
        taxon_id: TaxId,
        minimizer_score: f64,
        mem_score: f64,
        mem_coverage: f64,
        mems: &[MEM],
    ) -> Self {
        let num_mems = mems.iter().filter(|m| m.taxon_id == taxon_id).count();
        let total_len: usize = mems
            .iter()
            .filter(|m| m.taxon_id == taxon_id)
            .map(|m| m.length)
            .sum();
        let avg_mem_length = if num_mems > 0 {
            total_len as f64 / num_mems as f64
        } else {
            0.0
        };

        // Combined confidence: weighted average of scores
        let confidence = 0.6 * minimizer_score + 0.4 * mem_score;

        KatkaClassification {
            taxon_id,
            minimizer_score,
            mem_score,
            confidence,
            mem_coverage,
            num_mems,
            avg_mem_length,
        }
    }

    /// Check if classification is confident
    pub fn is_confident(&self, threshold: f64) -> bool {
        self.confidence >= threshold
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // Suffix Array Tests
    // =========================================================================

    #[test]
    fn test_suffix_array_empty() {
        let sequences: Vec<(Vec<u8>, TaxId)> = vec![];
        let sa = SuffixArray::build(&sequences, b'$');
        assert_eq!(sa.text_len(), 0);
        assert_eq!(sa.num_suffixes(), 0);
    }

    #[test]
    fn test_suffix_array_single_sequence() {
        let sequences = vec![(b"ACGT".to_vec(), 1)];
        let sa = SuffixArray::build(&sequences, b'$');

        // Text should be "ACGT$"
        assert_eq!(sa.text_len(), 5);
        assert_eq!(sa.num_suffixes(), 5);
    }

    #[test]
    fn test_suffix_array_multiple_sequences() {
        let sequences = vec![
            (b"ACGT".to_vec(), 1),
            (b"TGCA".to_vec(), 2),
        ];
        let sa = SuffixArray::build(&sequences, b'$');

        // Text should be "ACGT$TGCA$"
        assert_eq!(sa.text_len(), 10);
    }

    #[test]
    fn test_suffix_array_sorted() {
        let sequences = vec![(b"BANANA".to_vec(), 1)];
        let sa = SuffixArray::build(&sequences, b'$');

        // Verify suffix array is sorted
        for i in 1..sa.sa.len() {
            let prev_suffix = &sa.text[sa.sa[i - 1]..];
            let curr_suffix = &sa.text[sa.sa[i]..];
            assert!(prev_suffix <= curr_suffix, "Suffix array not sorted at position {}", i);
        }
    }

    // =========================================================================
    // MEM Finding Tests
    // =========================================================================

    #[test]
    fn test_find_mems_empty_query() {
        let sequences = vec![(b"ACGTACGT".to_vec(), 1)];
        let sa = SuffixArray::build(&sequences, b'$');

        let mems = sa.find_mems(b"", 3);
        assert!(mems.is_empty());
    }

    #[test]
    fn test_find_mems_no_match() {
        let sequences = vec![(b"AAAAAAA".to_vec(), 1)];
        let sa = SuffixArray::build(&sequences, b'$');

        // Query with no common substring >= min_len
        let mems = sa.find_mems(b"TTTTTTT", 3);
        assert!(mems.is_empty());
    }

    #[test]
    fn test_find_mems_exact_match() {
        let sequences = vec![(b"ACGTACGT".to_vec(), 42)];
        let sa = SuffixArray::build(&sequences, b'$');

        // Query that matches exactly
        let mems = sa.find_mems(b"ACGT", 3);
        assert!(!mems.is_empty());

        // Should find matches of length 4
        let long_mems: Vec<_> = mems.iter().filter(|m| m.length >= 4).collect();
        assert!(!long_mems.is_empty());
        assert_eq!(long_mems[0].taxon_id, 42);
    }

    #[test]
    fn test_find_mems_partial_match() {
        let sequences = vec![(b"ACGTNNNNACGT".to_vec(), 1)];
        let sa = SuffixArray::build(&sequences, b'$');

        // Query that partially matches
        let mems = sa.find_mems(b"ACGTTGCA", 3);

        // Should find ACGT match
        let has_acgt_match = mems.iter().any(|m| m.length >= 4);
        assert!(has_acgt_match, "Should find ACGT match");
    }

    #[test]
    fn test_mem_structure() {
        let mem = MEM {
            query_start: 10,
            ref_start: 100,
            length: 50,
            taxon_id: 42,
        };

        assert_eq!(mem.query_end(), 60);
        assert_eq!(mem.ref_end(), 150);
    }

    // =========================================================================
    // MEM Scorer Tests
    // =========================================================================

    #[test]
    fn test_mem_scorer_default() {
        let scorer = MEMScorer::default();
        assert_eq!(scorer.min_mem_length(), 31);
        assert_eq!(scorer.min_coverage(), 0.5);
    }

    #[test]
    fn test_mem_scorer_custom() {
        let scorer = MEMScorer::new(20, 2.0, 0.7);
        assert_eq!(scorer.min_mem_length(), 20);
        assert_eq!(scorer.min_coverage(), 0.7);
    }

    #[test]
    fn test_score_mems_empty() {
        let scorer = MEMScorer::new(10, 1.5, 0.5);
        let scores = scorer.score_mems(&[], 100);
        assert!(scores.is_empty());
    }

    #[test]
    fn test_score_mems_single() {
        let scorer = MEMScorer::new(10, 1.0, 0.5); // Linear weight for easier testing

        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 50, taxon_id: 1 },
        ];

        let scores = scorer.score_mems(&mems, 100);
        assert!(scores.contains_key(&1));
        assert!(scores[&1] > 0.0);
    }

    #[test]
    fn test_score_mems_multiple_taxa() {
        let scorer = MEMScorer::new(10, 1.0, 0.5);

        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 50, taxon_id: 1 },
            MEM { query_start: 50, ref_start: 100, length: 30, taxon_id: 2 },
        ];

        let scores = scorer.score_mems(&mems, 100);
        assert!(scores.contains_key(&1));
        assert!(scores.contains_key(&2));
        // Taxon 1 should have higher score (longer MEM)
        assert!(scores[&1] > scores[&2]);
    }

    #[test]
    fn test_compute_coverage_empty() {
        let scorer = MEMScorer::new(10, 1.5, 0.5);
        let coverage = scorer.compute_coverage(&[], 100, 1);
        assert_eq!(coverage, 0.0);
    }

    #[test]
    fn test_compute_coverage_full() {
        let scorer = MEMScorer::new(10, 1.5, 0.5);

        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 100, taxon_id: 1 },
        ];

        let coverage = scorer.compute_coverage(&mems, 100, 1);
        assert!((coverage - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_compute_coverage_partial() {
        let scorer = MEMScorer::new(10, 1.5, 0.5);

        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 50, taxon_id: 1 },
        ];

        let coverage = scorer.compute_coverage(&mems, 100, 1);
        assert!((coverage - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_compute_coverage_overlapping() {
        let scorer = MEMScorer::new(10, 1.5, 0.5);

        // Two overlapping MEMs
        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 60, taxon_id: 1 },
            MEM { query_start: 40, ref_start: 100, length: 60, taxon_id: 1 },
        ];

        let coverage = scorer.compute_coverage(&mems, 100, 1);
        // Should cover 0-100, so full coverage
        assert!((coverage - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_classify_no_mems() {
        let scorer = MEMScorer::new(10, 1.5, 0.5);
        let result = scorer.classify(&[], 100);
        assert!(result.is_none());
    }

    #[test]
    fn test_classify_low_coverage() {
        let scorer = MEMScorer::new(10, 1.5, 0.8); // High coverage requirement

        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 20, taxon_id: 1 },
        ];

        let result = scorer.classify(&mems, 100);
        assert!(result.is_none()); // Coverage too low
    }

    #[test]
    fn test_classify_success() {
        let scorer = MEMScorer::new(10, 1.5, 0.5);

        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 60, taxon_id: 1 },
        ];

        let result = scorer.classify(&mems, 100);
        assert!(result.is_some());

        let (taxon, score, coverage) = result.unwrap();
        assert_eq!(taxon, 1);
        assert!(score > 0.0);
        assert!((coverage - 0.6).abs() < 0.001);
    }

    // =========================================================================
    // KATKA Classification Tests
    // =========================================================================

    #[test]
    fn test_katka_classification_creation() {
        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 50, taxon_id: 1 },
            MEM { query_start: 50, ref_start: 100, length: 30, taxon_id: 1 },
        ];

        let result = KatkaClassification::new(1, 0.8, 0.7, 0.9, &mems);

        assert_eq!(result.taxon_id, 1);
        assert_eq!(result.minimizer_score, 0.8);
        assert_eq!(result.mem_score, 0.7);
        assert_eq!(result.mem_coverage, 0.9);
        assert_eq!(result.num_mems, 2);
        assert!((result.avg_mem_length - 40.0).abs() < 0.001);
    }

    #[test]
    fn test_katka_classification_confidence() {
        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 50, taxon_id: 1 },
        ];

        let result = KatkaClassification::new(1, 0.8, 0.7, 0.9, &mems);

        // Confidence = 0.6 * 0.8 + 0.4 * 0.7 = 0.48 + 0.28 = 0.76
        assert!((result.confidence - 0.76).abs() < 0.001);

        assert!(result.is_confident(0.7));
        assert!(!result.is_confident(0.8));
    }

    #[test]
    fn test_katka_classification_no_mems_for_taxon() {
        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 50, taxon_id: 2 },
        ];

        let result = KatkaClassification::new(1, 0.8, 0.0, 0.0, &mems);

        assert_eq!(result.num_mems, 0);
        assert_eq!(result.avg_mem_length, 0.0);
    }

    #[test]
    fn test_katka_classification_clone() {
        let mems = vec![
            MEM { query_start: 0, ref_start: 0, length: 50, taxon_id: 1 },
        ];

        let result = KatkaClassification::new(1, 0.8, 0.7, 0.9, &mems);
        let cloned = result.clone();

        assert_eq!(cloned.taxon_id, result.taxon_id);
        assert_eq!(cloned.confidence, result.confidence);
    }

    // =========================================================================
    // Integration Tests
    // =========================================================================

    #[test]
    fn test_full_pipeline() {
        // Build index from reference sequences
        let sequences = vec![
            (b"ACGTACGTACGTACGT".to_vec(), 1),
            (b"TGCATGCATGCATGCA".to_vec(), 2),
        ];
        let sa = SuffixArray::build(&sequences, b'$');

        // Query that matches taxon 1
        let query = b"ACGTACGT";
        let mems = sa.find_mems(query, 4);

        // Score MEMs
        let scorer = MEMScorer::new(4, 1.5, 0.3);
        let scores = scorer.score_mems(&mems, query.len());

        // Should have score for taxon 1
        assert!(scores.contains_key(&1));
        // Taxon 1 score should be higher than taxon 2 (if any)
        if scores.contains_key(&2) {
            assert!(scores[&1] >= scores[&2]);
        }
    }

    #[test]
    fn test_lcp_array_values() {
        // Simple test to verify LCP array correctness
        let sequences = vec![(b"ABAB".to_vec(), 1)];
        let sa = SuffixArray::build(&sequences, b'$');

        // LCP should be computed without panicking
        assert_eq!(sa.lcp.len(), sa.sa.len());
    }
}

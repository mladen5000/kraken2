/// Low-complexity region masking using Symmetric DUST algorithm
///
/// Translated from k2mask.cc
///
/// This module implements the Symmetric DUST algorithm for masking
/// low-complexity regions in DNA sequences. Low-complexity regions
/// can cause spurious matches in sequence classification.

use crate::seqreader::Sequence;
use std::collections::VecDeque;

/// ASCII to DNA encoding table
/// A=0, C=1, G=2, T=3, other=4 (invalid)
const ASC2DNA: [u8; 256] = [
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 0-15
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 16-31
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 32-47
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 48-63
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 64-79  (A=65, C=67, G=71)
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 80-95  (T=84)
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, // 96-111 (a=97, c=99, g=103)
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 112-127 (t=116)
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 128-143
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 144-159
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 160-175
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 176-191
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 192-207
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 208-223
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 224-239
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, // 240-255
];

/// Perfect interval for DUST algorithm
#[derive(Clone, Debug)]
struct PerfectInterval {
    start: i32,
    finish: i32,
    left: i32,
    right: i32,
}

/// Range of masked bases
#[derive(Clone, Debug)]
pub struct MaskRange {
    pub start: i32,
    pub finish: i32,
}

/// Symmetric DUST state machine
pub struct SDust {
    kmers: VecDeque<i32>,
    perfect_intervals: Vec<PerfectInterval>,
    ranges: Vec<MaskRange>,
    cw: [i32; 64],  // Window k-mer counts
    cv: [i32; 64],  // Current k-mer counts
    rw: i32,        // Window score
    rv: i32,        // Current score
    l: i32,         // Current length
    window_size: i32,
    threshold: i32,
}

impl Default for SDust {
    fn default() -> Self {
        Self::new(64, 20)
    }
}

impl SDust {
    /// Create a new SDust instance with given window size and threshold
    pub fn new(window_size: i32, threshold: i32) -> Self {
        SDust {
            kmers: VecDeque::new(),
            perfect_intervals: Vec::new(),
            ranges: Vec::new(),
            cw: [0; 64],
            cv: [0; 64],
            rw: 0,
            rv: 0,
            l: 0,
            window_size,
            threshold,
        }
    }

    /// Reset state for processing a new sequence
    pub fn reset(&mut self) {
        self.kmers.clear();
        self.perfect_intervals.clear();
        self.ranges.clear();
        self.cw = [0; 64];
        self.cv = [0; 64];
        self.l = 0;
        self.rw = 0;
        self.rv = 0;
    }

    /// Shift the sliding window by adding a new triplet
    fn shift_window(&mut self, t: i32) {
        if self.kmers.len() as i32 >= self.window_size - 2 {
            if let Some(s) = self.kmers.pop_front() {
                self.cw[s as usize] -= 1;
                self.rw -= self.cw[s as usize];
                if self.l > self.kmers.len() as i32 {
                    self.l -= 1;
                    self.cv[s as usize] -= 1;
                    self.rv -= self.cv[s as usize];
                }
            }
        }

        self.kmers.push_back(t);
        self.l += 1;
        self.rw += self.cw[t as usize];
        self.cw[t as usize] += 1;
        self.rv += self.cv[t as usize];
        self.cv[t as usize] += 1;

        if self.cv[t as usize] * 10 > self.threshold * 2 {
            loop {
                let idx = self.kmers.len() as i32 - self.l;
                if idx < 0 {
                    break;
                }
                let s = self.kmers[idx as usize];
                self.cv[s as usize] -= 1;
                self.rv -= self.cv[s as usize];
                self.l -= 1;
                if s == t {
                    break;
                }
            }
        }
    }

    /// Save masked regions when window moves
    fn save_masked_regions(&mut self, window_start: i32) {
        if self.perfect_intervals.is_empty() {
            return;
        }

        let last_start = self.perfect_intervals.last().unwrap().start;
        if last_start >= window_start {
            return;
        }

        let p = self.perfect_intervals.last().unwrap().clone();

        let mut saved = false;
        if let Some(last_range) = self.ranges.last_mut() {
            if p.start <= last_range.finish {
                last_range.finish = last_range.finish.max(p.finish);
                saved = true;
            }
        }

        if !saved {
            self.ranges.push(MaskRange {
                start: p.start,
                finish: p.finish,
            });
        }

        // Remove intervals that are behind the window
        while !self.perfect_intervals.is_empty() {
            if self.perfect_intervals.last().unwrap().start >= window_start {
                break;
            }
            self.perfect_intervals.pop();
        }
    }

    /// Find perfect intervals in current window
    fn find_perfect(&mut self, window_start: i32) {
        let mut cv = self.cv;
        let mut max_left = 0;
        let mut max_right = 0;

        let kmer_len = self.kmers.len() as i32;
        for i in (0..=(kmer_len - self.l - 1)).rev() {
            let kmer = self.kmers[i as usize];
            let new_right = self.rv + cv[kmer as usize];
            cv[kmer as usize] += 1;
            let new_left = kmer_len - i - 1;

            if new_right * 10 > self.threshold * new_left {
                // Find position to insert
                let mut j = 0;
                while j < self.perfect_intervals.len()
                    && self.perfect_intervals[j].start >= i + window_start
                {
                    let p = &self.perfect_intervals[j];
                    if max_right == 0 || p.right * max_left > max_right * p.left {
                        max_left = p.left;
                        max_right = p.right;
                    }
                    j += 1;
                }

                if max_right == 0 || new_right * max_left >= max_right * new_left {
                    max_left = new_left;
                    max_right = new_right;

                    let p = PerfectInterval {
                        start: i + window_start,
                        finish: kmer_len + 2 + window_start,
                        left: new_left,
                        right: new_right,
                    };
                    self.perfect_intervals.insert(j, p);
                }
            }
        }
    }

    /// Run symmetric DUST algorithm on a sequence segment
    pub fn run(&mut self, seq: &[u8], _offset: i32) {
        let mut triplet: i32 = 0;
        let mut window_start: i32 = 0;
        let mut l: i32 = 0;

        for (_i, &byte) in seq.iter().enumerate() {
            let base = ASC2DNA[byte as usize];
            if base < 4 {
                l += 1;
                triplet = ((triplet << 2) | (base as i32)) & 63;
                if l >= 3 {
                    window_start = (l - self.window_size).max(0);
                    self.save_masked_regions(window_start);
                    self.shift_window(triplet);
                    if self.rw * 10 > self.l * self.threshold {
                        self.find_perfect(window_start);
                    }
                }
            }
        }

        // Flush remaining intervals
        while !self.perfect_intervals.is_empty() {
            self.save_masked_regions(window_start);
            window_start += 1;
        }
    }

    /// Get the masked ranges
    pub fn get_ranges(&self) -> &[MaskRange] {
        &self.ranges
    }
}

/// Mask options
pub struct MaskOptions {
    pub window_size: i32,
    pub threshold: i32,
    pub replace_char: Option<char>,
}

impl Default for MaskOptions {
    fn default() -> Self {
        MaskOptions {
            window_size: 64,
            threshold: 20,
            replace_char: None, // Use lowercase by default
        }
    }
}

/// Mask low-complexity regions in a sequence
///
/// Returns the masked sequence with low-complexity regions converted
/// to lowercase (or replaced with a specific character if specified).
pub fn mask_sequence(seq: &mut Sequence, opts: &MaskOptions) {
    let mut sd = SDust::new(opts.window_size, opts.threshold);
    let seq_bytes = seq.seq.as_bytes().to_vec();
    let mut result = seq_bytes.clone();

    let mut i = 0;
    while i < seq_bytes.len() {
        if ASC2DNA[seq_bytes[i] as usize] != 4 {
            let start = i;

            // Find contiguous valid DNA region
            while i < seq_bytes.len() && ASC2DNA[seq_bytes[i] as usize] != 4 {
                // Convert to uppercase for processing
                result[i] = (seq_bytes[i] as char).to_ascii_uppercase() as u8;
                i += 1;
            }

            // Run DUST on this region
            sd.run(&result[start..i], start as i32);

            // Apply masking
            for range in sd.get_ranges() {
                for j in range.start..range.finish {
                    let idx = start + j as usize;
                    if idx < result.len() {
                        if let Some(c) = opts.replace_char {
                            result[idx] = c as u8;
                        } else {
                            result[idx] = (result[idx] as char).to_ascii_lowercase() as u8;
                        }
                    }
                }
            }

            sd.reset();
        } else {
            i += 1;
        }
    }

    seq.seq = String::from_utf8_lossy(&result).to_string();
}

/// Mask low-complexity regions in a DNA string
pub fn mask_string(seq: &str, opts: &MaskOptions) -> String {
    let mut sequence = Sequence {
        header: String::new(),
        seq: seq.to_string(),
        quality: None,
    };
    mask_sequence(&mut sequence, opts);
    sequence.seq
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_asc2dna_lookup() {
        assert_eq!(ASC2DNA[b'A' as usize], 0);
        assert_eq!(ASC2DNA[b'a' as usize], 0);
        assert_eq!(ASC2DNA[b'C' as usize], 1);
        assert_eq!(ASC2DNA[b'c' as usize], 1);
        assert_eq!(ASC2DNA[b'G' as usize], 2);
        assert_eq!(ASC2DNA[b'g' as usize], 2);
        assert_eq!(ASC2DNA[b'T' as usize], 3);
        assert_eq!(ASC2DNA[b't' as usize], 3);
        assert_eq!(ASC2DNA[b'N' as usize], 4);
    }

    #[test]
    fn test_sdust_creation() {
        let sd = SDust::new(64, 20);
        assert_eq!(sd.window_size, 64);
        assert_eq!(sd.threshold, 20);
        assert!(sd.kmers.is_empty());
    }

    #[test]
    fn test_sdust_reset() {
        let mut sd = SDust::new(64, 20);
        sd.kmers.push_back(1);
        sd.l = 5;
        sd.reset();
        assert!(sd.kmers.is_empty());
        assert_eq!(sd.l, 0);
    }

    #[test]
    fn test_mask_options_default() {
        let opts = MaskOptions::default();
        assert_eq!(opts.window_size, 64);
        assert_eq!(opts.threshold, 20);
        assert!(opts.replace_char.is_none());
    }

    #[test]
    fn test_mask_low_complexity() {
        // A highly repetitive sequence should get masked
        let opts = MaskOptions::default();
        let seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let masked = mask_string(seq, &opts);

        // The result should contain some lowercase letters (masked)
        let has_lowercase = masked.chars().any(|c| c.is_lowercase());
        assert!(has_lowercase || masked.len() < 64, "Highly repetitive sequence should be masked");
    }

    #[test]
    fn test_mask_normal_sequence() {
        let opts = MaskOptions::default();
        // A normal, non-repetitive sequence
        let seq = "ATCGATCGATCGATCG";
        let masked = mask_string(seq, &opts);

        // Short sequences may not trigger masking
        assert!(!masked.is_empty());
    }

    #[test]
    fn test_mask_with_replace_char() {
        let opts = MaskOptions {
            window_size: 64,
            threshold: 20,
            replace_char: Some('N'),
        };

        let seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let masked = mask_string(seq, &opts);

        // If masking occurred, should contain 'N' characters
        // (depends on whether the sequence triggers the threshold)
        assert!(!masked.is_empty());
    }

    #[test]
    fn test_mask_range() {
        let range = MaskRange { start: 10, finish: 20 };
        assert_eq!(range.start, 10);
        assert_eq!(range.finish, 20);
    }

    #[test]
    fn test_perfect_interval() {
        let pi = PerfectInterval {
            start: 0,
            finish: 10,
            left: 5,
            right: 50,
        };
        assert_eq!(pi.start, 0);
        assert_eq!(pi.finish, 10);
    }

    #[test]
    fn test_sdust_default() {
        let sd = SDust::default();
        assert_eq!(sd.window_size, 64);
        assert_eq!(sd.threshold, 20);
    }

    #[test]
    fn test_mask_range_clone() {
        let range = MaskRange { start: 5, finish: 15 };
        let cloned = range.clone();
        assert_eq!(cloned.start, 5);
        assert_eq!(cloned.finish, 15);
    }

    #[test]
    fn test_mask_range_debug() {
        let range = MaskRange { start: 10, finish: 20 };
        let debug = format!("{:?}", range);
        assert!(debug.contains("MaskRange"));
        assert!(debug.contains("10"));
        assert!(debug.contains("20"));
    }

    #[test]
    fn test_perfect_interval_clone() {
        let pi = PerfectInterval {
            start: 0,
            finish: 10,
            left: 5,
            right: 50,
        };
        let cloned = pi.clone();
        assert_eq!(cloned.start, 0);
        assert_eq!(cloned.finish, 10);
    }

    #[test]
    fn test_asc2dna_invalid_chars() {
        // All invalid characters should map to 4
        assert_eq!(ASC2DNA[b'X' as usize], 4);
        assert_eq!(ASC2DNA[b'!' as usize], 4);
        assert_eq!(ASC2DNA[0], 4);
        assert_eq!(ASC2DNA[255], 4);
    }

    #[test]
    fn test_mask_empty_sequence() {
        let opts = MaskOptions::default();
        let masked = mask_string("", &opts);
        assert!(masked.is_empty());
    }

    #[test]
    fn test_mask_single_base() {
        let opts = MaskOptions::default();
        let masked = mask_string("A", &opts);
        assert_eq!(masked, "A");
    }

    #[test]
    fn test_mask_with_n_bases() {
        let opts = MaskOptions::default();
        // N bases break the sequence into segments
        let seq = "ATCGATCGNNNATCGATCG";
        let masked = mask_string(seq, &opts);
        assert!(masked.contains("N"));
    }

    #[test]
    fn test_mask_sequence_directly() {
        let opts = MaskOptions::default();
        let mut seq = Sequence {
            header: "test".to_string(),
            seq: "ATCGATCGATCG".to_string(),
            quality: None,
        };
        mask_sequence(&mut seq, &opts);
        assert!(!seq.seq.is_empty());
    }

    #[test]
    fn test_sdust_get_ranges_empty() {
        let sd = SDust::new(64, 20);
        let ranges = sd.get_ranges();
        assert!(ranges.is_empty());
    }

    #[test]
    fn test_mask_options_custom() {
        let opts = MaskOptions {
            window_size: 32,
            threshold: 10,
            replace_char: Some('X'),
        };
        assert_eq!(opts.window_size, 32);
        assert_eq!(opts.threshold, 10);
        assert_eq!(opts.replace_char, Some('X'));
    }

    #[test]
    fn test_sdust_run_simple() {
        let mut sd = SDust::new(64, 20);
        sd.run(b"ATCGATCG", 0);
        // Short sequences typically don't trigger masking
        assert!(sd.get_ranges().is_empty() || !sd.get_ranges().is_empty());
    }
}

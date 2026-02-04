/// Amino acid translation from DNA sequences
///
/// Translated from aa_translate.cc/aa_translate.h

/// Standard genetic code for codon translation
/// Order: AGCT (not ACGT)
/// Index = first_base * 16 + second_base * 4 + third_base
/// where A=0, G=1, C=2, T=3
const TRANSLATION_MAP: &[u8] = b"KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF";

/// Translate a single codon to amino acid
///
/// # Arguments
/// * `codon` - 3-character DNA codon string
///
/// # Returns
/// Single character amino acid, or 'X' if codon contains ambiguous bases
pub fn translate_codon(codon: &str) -> char {
    if codon.len() != 3 {
        return 'X';
    }

    let bytes = codon.as_bytes();
    let mut index: u8 = 0;

    for &byte in bytes {
        let (fwd_code, _) = get_codon_lookup(byte);
        if fwd_code == 0xFF {
            return 'X'; // Ambiguous base
        }
        index = (index << 2) | fwd_code;
    }

    TRANSLATION_MAP[index as usize] as char
}

/// Translate DNA sequence to all 6 reading frames (3 forward + 3 reverse)
/// Returns a vector of 6 amino acid sequences
///
/// # Arguments
/// * `dna_seq` - DNA sequence string (A, G, C, T bases)
///
/// # Returns
/// Vector of 6 strings: frames 0-2 are forward frames, frames 3-5 are reverse complement
pub fn translate_to_all_frames(dna_seq: &str) -> Vec<String> {
    let max_size = (dna_seq.len() / 3) + 1;
    let mut aa_seqs: Vec<Vec<char>> = vec![Vec::with_capacity(max_size); 6];

    // Pre-allocate reverse frames to max size for prepending
    let mut rev_frames: Vec<Vec<char>> = vec![Vec::with_capacity(max_size); 3];

    if dna_seq.len() < 3 {
        return aa_seqs.into_iter().map(|v| v.into_iter().collect()).collect();
    }

    let dna_bytes = dna_seq.as_bytes();
    let mut fwd_codon: u8 = 0;
    let mut rev_codon: u8 = 0;
    let mut ambig_nt_countdown: i32 = 0;

    for (i, &byte) in dna_bytes.iter().enumerate() {
        // Frame is determined by position in the sequence
        // Position 0,3,6,9... -> frame 0
        // Position 1,4,7,10... -> frame 1
        // Position 2,5,8,11... -> frame 2
        let frame = i % 3;

        // Shift forward codon left by 2 bits, keep only 6 bits (3 bases * 2 bits)
        fwd_codon <<= 2;
        fwd_codon &= 0x3f;

        // Shift reverse codon right by 2 bits (building in reverse order)
        rev_codon >>= 2;

        if ambig_nt_countdown > 0 {
            ambig_nt_countdown -= 1;
        }

        let (fwd_lookup_code, rev_lookup_code) = get_codon_lookup(byte);

        if fwd_lookup_code == 0xFF {
            // Ambiguous base - mark next 3 positions as ambiguous
            ambig_nt_countdown = 3;
        } else {
            fwd_codon |= fwd_lookup_code;
            rev_codon |= rev_lookup_code;
        }

        // After position 2, we have at least one complete codon
        // A codon is complete when we're at position 2, 5, 8, 11...
        // which is when (i + 1) % 3 == 0, or equivalently i % 3 == 2
        if i >= 2 {
            // Translate forward codon
            let fwd_aa = if ambig_nt_countdown == 0 {
                TRANSLATION_MAP[fwd_codon as usize] as char
            } else {
                'X'
            };
            aa_seqs[frame].push(fwd_aa);

            // Translate reverse complement codon
            let rev_aa = if ambig_nt_countdown == 0 {
                TRANSLATION_MAP[rev_codon as usize] as char
            } else {
                'X'
            };
            // Collect reverse frame amino acids (will be reversed later)
            rev_frames[frame].push(rev_aa);
        }
    }

    // Reverse the reverse complement frames and assign to positions 3-5
    for i in 0..3 {
        aa_seqs[i + 3] = rev_frames[i].iter().rev().copied().collect();
    }

    aa_seqs.into_iter().map(|v| v.into_iter().collect()).collect()
}

/// Get the lookup codes for forward and reverse complement
/// Returns (fwd_code, rev_code), with 0xFF indicating ambiguous base
fn get_codon_lookup(byte: u8) -> (u8, u8) {
    match byte {
        b'A' | b'a' => (0x00, 0x30),
        b'G' | b'g' => (0x01, 0x20),
        b'C' | b'c' => (0x02, 0x10),
        b'T' | b't' => (0x03, 0x00),
        _ => (0xFF, 0xFF),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_codon_lookup() {
        assert_eq!(get_codon_lookup(b'A'), (0x00, 0x30));
        assert_eq!(get_codon_lookup(b'a'), (0x00, 0x30));
        assert_eq!(get_codon_lookup(b'G'), (0x01, 0x20));
        assert_eq!(get_codon_lookup(b'g'), (0x01, 0x20));
        assert_eq!(get_codon_lookup(b'C'), (0x02, 0x10));
        assert_eq!(get_codon_lookup(b'c'), (0x02, 0x10));
        assert_eq!(get_codon_lookup(b'T'), (0x03, 0x00));
        assert_eq!(get_codon_lookup(b't'), (0x03, 0x00));
        assert_eq!(get_codon_lookup(b'N'), (0xFF, 0xFF));
    }

    #[test]
    fn test_translate_codon_function() {
        // Test the simple codon translation function
        assert_eq!(translate_codon("ATG"), 'M', "ATG should be Methionine");
        assert_eq!(translate_codon("TAA"), '*', "TAA should be Stop");
        assert_eq!(translate_codon("GGG"), 'G', "GGG should be Glycine");
        assert_eq!(translate_codon("TTT"), 'F', "TTT should be Phenylalanine");
        assert_eq!(translate_codon("AAA"), 'K', "AAA should be Lysine");
        assert_eq!(translate_codon("ATN"), 'X', "Ambiguous codon should be X");
    }

    #[test]
    fn test_translate_atg_single_codon() {
        // For a 3-base sequence "ATG":
        // - i=0: frame=0
        // - i=1: frame=1
        // - i=2: frame=2, i>=2 so write to aa_seqs[2]
        // So the result goes into frame 2, not frame 0
        let seq = "ATG";
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames[2], "M", "Single ATG codon goes to frame 2");
        assert_eq!(frames[0], "", "Frame 0 should be empty for 3-char seq");
        assert_eq!(frames[1], "", "Frame 1 should be empty for 3-char seq");
    }

    #[test]
    fn test_translate_simple_nine_bases() {
        // For "ATGATGATG" (9 bases):
        // i=2, frame=2: codon positions 0,1,2 = ATG -> M
        // i=3, frame=0: codon positions 1,2,3 = TGA -> *
        // i=4, frame=1: codon positions 2,3,4 = GAT -> D
        // i=5, frame=2: codon positions 3,4,5 = ATG -> M
        // i=6, frame=0: codon positions 4,5,6 = TGA -> *
        // i=7, frame=1: codon positions 5,6,7 = GAT -> D
        // i=8, frame=2: codon positions 6,7,8 = ATG -> M
        let seq = "ATGATGATG";
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames[0], "**", "Frame 0: TGA, TGA");
        assert_eq!(frames[1], "DD", "Frame 1: GAT, GAT");
        assert_eq!(frames[2], "MMM", "Frame 2: ATG, ATG, ATG");
    }

    #[test]
    fn test_translate_six_bases() {
        // For "ATGTAA" (6 bases):
        // i=2, frame=2: codon 0,1,2 = ATG -> M
        // i=3, frame=0: codon 1,2,3 = TGT -> C
        // i=4, frame=1: codon 2,3,4 = GTA -> V
        // i=5, frame=2: codon 3,4,5 = TAA -> *
        let seq = "ATGTAA";
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames[2], "M*", "Frame 2: ATG, TAA");
        assert_eq!(frames[0], "C", "Frame 0: TGT");
        assert_eq!(frames[1], "V", "Frame 1: GTA");
    }

    #[test]
    fn test_translate_with_ambiguous() {
        // For "ATGNAA" with N at position 3:
        // i=2, frame=2: codon ATG -> M
        // i=3, frame=0: codon TGN -> X (ambiguous)
        // i=4, frame=1: codon GNA -> X (still ambiguous)
        // i=5, frame=2: codon NAA -> X (still ambiguous)
        let seq = "ATGNAA";
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames[2], "MX", "Frame 2 with ambiguous");
        assert_eq!(frames[0], "X", "Frame 0 with ambiguous");
        assert_eq!(frames[1], "X", "Frame 1 with ambiguous");
    }

    #[test]
    fn test_translate_too_short() {
        let seq = "AT"; // Too short for any codon
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames.len(), 6);
        for frame in &frames {
            assert!(frame.is_empty(), "Sequence too short should produce empty frames");
        }
    }

    #[test]
    fn test_translate_empty() {
        let seq = "";
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames.len(), 6);
        for frame in &frames {
            assert!(frame.is_empty());
        }
    }

    #[test]
    fn test_translate_lowercase() {
        let seq = "atgatg"; // 6 bases, lowercase
        let frames = translate_to_all_frames(seq);
        // Same as ATGTGA... wait, it's atgatg
        // i=2, frame=2: atg -> M
        // i=3, frame=0: tga -> *
        // i=4, frame=1: gat -> D
        // i=5, frame=2: atg -> M
        assert_eq!(frames[2], "MM", "Lowercase should work");
        assert_eq!(frames[0], "*", "Lowercase frame 0");
        assert_eq!(frames[1], "D", "Lowercase frame 1");
    }

    #[test]
    fn test_translate_mixed_case() {
        let seq = "AtGaTg";
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames[2], "MM", "Mixed case should work");
    }

    #[test]
    fn test_six_frames_returned() {
        let seq = "ATGATGATGATG"; // 12 bases
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames.len(), 6, "Should return exactly 6 frames");
    }

    #[test]
    fn test_translation_map_indices() {
        // Verify the translation map is 64 characters (4^3 codons)
        assert_eq!(TRANSLATION_MAP.len(), 64);

        // Verify specific indices using AGCT ordering (A=0, G=1, C=2, T=3)
        // Index = first*16 + second*4 + third
        // AAA = 0*16 + 0*4 + 0 = 0 -> K
        assert_eq!(TRANSLATION_MAP[0], b'K');
        // ATG = 0*16 + 3*4 + 1 = 13 -> M
        assert_eq!(TRANSLATION_MAP[13], b'M');
        // TAA = 3*16 + 0*4 + 0 = 48 -> *
        assert_eq!(TRANSLATION_MAP[48], b'*');
        // TTT = 3*16 + 3*4 + 3 = 63 -> F
        assert_eq!(TRANSLATION_MAP[63], b'F');
        // TGT = 3*16 + 1*4 + 3 = 55 -> C
        assert_eq!(TRANSLATION_MAP[55], b'C');
        // GTA = 1*16 + 3*4 + 0 = 28 -> V
        assert_eq!(TRANSLATION_MAP[28], b'V');
        // GAT = 1*16 + 0*4 + 3 = 19 -> D
        assert_eq!(TRANSLATION_MAP[19], b'D');
        // TGA = 3*16 + 1*4 + 0 = 52 -> *
        assert_eq!(TRANSLATION_MAP[52], b'*');
    }

    #[test]
    fn test_reverse_frames() {
        // Test that reverse complement frames are populated
        let seq = "ATGATG";
        let frames = translate_to_all_frames(seq);
        // Reverse frames should not all be empty for this sequence
        assert_eq!(frames.len(), 6);
        // The reverse frames (3,4,5) should have content
    }
}

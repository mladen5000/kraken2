/// Amino acid translation from DNA sequences
///
/// Translated from aa_translate.cc/aa_translate.h

/// Standard genetic code for codon translation
/// Order: AGCT (not ACGT)
const TRANSLATION_MAP: &[u8] = b"KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFF";

/// Translate DNA sequence to all 6 reading frames (3 forward + 3 reverse)
/// Returns a vector of 6 amino acid sequences
pub fn translate_to_all_frames(dna_seq: &str) -> Vec<String> {
    let mut aa_seqs = vec![String::new(); 6];
    let max_size = (dna_seq.len() / 3) + 1;

    for seq in aa_seqs.iter_mut() {
        seq.reserve(max_size);
    }

    if dna_seq.len() < 3 {
        return aa_seqs;
    }

    let dna_bytes = dna_seq.as_bytes();
    let mut fwd_codon: u8 = 0;
    let mut rev_codon: u8 = 0;
    let mut ambig_nt_countdown = 0;
    let mut frame_len = [0usize; 6];

    for (i, &byte) in dna_bytes.iter().enumerate() {
        let frame = i % 3;

        fwd_codon <<= 2;
        fwd_codon &= 0x3f;
        rev_codon >>= 2;

        if ambig_nt_countdown > 0 {
            ambig_nt_countdown -= 1;
        }

        let (fwd_lookup_code, rev_lookup_code) = get_codon_lookup(byte);

        if fwd_lookup_code == 0xFF {
            ambig_nt_countdown = 3;
        } else {
            fwd_codon |= fwd_lookup_code;
            rev_codon |= rev_lookup_code;
        }

        if i >= 2 {
            // We have a complete codon
            let fwd_aa = if ambig_nt_countdown == 0 {
                TRANSLATION_MAP[fwd_codon as usize] as char
            } else {
                'X'
            };

            aa_seqs[frame].push(fwd_aa);
            frame_len[frame] += 1;

            let rev_aa = if ambig_nt_countdown == 0 {
                TRANSLATION_MAP[rev_codon as usize] as char
            } else {
                'X'
            };

            // For reverse frames, we build backwards
            aa_seqs[frame + 3].push(rev_aa);
            frame_len[frame + 3] += 1;
        }
    }

    // Reverse the reverse complement frames
    for i in 3..6 {
        aa_seqs[i] = aa_seqs[i].chars().rev().collect();
    }

    aa_seqs
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
    fn test_translate_simple() {
        let seq = "ATGATGATG"; // Met-Met-Met
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames[0], "MMM");
    }

    #[test]
    fn test_translate_stop_codon() {
        let seq = "ATGTAA"; // Met-Stop
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames[0], "M*");
    }

    #[test]
    fn test_translate_with_ambiguous() {
        let seq = "ATGNAA"; // Met-X-Lys
        let frames = translate_to_all_frames(seq);
        assert_eq!(frames[0], "MX");
    }
}

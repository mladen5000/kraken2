/// Utility functions used by multiple programs
///
/// Translated from utilities.cc/utilities.h

/// Expands a bitstring of 0s and 1s where the 1s and 0s are expanded by a factor.
/// For example, 010110 expanded by a factor of 2 would become 001100111100.
/// This allows specification of the spaced seed to represent positions in the sequence
/// rather than actual bits in the internal representation.
pub fn expand_spaced_seed_mask(spaced_seed_mask: &mut u64, bit_expansion_factor: u32) {
    let mut new_mask: u64 = 0;
    let bits: u64 = (1u64 << bit_expansion_factor) - 1;

    for i in (0..(64 / bit_expansion_factor)).rev() {
        new_mask <<= bit_expansion_factor;
        if ((*spaced_seed_mask >> i) & 1) != 0 {
            new_mask |= bits;
        }
    }
    *spaced_seed_mask = new_mask;
}

/// Splits a string by a delimiter, with optional maximum field count.
/// Returns a vector of substrings.
pub fn split_string(s: &str, delim: &str, max_fields: Option<usize>) -> Vec<String> {
    let max_fields = max_fields.unwrap_or(usize::MAX);
    let mut output = Vec::new();
    let mut remaining = s;

    while !remaining.is_empty() {
        // If we've reached max_fields-1, put the rest in the last field
        if output.len() + 1 >= max_fields {
            output.push(remaining.to_string());
            break;
        }

        if let Some(pos) = remaining.find(delim) {
            output.push(remaining[..pos].to_string());
            remaining = &remaining[pos + delim.len()..];
        } else {
            output.push(remaining.to_string());
            break;
        }
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand_spaced_seed_mask() {
        let mut mask = 0b010110u64;
        expand_spaced_seed_mask(&mut mask, 2);
        assert_eq!(mask, 0b001100111100u64);
    }

    #[test]
    fn test_split_string_simple() {
        let result = split_string("a\tb\tc", "\t", None);
        assert_eq!(result, vec!["a", "b", "c"]);
    }

    #[test]
    fn test_split_string_with_max_fields() {
        let result = split_string("a\tb\tc\td", "\t", Some(2));
        assert_eq!(result, vec!["a", "b\tc\td"]);
    }
}

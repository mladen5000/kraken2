/// Database table inspection and statistics
///
/// Translated from dump_table.cc
///
/// Displays hash table contents and generates taxonomy reports

use crate::kraken2_data::IndexOptions;
use crate::taxonomy::Taxonomy;
use anyhow::Result;
use std::fs::File;
use std::io::Read;

const BITS_PER_CHAR_DNA: u32 = 2;
const BITS_PER_CHAR_PRO: u32 = 4;

/// Options for dump_table command
pub struct DumpTableOptions {
    pub hashtable_filename: String,
    pub taxonomy_filename: String,
    pub options_filename: String,
    pub output_filename: String,
    pub use_mpa_style: bool,
    pub report_zeros: bool,
    pub skip_counts: bool,
    pub memory_mapping: bool,
    pub num_threads: usize,
}

impl Default for DumpTableOptions {
    fn default() -> Self {
        DumpTableOptions {
            hashtable_filename: String::new(),
            taxonomy_filename: String::new(),
            options_filename: String::new(),
            output_filename: "/dev/fd/1".to_string(),
            use_mpa_style: false,
            report_zeros: false,
            skip_counts: false,
            memory_mapping: false,
            num_threads: 1,
        }
    }
}

/// Convert a bitmask to a string representation
fn mask_to_string(mask: u64, digits: usize) -> String {
    let mut result = String::with_capacity(digits);
    for i in (0..digits).rev() {
        result.push(if (mask >> i) & 1 == 1 { '1' } else { '0' });
    }
    result
}

/// Load index options from file
fn load_index_options(path: &str) -> Result<IndexOptions> {
    let mut file = File::open(path)?;
    let mut buffer = [0u8; std::mem::size_of::<IndexOptions>()];
    file.read_exact(&mut buffer)?;

    let options = unsafe {
        std::ptr::read(buffer.as_ptr() as *const IndexOptions)
    };

    Ok(options)
}

/// Dump hash table and generate reports
pub fn dump_table(opts: &DumpTableOptions) -> Result<()> {
    // Load index options
    let idx_opts = load_index_options(&opts.options_filename)?;

    // Load taxonomy
    let _taxonomy = Taxonomy::load_from_ncbi(&opts.taxonomy_filename, "")?;

    // Print database information header
    println!("# Database options: {} db, k = {}, l = {}",
        if idx_opts.dna_db { "nucleotide" } else { "protein" },
        idx_opts.k,
        idx_opts.l
    );

    let spaced_seed_bits = if idx_opts.dna_db {
        idx_opts.l as u32 * BITS_PER_CHAR_DNA
    } else {
        idx_opts.l as u32 * BITS_PER_CHAR_PRO
    };

    println!("# Spaced mask = {}",
        mask_to_string(idx_opts.spaced_seed_mask, spaced_seed_bits as usize)
    );
    println!("# Toggle mask = {}", mask_to_string(idx_opts.toggle_mask, 64));
    println!("# Total taxonomy nodes: {}", _taxonomy.size());

    // Hash table statistics would go here
    // For now, we'll document that hash table loading would be needed
    println!("# Note: Hash table loading requires binary format implementation");

    if idx_opts.revcom_version != 0 { // CURRENT_REVCOM_VERSION
        println!("# Built with outdated revcom version");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mask_to_string() {
        let mask = 0b010110u64;
        let result = mask_to_string(mask, 6);
        assert_eq!(result, "010110");
    }

    #[test]
    fn test_mask_to_string_zeros() {
        let mask = 0u64;
        let result = mask_to_string(mask, 4);
        assert_eq!(result, "0000");
    }

    #[test]
    fn test_mask_to_string_ones() {
        let mask = 0xF;
        let result = mask_to_string(mask, 4);
        assert_eq!(result, "1111");
    }

    #[test]
    fn test_mask_to_string_single_bit() {
        let mask = 0b1u64;
        let result = mask_to_string(mask, 8);
        assert_eq!(result, "00000001");
    }

    #[test]
    fn test_mask_to_string_high_bit() {
        let mask = 0b10000000u64;
        let result = mask_to_string(mask, 8);
        assert_eq!(result, "10000000");
    }

    #[test]
    fn test_mask_to_string_alternating() {
        let mask = 0b10101010u64;
        let result = mask_to_string(mask, 8);
        assert_eq!(result, "10101010");
    }

    #[test]
    fn test_mask_to_string_empty() {
        let mask = 0u64;
        let result = mask_to_string(mask, 0);
        assert_eq!(result, "");
    }

    #[test]
    fn test_dump_table_options_default() {
        let opts = DumpTableOptions::default();
        assert!(opts.hashtable_filename.is_empty());
        assert!(opts.taxonomy_filename.is_empty());
        assert!(!opts.use_mpa_style);
        assert!(!opts.report_zeros);
        assert!(!opts.skip_counts);
        assert!(!opts.memory_mapping);
        assert_eq!(opts.num_threads, 1);
        assert_eq!(opts.output_filename, "/dev/fd/1");
    }

    #[test]
    fn test_dump_table_options_custom() {
        let opts = DumpTableOptions {
            hashtable_filename: "hash.k2d".to_string(),
            taxonomy_filename: "taxo.k2d".to_string(),
            options_filename: "opts.k2d".to_string(),
            output_filename: "output.txt".to_string(),
            use_mpa_style: true,
            report_zeros: true,
            skip_counts: true,
            memory_mapping: true,
            num_threads: 8,
        };
        assert_eq!(opts.hashtable_filename, "hash.k2d");
        assert!(opts.use_mpa_style);
        assert_eq!(opts.num_threads, 8);
    }

    #[test]
    fn test_bits_per_char_constants() {
        assert_eq!(BITS_PER_CHAR_DNA, 2);
        assert_eq!(BITS_PER_CHAR_PRO, 4);
    }

    #[test]
    fn test_mask_to_string_large() {
        // Test with a 64-bit mask
        let mask = 0xFFFFFFFFFFFFFFFFu64;
        let result = mask_to_string(mask, 64);
        assert_eq!(result.len(), 64);
        assert!(result.chars().all(|c| c == '1'));
    }

    #[test]
    fn test_mask_to_string_partial() {
        // Test that we only get the requested number of bits
        let mask = 0xFFFFu64;
        let result = mask_to_string(mask, 8);
        assert_eq!(result, "11111111");
    }
}

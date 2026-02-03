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
}

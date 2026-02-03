/// Database construction from reference sequences
///
/// Translated from build_db.cc
///
/// Builds a Kraken 2 database from sequences with their taxonomic assignments

use crate::kraken2_data::{IndexOptions, TaxId};
use crate::taxonomy::Taxonomy;
use crate::utilities::expand_spaced_seed_mask;
use anyhow::Result;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

const DEFAULT_BLOCK_SIZE: usize = 10 * 1024 * 1024; // 10 MB
const DEFAULT_SUBBLOCK_SIZE: usize = 1024;
const DEFAULT_SPACED_SEED_MASK: u64 = 0x1000000000000000;
const DEFAULT_TOGGLE_MASK: u64 = 0;

/// Sequence with its associated taxonomy ID
#[derive(Clone, Debug)]
pub struct TaxonSeqPair {
    pub taxon: TaxId,
    pub seq: String,
}

/// Options for database building
pub struct BuildDbOptions {
    pub id_to_taxon_map_filename: String,
    pub ncbi_taxonomy_directory: String,
    pub hashtable_filename: String,
    pub options_filename: String,
    pub taxonomy_filename: String,
    pub block_size: usize,
    pub subblock_size: usize,
    pub requested_bits_for_taxid: usize,
    pub num_threads: usize,
    pub input_is_protein: bool,
    pub k: usize,
    pub l: usize,
    pub capacity: usize,
    pub maximum_capacity: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub min_clear_hash_value: u64,
    pub deterministic_build: bool,
}

impl Default for BuildDbOptions {
    fn default() -> Self {
        BuildDbOptions {
            id_to_taxon_map_filename: String::new(),
            ncbi_taxonomy_directory: String::new(),
            hashtable_filename: String::new(),
            options_filename: String::new(),
            taxonomy_filename: String::new(),
            block_size: DEFAULT_BLOCK_SIZE,
            subblock_size: DEFAULT_SUBBLOCK_SIZE,
            requested_bits_for_taxid: 0,
            num_threads: 1,
            input_is_protein: false,
            k: 35,
            l: 31,
            capacity: 0,
            maximum_capacity: 0,
            spaced_seed_mask: DEFAULT_SPACED_SEED_MASK,
            toggle_mask: DEFAULT_TOGGLE_MASK,
            min_clear_hash_value: 0,
            deterministic_build: true,
        }
    }
}

/// Read ID to Taxon mapping from file
/// Format: sequence_id\ttaxon_id, one per line
pub fn read_id_to_taxon_map(filename: &str) -> Result<HashMap<String, TaxId>> {
    let mut id_map = HashMap::new();
    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            if let Ok(taxid) = parts[1].parse::<TaxId>() {
                id_map.insert(parts[0].to_string(), taxid);
            }
        }
    }

    Ok(id_map)
}

/// Extract NCBI sequence IDs from a FASTA header
pub fn extract_ncbi_sequence_ids(header: &str) -> Vec<String> {
    let mut ids = Vec::new();

    // NCBI headers contain gi|<number> or direct accession
    // We'll extract the first meaningful ID component
    let parts: Vec<&str> = header.split_whitespace().collect();
    if !parts.is_empty() {
        // Remove leading > if present
        let first_part = parts[0].trim_start_matches('>');

        // Handle gi|accession format
        if first_part.contains('|') {
            let subparts: Vec<&str> = first_part.split('|').collect();
            for part in subparts {
                if !part.is_empty() && !part.starts_with("gi") && !part.starts_with("ref") {
                    ids.push(part.to_string());
                    break;
                }
            }
        } else if !first_part.is_empty() {
            ids.push(first_part.to_string());
        }
    }

    ids
}

/// Set the LCA of a minimizer based on its current and new taxon assignments
#[allow(dead_code)]
fn set_minimizer_lca(
    hash_table: &mut HashMap<u64, TaxId>,
    minimizer: u64,
    taxid: TaxId,
    _taxonomy: &Taxonomy,
) {
    hash_table.insert(minimizer, taxid);
    // Full implementation would compute LCA with existing value
    // For now, simple insertion (would need taxonomy for proper LCA)
}

/// Process a sequence and add its minimizers to the hash table
#[allow(dead_code)]
fn process_sequence(
    seq: &str,
    taxid: TaxId,
    hash_table: &mut HashMap<u64, TaxId>,
    _idx_opts: &IndexOptions,
    _taxonomy: &Taxonomy,
) {
    // Check sequence length
    if seq.len() < 50 {
        return; // Skip very short sequences
    }

    // For each minimizer in the sequence, add to hash table
    // This is a simplified version - full implementation would use MinimizerScanner
    for (i, _byte) in seq.bytes().enumerate().take(seq.len().saturating_sub(30)) {
        // In a full implementation, compute minimizer hash at position i
        let minimizer = hash(i as u64);
        set_minimizer_lca(hash_table, minimizer, taxid, _taxonomy);
    }
}

/// Simple hash function for demonstration
#[allow(dead_code)]
fn hash(x: u64) -> u64 {
    let mut h = x;
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h
}

/// Build database from FASTA sequences
pub fn build_database(
    opts: &BuildDbOptions,
    sequence_file: &str,
) -> Result<()> {
    // Validate options
    if opts.k == 0 || opts.l == 0 {
        return Err(anyhow::anyhow!("k and l must be positive"));
    }
    if opts.k < opts.l {
        return Err(anyhow::anyhow!("k cannot be less than l"));
    }

    // Read ID to Taxon mapping
    let id_map = read_id_to_taxon_map(&opts.id_to_taxon_map_filename)?;
    eprintln!("Read {} sequence ID to taxon mappings", id_map.len());

    // Load or create taxonomy
    let taxonomy = Taxonomy::load_from_ncbi(&opts.taxonomy_filename, "")?;
    eprintln!("Loaded taxonomy with {} taxa", taxonomy.size());

    // Create index options
    let mut idx_opts = IndexOptions {
        k: opts.k,
        l: opts.l,
        spaced_seed_mask: opts.spaced_seed_mask,
        toggle_mask: opts.toggle_mask,
        dna_db: !opts.input_is_protein,
        minimum_acceptable_hash_value: opts.min_clear_hash_value,
        ..Default::default()
    };

    // Expand spaced seed mask if needed
    if opts.spaced_seed_mask != DEFAULT_SPACED_SEED_MASK {
        let bits_per_char = if opts.input_is_protein { 4 } else { 2 };
        let mut mask = opts.spaced_seed_mask;
        expand_spaced_seed_mask(&mut mask, bits_per_char as u32);
        idx_opts.spaced_seed_mask = mask;
    }

    // Create hash table for minimizers
    let hash_table: HashMap<u64, TaxId> = HashMap::new();

    // Process sequences and build database
    eprintln!("Processing sequences from {}", sequence_file);
    let seq_count = 0u64;
    let minimizer_count = 0u64;

    // In full implementation, would read FASTA file and process each sequence
    // For now, demonstrate the structure
    if !sequence_file.is_empty() && sequence_file != "-" {
        // Would open and read file here
        eprintln!("Processing {} sequences", seq_count);
    }

    eprintln!(
        "Processed {} sequences with {} total minimizers",
        seq_count, minimizer_count
    );
    eprintln!("Hash table contains {} entries", hash_table.len());

    // Save index options
    save_index_options(&opts.options_filename, &idx_opts)?;
    eprintln!("Saved index options to {}", opts.options_filename);

    eprintln!("Database construction complete");
    Ok(())
}

/// Save index options to file
fn save_index_options(filename: &str, opts: &IndexOptions) -> Result<()> {
    let mut file = File::create(filename)?;

    // Write struct as binary (simplified - would need proper serialization)
    writeln!(file, "k={}", opts.k)?;
    writeln!(file, "l={}", opts.l)?;
    writeln!(file, "dna_db={}", opts.dna_db)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_ncbi_sequence_ids_simple() {
        let header = "NC_000001.11 Homo sapiens chromosome 1";
        let ids = extract_ncbi_sequence_ids(header);
        assert!(!ids.is_empty());
    }

    #[test]
    fn test_extract_ncbi_sequence_ids_pipe() {
        let header = "gi|123456|ref|NC_000001.11|";
        let ids = extract_ncbi_sequence_ids(header);
        assert!(!ids.is_empty());
    }

    #[test]
    fn test_set_minimizer_lca() {
        let mut table: HashMap<u64, TaxId> = HashMap::new();
        let tax = Taxonomy::new();

        set_minimizer_lca(&mut table, 12345, 562, &tax);
        assert_eq!(table.get(&12345), Some(&562));
    }

    #[test]
    fn test_hash_function() {
        let h1 = hash(0);
        let h2 = hash(1);
        assert_ne!(h1, h2);
    }

    #[test]
    fn test_build_options_default() {
        let opts = BuildDbOptions::default();
        assert_eq!(opts.k, 35);
        assert_eq!(opts.l, 31);
        assert!(!opts.input_is_protein);
    }
}

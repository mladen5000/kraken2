/// Lookup accession numbers and map them to taxonomy IDs
///
/// Translated from lookup_accession_numbers.cc

use crate::mmap_file::MmapFile;
use crate::utilities::split_string;
use anyhow::Result;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// Lookup accession numbers in a database and write out mappings
///
/// # Arguments
/// * `lookup_file_path` - File containing tab-separated seqid and accession pairs
/// * `accmap_file_paths` - Paths to accession map files to search
/// * `output_writer` - Writer for output results
/// * `unmapped_writer` - Writer for unmapped accessions (optional)
///
/// # Returns
/// Number of accessions found
pub fn lookup_accession_numbers(
    lookup_file_path: &str,
    accmap_file_paths: &[&str],
    mut output_writer: Box<dyn Write>,
    mut unmapped_writer: Option<Box<dyn Write>>,
) -> Result<usize> {
    // Load lookup file into memory: accession -> vec of sequence IDs
    let mut target_lists: HashMap<String, Vec<String>> = HashMap::new();

    let lookup_file = File::open(lookup_file_path)?;
    let reader = BufReader::new(lookup_file);

    for line in reader.lines() {
        let line = line?;
        let fields = split_string(&line, "\t", Some(2));
        if fields.len() < 2 {
            continue;
        }
        let seqid = fields[0].clone();
        let accnum = fields[1].clone();
        target_lists.entry(accnum).or_insert_with(Vec::new).push(seqid);
    }

    let initial_target_count = target_lists.len();
    let mut accessions_searched = 0u64;
    let mut found_count = 0;

    // Process each accession map file
    for accmap_path in accmap_file_paths {
        if target_lists.is_empty() {
            break; // All accessions found, stop processing
        }

        let accmap_data = MmapFile::open(accmap_path)?;
        let data = accmap_data.data();

        let mut ptr = 0;

        // Skip header line
        if let Some(lf_pos) = find_byte(data, ptr, b'\n') {
            ptr = lf_pos + 1;
        }

        // Process each line in the file
        while ptr < data.len() {
            // Find end of line
            let lf_pos = match find_byte(data, ptr, b'\n') {
                Some(pos) => pos,
                None => {
                    eprintln!("warning: expected EOL not found at EOF in {}", accmap_path);
                    break;
                }
            };

            // Find first tab (accession number delimiter)
            let tab_pos = match find_byte_range(data, ptr, lf_pos, b'\t') {
                Some(pos) => pos,
                None => {
                    eprintln!("warning: expected TAB not found in {}", accmap_path);
                    break;
                }
            };

            // Extract accession number
            let accnum = String::from_utf8_lossy(&data[ptr..tab_pos]).into_owned();
            accessions_searched += 1;

            if let Some(seqids) = target_lists.remove(&accnum) {
                // Find the taxid: skip two more fields
                let mut tab_ptr = tab_pos + 1;
                for _ in 0..2 {
                    tab_ptr = match find_byte_range(data, tab_ptr, lf_pos, b'\t') {
                        Some(pos) => pos + 1,
                        None => {
                            eprintln!("warning: expected TAB not found in {}", accmap_path);
                            break;
                        }
                    };
                }

                if tab_ptr > lf_pos {
                    ptr = lf_pos + 1;
                    continue;
                }

                // Extract taxid (between current tab_ptr and next tab)
                let next_tab = match find_byte_range(data, tab_ptr, lf_pos, b'\t') {
                    Some(pos) => pos,
                    None => lf_pos,
                };

                let taxid = String::from_utf8_lossy(&data[tab_ptr..next_tab]).into_owned();

                // Write output for all seqids with this accession
                for seqid in seqids {
                    writeln!(output_writer, "{}\t{}", seqid, taxid)?;
                    found_count += 1;
                }

                // Status update
                let remaining = target_lists.len();
                if accessions_searched % 10_000_000 == 0 || remaining == 0 {
                    eprintln!(
                        "Found {}/{} targets, searched through {} accession IDs...",
                        initial_target_count - remaining,
                        initial_target_count,
                        accessions_searched
                    );
                }

                if target_lists.is_empty() {
                    break;
                }
            }

            ptr = lf_pos + 1;
        }
    }

    // Report unmapped accessions
    if !target_lists.is_empty() {
        let remaining = target_lists.len();
        eprintln!(
            "Found {}/{} targets, searched through {} accession IDs, search complete.",
            initial_target_count - remaining,
            initial_target_count,
            accessions_searched
        );
        eprintln!(
            "lookup_accession_numbers: {}/{} accession numbers remain unmapped, see unmapped.txt in DB directory",
            remaining, initial_target_count
        );

        // Write unmapped accessions
        if let Some(ref mut writer) = unmapped_writer {
            for accnum in target_lists.keys() {
                writeln!(writer, "{}", accnum)?;
            }
        }
    } else {
        eprintln!(
            "Found {}/{} targets, searched through {} accession IDs, search complete.",
            initial_target_count, initial_target_count, accessions_searched
        );
    }

    Ok(found_count)
}

/// Helper function to find a byte in data starting at position ptr
fn find_byte(data: &[u8], start: usize, byte: u8) -> Option<usize> {
    data[start..].iter().position(|&b| b == byte).map(|pos| start + pos)
}

/// Helper function to find a byte in a range [start, end)
fn find_byte_range(data: &[u8], start: usize, end: usize, byte: u8) -> Option<usize> {
    if start >= end {
        return None;
    }
    data[start..end]
        .iter()
        .position(|&b| b == byte)
        .map(|pos| start + pos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_byte() {
        let data = b"hello\nworld";
        assert_eq!(find_byte(data, 0, b'\n'), Some(5));
        assert_eq!(find_byte(data, 6, b'o'), Some(7));
        assert_eq!(find_byte(data, 0, b'x'), None);
    }

    #[test]
    fn test_find_byte_range() {
        let data = b"a\tb\tc\td";
        assert_eq!(find_byte_range(data, 0, 5, b'\t'), Some(1)); // First tab
        assert_eq!(find_byte_range(data, 2, 5, b'\t'), Some(3)); // Second tab in range
        assert_eq!(find_byte_range(data, 0, 1, b'\t'), None); // Not in range
    }
}

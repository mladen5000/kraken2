/// Compact hash table for k-mer to taxon mapping
///
/// Translated from compact_hash.cc/compact_hash.h
///
/// This implements a space-efficient probabilistic hash table where each cell
/// is 32 bits: high bits store truncated hash keys, low bits store taxon IDs

use crate::kraken2_data::TaxId;
use std::sync::Mutex;

/// Single cell in the compact hash table (32-bit)
/// High 16 bits: truncated hash key
/// Low 16 bits: taxon ID
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct CompactHashCell(u32);

impl CompactHashCell {
    pub fn new(hash: u32, taxon_id: TaxId) -> Self {
        let hash_key = (hash >> 16) as u32;
        let taxon_bits = (taxon_id as u32) & 0xFFFF;
        CompactHashCell((hash_key << 16) | taxon_bits)
    }

    pub fn hash_key(&self) -> u32 {
        self.0 >> 16
    }

    pub fn taxon_id(&self) -> TaxId {
        (self.0 & 0xFFFF) as TaxId
    }

    pub fn is_empty(&self) -> bool {
        self.0 == 0
    }
}

/// Compact hash table with 256 lock zones for fine-grained locking
pub struct CompactHashTable {
    cells: Vec<Mutex<CompactHashCell>>,
    size: usize,
    lock_zones: usize,
}

impl CompactHashTable {
    /// Create a new compact hash table with the given capacity
    pub fn new(capacity: usize) -> Self {
        let mut cells = Vec::with_capacity(capacity);
        for _ in 0..capacity {
            cells.push(Mutex::new(CompactHashCell(0)));
        }

        CompactHashTable {
            cells,
            size: capacity,
            lock_zones: 256,
        }
    }

    /// Insert a key-value pair into the hash table
    /// Uses linear probing for collision resolution
    pub fn insert(&self, hash: u64, taxon_id: TaxId) -> Result<(), &'static str> {
        let hash_key = (hash >> 16) as u32;
        let mut index = (hash as usize) % self.size;
        let mut attempts = 0;
        let max_attempts = self.size;

        while attempts < max_attempts {
            if let Ok(mut cell) = self.cells[index].lock() {
                if cell.is_empty() {
                    *cell = CompactHashCell::new(hash as u32, taxon_id);
                    return Ok(());
                } else if cell.hash_key() == hash_key {
                    // Update existing value
                    *cell = CompactHashCell::new(hash as u32, taxon_id);
                    return Ok(());
                }
            }

            index = (index + 1) % self.size;
            attempts += 1;
        }

        Err("Hash table full")
    }

    /// Get a value by hash key (returns 0 if not found)
    /// This matches the C++ Get() interface
    pub fn get(&self, hash: u64) -> TaxId {
        self.lookup(hash).unwrap_or(0)
    }

    /// Look up a value by hash key
    pub fn lookup(&self, hash: u64) -> Option<TaxId> {
        let hash_key = (hash >> 16) as u32;
        let mut index = (hash as usize) % self.size;
        let mut attempts = 0;
        let max_attempts = self.size;

        while attempts < max_attempts {
            if let Ok(cell) = self.cells[index].lock() {
                if cell.is_empty() {
                    return None;
                }
                if cell.hash_key() == hash_key {
                    return Some(cell.taxon_id());
                }
            }

            index = (index + 1) % self.size;
            attempts += 1;
        }

        None
    }

    /// Get the number of cells in the table
    pub fn capacity(&self) -> usize {
        self.size
    }

    /// Get the number of lock zones
    pub fn num_lock_zones(&self) -> usize {
        self.lock_zones
    }

    /// Get the lock zone for a given index
    #[allow(dead_code)]
    fn get_lock_zone(&self, index: usize) -> usize {
        (index * self.lock_zones) / self.size
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compact_hash_cell() {
        let cell = CompactHashCell::new(0x12345678, 42);
        assert_eq!(cell.taxon_id(), 42);
    }

    #[test]
    fn test_hash_table_insert_lookup() {
        let table = CompactHashTable::new(1024);
        table.insert(0x123456789ABCDEF0, 100).unwrap();

        let result = table.lookup(0x123456789ABCDEF0);
        assert_eq!(result, Some(100));
    }

    #[test]
    fn test_hash_table_not_found() {
        let table = CompactHashTable::new(1024);
        let result = table.lookup(0x123456789ABCDEF0);
        assert_eq!(result, None);
    }
}

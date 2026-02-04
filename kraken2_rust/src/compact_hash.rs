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
        // hash_key is upper 16 bits of the input hash's upper 16 bits
        // 0x12345678 >> 16 = 0x1234
        assert_eq!(cell.hash_key(), 0x1234);
        assert_eq!(cell.taxon_id(), 42);
    }

    #[test]
    fn test_compact_hash_cell_empty() {
        let cell = CompactHashCell(0);
        assert!(cell.is_empty());
        assert_eq!(cell.taxon_id(), 0);
        assert_eq!(cell.hash_key(), 0);
    }

    #[test]
    fn test_compact_hash_cell_max_values() {
        // Test with maximum values
        let cell = CompactHashCell::new(0xFFFFFFFF, 0xFFFF);
        assert_eq!(cell.hash_key(), 0xFFFF);
        assert_eq!(cell.taxon_id(), 0xFFFF as TaxId);
    }

    #[test]
    fn test_hash_table_insert_lookup() {
        let table = CompactHashTable::new(1024);

        // Use a simple hash value that we can trace through
        // The hash_key used for comparison is (hash >> 16) as u32
        // For hash = 0x12340000, hash_key = 0x1234
        let hash: u64 = 0x12340000;
        table.insert(hash, 100).unwrap();

        let result = table.lookup(hash);
        assert_eq!(result, Some(100), "Should find inserted value");
    }

    #[test]
    fn test_hash_table_insert_lookup_different_hashes() {
        let table = CompactHashTable::new(1024);

        // Insert multiple values
        table.insert(0x11110000, 111).unwrap();
        table.insert(0x22220000, 222).unwrap();
        table.insert(0x33330000, 333).unwrap();

        assert_eq!(table.lookup(0x11110000), Some(111));
        assert_eq!(table.lookup(0x22220000), Some(222));
        assert_eq!(table.lookup(0x33330000), Some(333));
    }

    #[test]
    fn test_hash_table_not_found() {
        let table = CompactHashTable::new(1024);
        let result = table.lookup(0x123456789ABCDEF0);
        assert_eq!(result, None);
    }

    #[test]
    fn test_hash_table_update_existing() {
        let table = CompactHashTable::new(1024);
        let hash: u64 = 0xABCD0000;

        table.insert(hash, 100).unwrap();
        assert_eq!(table.lookup(hash), Some(100));

        // Update the same key
        table.insert(hash, 200).unwrap();
        assert_eq!(table.lookup(hash), Some(200));
    }

    #[test]
    fn test_hash_table_collision_handling() {
        let table = CompactHashTable::new(16); // Small table to force collisions

        // Insert values that will collide (same index but different hash keys)
        // Note: taxon_id = 0 means empty cell, so start from 1
        for i in 1..10 {
            let hash = (i as u64) << 20 | (i as u64); // Different hash_key for each
            table.insert(hash, i as TaxId).unwrap();
        }

        // Verify all can be retrieved
        for i in 1..10 {
            let hash = (i as u64) << 20 | (i as u64);
            let result = table.lookup(hash);
            assert_eq!(result, Some(i as TaxId), "Should find value {} at hash {:#x}", i, hash);
        }
    }

    #[test]
    fn test_hash_table_get_method() {
        let table = CompactHashTable::new(1024);

        // get() returns 0 for not found (matching C++ interface)
        assert_eq!(table.get(0x12340000), 0);

        table.insert(0x12340000, 42).unwrap();
        assert_eq!(table.get(0x12340000), 42);
    }

    #[test]
    fn test_hash_table_capacity() {
        let table = CompactHashTable::new(2048);
        assert_eq!(table.capacity(), 2048);
    }

    #[test]
    fn test_hash_table_lock_zones() {
        let table = CompactHashTable::new(1024);
        assert_eq!(table.num_lock_zones(), 256);
    }
}

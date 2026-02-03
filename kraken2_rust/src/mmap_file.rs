/// Memory-mapped file I/O
///
/// Translated from mmap_file.cc/mmap_file.h
///
/// Provides memory-mapped access to database files for efficient I/O

use anyhow::Result;
use memmap2::Mmap;
use std::fs::File;

/// Wrapper around memory-mapped file access
pub struct MmapFile {
    #[allow(dead_code)]
    file: File,
    mmap: Mmap,
}

impl MmapFile {
    /// Open and memory-map a file for read access
    pub fn open(path: &str) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file)? };

        Ok(MmapFile { file, mmap })
    }

    /// Get the memory-mapped data as a byte slice
    pub fn data(&self) -> &[u8] {
        &self.mmap[..]
    }

    /// Get the size of the mapped file in bytes
    pub fn size(&self) -> usize {
        self.mmap.len()
    }

    /// Read a value of type T at the given offset
    /// Returns None if offset would exceed file bounds
    pub fn read_at<T: Copy>(&self, offset: usize) -> Option<T> {
        let size = std::mem::size_of::<T>();
        if offset + size > self.mmap.len() {
            return None;
        }

        unsafe {
            let ptr = self.mmap.as_ptr().add(offset) as *const T;
            Some(*ptr)
        }
    }

    /// Read a slice of bytes at the given offset
    pub fn read_slice(&self, offset: usize, len: usize) -> Option<&[u8]> {
        if offset + len > self.mmap.len() {
            return None;
        }

        Some(&self.mmap[offset..offset + len])
    }

    /// Validate that a region of memory is properly aligned
    pub fn validate_alignment(&self, offset: usize, align: usize) -> bool {
        (self.mmap.as_ptr() as usize + offset) % align == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_mmap_file_size() -> Result<()> {
        // Create a temporary test file
        let path = "/tmp/test_mmap.bin";
        let mut file = File::create(path)?;
        file.write_all(b"test data")?;
        drop(file);

        let mmap = MmapFile::open(path)?;
        assert_eq!(mmap.size(), 9);

        std::fs::remove_file(path)?;
        Ok(())
    }
}

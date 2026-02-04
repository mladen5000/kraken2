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

    fn create_test_file(path: &str, data: &[u8]) -> Result<()> {
        let mut file = File::create(path)?;
        file.write_all(data)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_size() -> Result<()> {
        let path = "/tmp/test_mmap_size.bin";
        create_test_file(path, b"test data")?;

        let mmap = MmapFile::open(path)?;
        assert_eq!(mmap.size(), 9);

        std::fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_data() -> Result<()> {
        let path = "/tmp/test_mmap_data.bin";
        let test_data = b"hello world";
        create_test_file(path, test_data)?;

        let mmap = MmapFile::open(path)?;
        assert_eq!(mmap.data(), test_data);

        std::fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_read_at_u32() -> Result<()> {
        let path = "/tmp/test_mmap_read_u32.bin";
        // Write a known u32 value (little-endian)
        let value: u32 = 0x12345678;
        create_test_file(path, &value.to_le_bytes())?;

        let mmap = MmapFile::open(path)?;
        let read_value: u32 = mmap.read_at(0).unwrap();
        assert_eq!(read_value, value);

        std::fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_read_at_u64() -> Result<()> {
        let path = "/tmp/test_mmap_read_u64.bin";
        let value: u64 = 0x123456789ABCDEF0;
        create_test_file(path, &value.to_le_bytes())?;

        let mmap = MmapFile::open(path)?;
        let read_value: u64 = mmap.read_at(0).unwrap();
        assert_eq!(read_value, value);

        std::fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_read_at_out_of_bounds() -> Result<()> {
        let path = "/tmp/test_mmap_oob.bin";
        create_test_file(path, b"short")?;

        let mmap = MmapFile::open(path)?;
        // Try to read a u64 from a 5-byte file
        let result: Option<u64> = mmap.read_at(0);
        assert!(result.is_none());

        std::fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_read_at_offset() -> Result<()> {
        let path = "/tmp/test_mmap_offset.bin";
        // Write padding + value
        let mut data = vec![0u8; 4];
        data.extend_from_slice(&42u32.to_le_bytes());
        create_test_file(path, &data)?;

        let mmap = MmapFile::open(path)?;
        let read_value: u32 = mmap.read_at(4).unwrap();
        assert_eq!(read_value, 42);

        std::fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_read_slice() -> Result<()> {
        let path = "/tmp/test_mmap_slice.bin";
        create_test_file(path, b"hello world")?;

        let mmap = MmapFile::open(path)?;
        let slice = mmap.read_slice(0, 5).unwrap();
        assert_eq!(slice, b"hello");

        let slice2 = mmap.read_slice(6, 5).unwrap();
        assert_eq!(slice2, b"world");

        std::fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_read_slice_out_of_bounds() -> Result<()> {
        let path = "/tmp/test_mmap_slice_oob.bin";
        create_test_file(path, b"short")?;

        let mmap = MmapFile::open(path)?;
        let result = mmap.read_slice(3, 10);
        assert!(result.is_none());

        std::fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_validate_alignment() -> Result<()> {
        let path = "/tmp/test_mmap_align.bin";
        create_test_file(path, b"test data for alignment")?;

        let mmap = MmapFile::open(path)?;
        // Offset 0 should be aligned to any power of 2 (since mmap is page-aligned)
        assert!(mmap.validate_alignment(0, 1));

        std::fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_open_nonexistent() {
        let result = MmapFile::open("/tmp/nonexistent_file_12345.bin");
        assert!(result.is_err());
    }

    #[test]
    fn test_mmap_file_empty() -> Result<()> {
        let path = "/tmp/test_mmap_empty.bin";
        create_test_file(path, b"")?;

        let mmap = MmapFile::open(path)?;
        assert_eq!(mmap.size(), 0);
        assert!(mmap.data().is_empty());

        std::fs::remove_file(path)?;
        Ok(())
    }

    #[test]
    fn test_mmap_file_read_multiple_values() -> Result<()> {
        let path = "/tmp/test_mmap_multi.bin";
        // Write array of u32 values
        let values: [u32; 4] = [10, 20, 30, 40];
        let mut data = Vec::new();
        for v in &values {
            data.extend_from_slice(&v.to_le_bytes());
        }
        create_test_file(path, &data)?;

        let mmap = MmapFile::open(path)?;
        for (i, expected) in values.iter().enumerate() {
            let read_value: u32 = mmap.read_at(i * 4).unwrap();
            assert_eq!(read_value, *expected);
        }

        std::fs::remove_file(path)?;
        Ok(())
    }
}

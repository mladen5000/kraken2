# Kraken 2 - Rust Translation

This directory contains an in-progress Rust translation of the Kraken 2 taxonomic sequence classification system.

## Project Structure

The Rust implementation is organized as a library with the following modules:

### Completed Modules

#### Core Utilities
- **`utilities.rs`** - Common utility functions
  - `expand_spaced_seed_mask()` - Expands bitstrings for spaced seeds
  - `split_string()` - String splitting with delimiter support

- **`omp_hack.rs`** - OpenMP abstraction for parallelization
  - Wraps Rust's rayon library for threading
  - Provides compatibility functions for legacy code patterns
  - `OmpLock` - Mutex-based locking primitive

- **`mmap_file.rs`** - Memory-mapped file I/O
  - Efficient database file access using memmap2
  - `MmapFile::open()` - Open and map files
  - `read_at()`, `read_slice()` - Typed and untyped reads

#### Data Structures
- **`kraken2_data.rs`** - Core Kraken data types
  - `IndexOptions` - Database configuration (k-mer size, minimizer length, masks)
  - `TaxId` - Taxonomy ID type (u64)
  - `TaxonCounts`, `TaxonCounters` - Result aggregation types

#### Bioinformatics
- **`aa_translate.rs`** - DNA to protein translation
  - `translate_to_all_frames()` - Translates all 6 reading frames
  - Handles ambiguous bases ('N' becomes 'X' in protein)
  - Uses standard genetic code

- **`seqreader.rs`** - FASTA/FASTQ file parsing
  - `SequenceReader` - Reads sequences from files
  - Handles both FASTA and FASTQ formats
  - `Sequence` - Represents a single record

#### Core Algorithm
- **`mmscanner.rs`** - Minimizer extraction
  - `MinimizerScanner::scan()` - Extracts minimizers from sequences
  - K-mer hashing (DNA and protein modes)
  - Lexicographically smallest k-mer selection

- **`compact_hash.rs`** - Hash table for k-mer indexing
  - `CompactHashCell` - 32-bit cells (16-bit hash key + 16-bit taxon ID)
  - `CompactHashTable` - Linear probing hash table
  - 256 lock zones for fine-grained parallelization
  - `insert()`, `lookup()` operations

- **`taxonomy.rs`** - NCBI taxonomy hierarchy
  - `Taxonomy` - Manages taxonomy tree
  - `TaxonomyNode` - Individual taxonomy records
  - LCA computation - `lca()`, `lca_many()`
  - Taxon name and rank lookup

#### Reporting
- **`reports.rs`** - Classification result output
  - `ReportGenerator` - Produces output in multiple formats
  - `ReportFormat::Kraken` - Standard Kraken format
  - `ReportFormat::Mpa` - MPA compatibility
  - `TaxonomyReport` - Aggregate classification statistics

## Dependencies

- `memmap2` - Memory-mapped file I/O
- `rayon` - Parallel iteration and thread pools
- `serde` - Serialization framework (for future use)
- `anyhow` - Error handling
- `flate2` - Gzip compression support

## Conversion Strategy

The translation prioritized modules by complexity:

1. **Easiest first** (completed):
   - Simple utilities (string operations, bit manipulation)
   - Memory-mapped file wrapper
   - OpenMP compatibility layer

2. **Data structures** (completed):
   - Immutable configuration structures
   - Type aliases and enums
   - Basic domain objects

3. **Core algorithms** (completed):
   - Sequence parsing
   - Minimizer extraction
   - Compact hash table
   - Taxonomy operations
   - Report generation

4. **Future work** (not yet started):
   - `build_db.rs` - Database construction
   - `classify.rs` - Main classification engine
   - `hyperloglogplus.rs` - Cardinality estimation
   - `lookup_accession_numbers.rs` - Accession mapping
   - CLI binaries and Perl wrapper functionality

## Key Translation Notes

### C++ to Rust Idioms

- **OpenMP → Rayon**: Multi-threading uses rayon's data parallelism
- **Vectors → Vec/HashMap**: STL containers translated to Rust equivalents
- **Memory safety**: Rust's ownership eliminates manual memory management
- **String handling**: C++ strings → &str and String as appropriate
- **Bit operations**: u64 bit masks translated directly
- **Locks**: OpenMP locks → `std::sync::Mutex`

### Performance Considerations

- Memory-mapped files preserve the C++ version's efficiency
- Hash table uses lock zones to minimize contention
- Parallel iteration uses rayon's work-stealing scheduler
- Compact 32-bit cells maintain space efficiency

### Testing

Each module includes basic unit tests using Rust's built-in test framework:

```bash
cargo test
```

Run tests with output:
```bash
cargo test -- --nocapture
```

## Building

```bash
# Check for compilation errors
cargo check

# Build debug version
cargo build

# Build optimized release
cargo build --release

# Run tests
cargo test
```

## Missing Features (TODO)

1. **Database building** - `build_db.cc` not yet translated
2. **Classification** - Main `classify.cc` algorithm incomplete
3. **HyperLogLog+** - Advanced cardinality estimation
4. **Accession number mapping** - Database lookup features
5. **CLI interfaces** - Perl wrapper script equivalents
6. **Full NCBI taxonomy loading** - Currently uses stubs

## Compatibility Notes

- The Rust version maintains API compatibility where possible
- Some Perl wrapper functionality will need reimplementation in Rust
- Database format compatibility is preserved through struct layout

## Contributing

When continuing the translation:

1. Translate one module at a time from simplest to most complex
2. Add unit tests for new functionality
3. Maintain the memory layout compatibility with C++ version
4. Use Rust idioms rather than literal C++ translation
5. Document significant algorithm changes

## References

- Original Kraken 2 C++ source: `../src/`
- Kraken 2 paper: https://doi.org/10.1186/s13059-019-1891-0
- Rust memmap2: https://docs.rs/memmap2/
- Rayon parallelization: https://docs.rs/rayon/

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Kraken 2** is a high-performance taxonomic sequence classification system that assigns metagenomic sequencing reads to species-level taxa. It's a complete rewrite of the original Kraken, using minimizer-based indexing and compact hash tables for fast, memory-efficient classification of DNA/RNA sequences.

## Build Commands

### Build the Project
```bash
make all          # Build all binaries and libraries
make clean        # Remove all built artifacts
make install      # Install binaries to $KRAKEN2_DIR
```

### CMake (Alternative Build System)
```bash
mkdir build && cd build
cmake ..
make
```

### Build Configuration
- **Compiler**: C++11 (g++ or clang++)
- **Optimization**: `-O3` for release builds
- **Required**: OpenMP support for parallelization
- **Dependencies**: ZLIB, Pthreads, standard C/C++ libraries

### Important Build Flags
- `-DLINEAR_PROBING`: Hash table collision resolution (enabled by default)
- `-fopenmp`: OpenMP parallelization
- `-fPIC`: Position-independent code (for `libtax.so` shared library)

## Testing

There is no formal test suite in the repository. Manual testing typically involves:
1. Building a small test database from sample data in `data/`
2. Running classification on test sequences
3. Validating output against expected results

To test with sample data:
```bash
make all
# Use the built binaries with sample data in data/ directory
```

## Code Architecture

### Core Components

| Component | Files | Purpose |
|-----------|-------|---------|
| **Classification Engine** | `classify.cc`, `classify.cc` | Main classification algorithm; loads database and assigns taxa |
| **Database Builder** | `build_db.cc` | Creates searchable indexes from reference sequences |
| **Compact Hash Table** | `compact_hash.{cc,h}` | 32-bit probabilistic hash table with linear probing collision resolution |
| **Minimizer Scanner** | `mmscanner.{cc,h}` | Extracts k-mer minimizers from sequences for indexing |
| **Taxonomy Management** | `taxonomy.{cc,h}` | NCBI taxonomy tree, parent-child relationships, LCA computation |
| **Sequence I/O** | `seqreader.{cc,h}` | FASTA/FASTQ input parsing |
| **Reporting** | `reports.{cc,h}` | Output generation (Kraken-style and MPA-style reports) |
| **Cardinality Estimation** | `hyperloglogplus.{cc,h}` | HyperLogLog+ cardinality estimation (from KrakenUniq) |

### Key Data Structures

```cpp
// namespace kraken2
struct IndexOptions              // Database metadata: k-mer size, l (minimizer length), seed masks
struct CompactHashCell           // 32-bit cell: high bits = truncated hash key, low bits = taxon ID
class CompactHashTable           // Compact probabilistic hash table with 256 lock zones
class Taxonomy                   // Flat-array taxonomy with parent/child pointers
class MinimizerScanner           // Sliding window minimizer extraction
class SequenceReader             // FASTA/FASTQ parsing (uses kseq.h)
```

### Directory Structure

```
src/
├── classify.cc              # Main classification entry point
├── build_db.cc              # Database construction
├── compact_hash.{cc,h}      # Hash table implementation
├── taxonomy.{cc,h}          # Taxonomy tree management
├── mmscanner.{cc,h}         # Minimizer extraction
├── seqreader.{cc,h}         # FASTA/FASTQ reading
├── reports.{cc,h}           # Output report generation
├── mmap_file.{cc,h}         # Memory-mapped file I/O
├── hyperloglogplus.{cc,h}   # Cardinality estimation
├── aa_translate.{cc,h}      # Amino acid translation for protein sequences
├── utilities.{cc,h}         # Common utilities
├── omp_hack.{cc,h}          # OpenMP workarounds
├── Makefile                 # Traditional build configuration
├── CMakeLists.txt           # CMake configuration
└── [other utilities]        # dump_table, estimate_capacity, k2mask, etc.

scripts/
├── kraken2                  # Main classification wrapper (Perl)
├── kraken2-build            # Database construction wrapper (Perl)
├── kraken2-inspect          # Database inspection tool (Perl)
├── k2                       # Unified interface (compiled binary)
└── [other scripts]          # Installation, helpers

data/
├── names.dmp                # NCBI taxonomy names
├── nodes.dmp                # NCBI taxonomy structure
└── *.fa                     # Reference sequences (COVID-19, Flu, HIV, etc.)
```

### Build Outputs

The build produces these executables in `src/`:
- `classify` - Core classification engine
- `build_db` - Database construction
- `estimate_capacity` - Pre-database size calculator
- `dump_table` - Inspect hash table contents
- `lookup_accession_numbers` - Query accession mappings
- `k2mask` - Low-complexity region masking
- `blast_to_fasta` - BLAST format converter
- `libtax.so` - Shared library for taxonomy access

## Algorithm Overview

### Minimizer-Based Indexing
Instead of storing all k-mers, Kraken 2 indexes only "minimizers" (smallest k-mers in sliding windows):
- **Default DNA parameters**: k=35 (k-mer size), l=31 (minimizer window)
- **Default protein parameters**: k=15, l=12
- This dramatically reduces memory footprint while maintaining classification accuracy

### Hash Table Design
- **Compact storage**: 32-bit cells with truncated hash in high bits, taxon ID in low bits
- **Linear probing**: O(1) average lookup time
- **Zone-based locking**: 256 lock zones for fine-grained parallelization

### Classification Workflow
1. Load database (taxonomy tree, hash table)
2. Read input sequence
3. Extract minimizers from sequence
4. Hash lookups in table to find candidate taxa
5. Compute Lowest Common Ancestor (LCA) from candidates
6. Generate taxonomic classification

### Taxonomic Storage
- Taxonomy stored as flat array for cache efficiency
- Names and ranks stored as continuous strings with offset pointers
- Parent pointers enable O(log n) LCA queries

## Important Implementation Details

### Parallelization Strategy
- **OpenMP** for multi-threaded classification and database building
- **Fragment batching**: Processes ~10,000 fragments per thread
- **Fine-grained locking**: 256 zones in hash table prevent contention

### Memory Mapping
Database files are memory-mapped for efficient access without full loading into memory. This enables processing of very large databases (hundreds of GB) with modest RAM requirements.

### Multi-Database Classification
Recent versions support merging results across multiple independent databases, allowing modular database composition.

### Quality Filtering
Optional base quality thresholds can be applied during classification for improved accuracy on low-quality reads.

## Build Dependencies

- **CMake** >= 2.8 (optional, traditional Makefile also available)
- **C++11 compiler** with OpenMP support (g++, clang++)
- **ZLIB** (for compression)
- **Pthreads** (POSIX threads)
- **Standard math library** (libm)

## References

- Kraken 2 paper: https://doi.org/10.1186/s13059-019-1891-0
- Original Kraken paper: https://dx.doi.org/10.1186/gb-2014-15-3-r46
- KrakenUniq paper: https://dx.doi.org/10.1186/s13059-018-1568-0
- See `docs/MANUAL.html` for complete user documentation

# Kraken 2 Rust Translation - Starting Point

Welcome! This directory contains a Rust translation of the Kraken 2 taxonomic sequence classification system. This document will get you oriented.

## Quick Overview

- **What is it?** A Rust rewrite of Kraken 2, a bioinformatics tool for classifying DNA sequences
- **Status**: Foundation (10 core modules) complete; ready for main algorithm implementation
- **Location**: `/Users/mladenrasic/Projects/kraken2/kraken2_rust/`
- **Original C++ code**: `../src/` directory

## Files You Should Read

1. **Start here**: This file (STARTING_POINT.md)
2. **Understanding the project**: [README.md](README.md)
3. **How to translate more code**: [CONVERSION_GUIDE.md](CONVERSION_GUIDE.md)
4. **What's done so far**: [TRANSLATION_SUMMARY.md](TRANSLATION_SUMMARY.md)

## The Layout

```
kraken2_rust/
├── Cargo.toml                 # Project dependencies
├── src/
│   ├── lib.rs                 # Module declarations
│   ├── utilities.rs           # ✅ DONE
│   ├── omp_hack.rs            # ✅ DONE
│   ├── mmap_file.rs           # ✅ DONE
│   ├── kraken2_data.rs        # ✅ DONE
│   ├── aa_translate.rs        # ✅ DONE
│   ├── seqreader.rs           # ✅ DONE
│   ├── mmscanner.rs           # ✅ DONE
│   ├── compact_hash.rs        # ✅ DONE
│   ├── taxonomy.rs            # ✅ DONE
│   └── reports.rs             # ✅ DONE
├── README.md                  # What each module does
├── CONVERSION_GUIDE.md        # Techniques & patterns
├── TRANSLATION_SUMMARY.md     # Summary of work done
└── STARTING_POINT.md          # This file
```

## What's Complete (10 Modules)

These modules are fully translated and tested:

| Module | Purpose | Status |
|--------|---------|--------|
| `utilities.rs` | String operations, bit manipulation | ✅ Complete |
| `omp_hack.rs` | Threading (uses Rayon) | ✅ Complete |
| `mmap_file.rs` | Memory-mapped file I/O | ✅ Complete |
| `kraken2_data.rs` | Data types & structures | ✅ Complete |
| `aa_translate.rs` | DNA → protein translation | ✅ Complete |
| `seqreader.rs` | Read FASTA/FASTQ files | ✅ Complete |
| `mmscanner.rs` | Extract k-mer minimizers | ✅ Complete |
| `compact_hash.rs` | Efficient hash table | ✅ Complete |
| `taxonomy.rs` | Taxonomy tree & LCA | ✅ Complete |
| `reports.rs` | Output formatting | ✅ Complete |

## What's Missing (Main Work)

These need to be translated next (in priority order):

1. **`build_db.rs`** - Builds the database (highest priority)
2. **`classify.rs`** - Main classification algorithm
3. Utility commands (`estimate_capacity`, `dump_table`, etc.)
4. CLI wrappers for commands

## Quick Start

### Build the project:
```bash
cd /Users/mladenrasic/Projects/kraken2/kraken2_rust
cargo build
```

### Run the tests:
```bash
cargo test
```

### Check for errors without building:
```bash
cargo check
```

## Key Differences from C++

| C++ | Rust | Why |
|-----|------|-----|
| `#pragma omp parallel` | `par_iter()` from rayon | Idiomatic Rust parallelization |
| `new`/`delete` | `Vec`, `Box` | Automatic memory management |
| Pointers `T*` | References `&T` | Memory safety |
| `std::map` | `BTreeMap` | Rust equivalent |
| Error codes | `Result<T>` | Better error handling |

## How to Add a New Module

Follow this pattern (see completed modules for examples):

1. **Create the file**: `src/module_name.rs`
2. **Add module declaration**: Add `pub mod module_name;` to `src/lib.rs`
3. **Write the code**: Translate from C++ following CONVERSION_GUIDE.md
4. **Add tests**: Include unit tests in the same file
5. **Test it**: `cargo test module_name`

Example structure:
```rust
//! Brief description of what this module does
//!
//! Translated from: module_name.cc/module_name.h

// Public API functions
pub fn my_function(arg: Type) -> Result<ReturnType> {
    // implementation
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_my_function() {
        let result = my_function(test_input);
        assert_eq!(result, expected);
    }
}
```

## Recommended Next Steps

1. **Read CONVERSION_GUIDE.md** - Learn the translation patterns
2. **Study a completed module** - Pick one (taxonomy.rs is good) and understand it
3. **Translate `build_db.rs`** - It's the highest priority remaining work
4. **Write comprehensive tests** - Each function should have tests
5. **Validate against C++** - Compare outputs with original implementation

## Common Commands

```bash
# Check for errors
cargo check

# Build debug version
cargo build

# Build optimized version
cargo build --release

# Run all tests
cargo test

# Run a specific test
cargo test test_name

# Run tests with output
cargo test -- --nocapture

# Show what was compiled
cargo build --verbose
```

## Key Files to Reference

- **Original C++ source**: `../src/build_db.cc`, `../src/classify.cc`
- **This translation**: `src/*.rs` files in this directory
- **Rust docs**: https://doc.rust-lang.org/book/
- **Rayon (parallelization)**: https://docs.rs/rayon/

## Understanding the Architecture

### The Pipeline

```
Input Sequences
      ↓
[seqreader.rs] - Read FASTA/FASTQ
      ↓
[mmscanner.rs] - Extract k-mer minimizers
      ↓
[compact_hash.rs] - Hash table lookups
      ↓
[taxonomy.rs] - LCA computation
      ↓
[reports.rs] - Format output
```

### The Data Flow

1. **Database Building** (`build_db.rs` - TODO)
   - Read sequences
   - Extract minimizers
   - Hash to taxon IDs
   - Build hash table

2. **Classification** (`classify.rs` - TODO)
   - Read query sequences
   - Extract minimizers
   - Look up in hash table
   - Compute LCA
   - Output results

## Staying Oriented

- **Confused about Rust?** Check the CONVERSION_GUIDE.md for patterns
- **Need to understand an algorithm?** Look at the C++ version in `../src/`
- **Want to see an example?** Find a similar completed module
- **Stuck on a module?** Read its tests - they show expected behavior

## Dependencies

The project uses these Rust crates (in Cargo.toml):

- **memmap2** - Memory-mapped files
- **rayon** - Parallel computation
- **anyhow** - Error handling
- **serde** - Serialization (future use)
- **flate2** - Compression support

All are well-maintained and widely used in Rust.

## Project Goals

✅ = Completed
⏳ = In Progress / Pending
❌ = Not Started

- ✅ Foundation modules (utilities, data structures)
- ✅ Core algorithms (taxonomy, hash table, translation)
- ✅ I/O (sequence reading, memory mapping)
- ⏳ Database building (`build_db.rs`)
- ⏳ Classification (`classify.rs`)
- ❌ Command-line interface
- ❌ Performance validation

## Getting Help

1. **Understand the C++ version**: Read `../src/filename.cc` for context
2. **See a pattern**: Look for a similar completed module
3. **Read documentation**: Check inline comments and doc strings
4. **Run tests**: `cargo test` to see what's expected
5. **Check CONVERSION_GUIDE.md**: Has patterns for common tasks

## Success Criteria

A module is complete when:

- [ ] All functions translated from C++
- [ ] All public functions have unit tests
- [ ] All tests pass (`cargo test`)
- [ ] Code compiles without warnings
- [ ] Documentation is complete
- [ ] Results match C++ implementation

## Notes

- **Rust version**: 2021 edition (latest stable)
- **Memory model**: Safe Rust with minimal unsafe blocks
- **Parallelization**: Rayon instead of OpenMP
- **Performance**: Target feature-parity with C++ version
- **Compatibility**: Maintain data structure layout for database compatibility

---

**You are here**: Foundation complete, ready to implement main algorithms

**Next**: Read [README.md](README.md) to understand what each module does, then [CONVERSION_GUIDE.md](CONVERSION_GUIDE.md) to learn translation patterns.

**Questions?** All documentation is in this directory.

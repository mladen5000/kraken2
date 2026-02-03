# Kraken 2 Rust Conversion Guide

This document provides guidance for continuing the Rust translation of Kraken 2.

## Translation Priority and Status

### ✅ Completed (7 modules)

1. **utilities.rs** - String operations, bit manipulation
2. **omp_hack.rs** - Threading abstraction (Rayon-based)
3. **mmap_file.rs** - Memory-mapped file access
4. **kraken2_data.rs** - Core data types and structures
5. **aa_translate.rs** - DNA → protein translation
6. **seqreader.rs** - FASTA/FASTQ file reading
7. **mmscanner.rs** - Minimizer extraction
8. **compact_hash.rs** - Space-efficient hash table
9. **taxonomy.rs** - Taxonomy tree management
10. **reports.rs** - Classification output formatting

### 🚧 In Progress / Blocked

None - all foundation modules completed.

### ⏳ Not Started (requires foundation)

1. **build_db.rs** - Database construction
   - Depends on: seqreader, mmscanner, compact_hash, taxonomy
   - ~500 lines of C++ code
   - Handles: sequence loading, minimizer extraction, hash table population

2. **classify.rs** - Main classification engine
   - Depends on: all core modules
   - ~900 lines of C++ code
   - Handles: sequence classification, LCA computation, result reporting

3. **hyperloglogplus.rs** - Advanced cardinality estimation
   - Depends on: nothing critical
   - ~600 lines of C++ code
   - Provides: minimizer-level abundance estimation

4. **lookup_accession_numbers.rs** - Accession mapping
   - Simple file-based lookup
   - ~100 lines of C++ code

5. **Command-line interfaces**
   - `cmd_classify.rs` - Wrapper around `classify()`
   - `cmd_build_db.rs` - Wrapper around `build_db()`
   - `cmd_inspect.rs` - Database inspection
   - ~200 lines total

## File-by-File C++ to Rust Mapping

### Source Files to Convert

```
src/
├── build_db.cc              → src/build_db.rs
├── classify.cc              → src/classify.rs
├── hyperloglogplus.cc       → src/hyperloglogplus.rs
├── lookup_accession_numbers.cc → src/lookup_accession_numbers.rs
├── dump_table.cc            → src/cmd_dump_table.rs
├── estimate_capacity.cc     → src/cmd_estimate_capacity.rs
├── k2mask.cc                → src/cmd_k2mask.rs
├── blast_to_fasta.c         → src/blast_to_fasta.rs
├── aa_translate.cc          ✅ src/aa_translate.rs
├── utilities.cc             ✅ src/utilities.rs
├── mmap_file.cc             ✅ src/mmap_file.rs
├── omp_hack.cc              ✅ src/omp_hack.rs
├── taxonomy.cc              ✅ src/taxonomy.rs
├── mmscanner.cc             ✅ src/mmscanner.rs
├── compact_hash.cc          ✅ src/compact_hash.rs
├── seqreader.cc             ✅ src/seqreader.rs
└── reports.cc               ✅ src/reports.rs
```

## Translation Techniques

### 1. C++ → Rust Data Structures

| C++ | Rust | Notes |
|-----|------|-------|
| `std::vector<T>` | `Vec<T>` | Growable array |
| `std::map<K,V>` | `BTreeMap<K,V>` | Ordered map |
| `std::unordered_map<K,V>` | `HashMap<K,V>` | Hash map |
| `std::string` | `String` or `&str` | Owned vs borrowed |
| `uint64_t` | `u64` | 64-bit unsigned |
| `omp_lock_t` | `Mutex<T>` | Thread synchronization |
| Pointers `T*` | References `&T` or `&mut T` | Borrow checking |

### 2. Memory Management

**C++:**
```cpp
void* buffer = malloc(size);
CompactHashCell* cells = new CompactHashCell[size];
// ... use ...
free(buffer);
delete[] cells;
```

**Rust:**
```rust
let mut buffer: Vec<u8> = vec![0; size];
let mut cells: Vec<CompactHashCell> = vec![CompactHashCell(0); size];
// ... use ...
// Automatically freed at end of scope
```

### 3. OpenMP Parallelization

**C++:**
```cpp
#pragma omp parallel for
for (int i = 0; i < n; i++) {
    process(data[i]);
}
```

**Rust:**
```rust
use rayon::prelude::*;

data.par_iter()
    .for_each(|item| {
        process(item);
    });
```

### 4. File I/O

**C++:**
```cpp
FILE* f = fopen(path, "rb");
fseek(f, offset, SEEK_SET);
fread(&value, sizeof(T), 1, f);
fclose(f);
```

**Rust:**
```rust
use memmap2::Mmap;

let file = File::open(path)?;
let mmap = unsafe { Mmap::map(&file)? };
let value: T = mmap.read_at(offset)?;
// Automatically unmapped
```

### 5. Error Handling

**C++:**
```cpp
void process_file(const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) {
        perror("Failed to open");
        return;
    }
    // ... process ...
    fclose(f);
}
```

**Rust:**
```rust
fn process_file(path: &str) -> Result<()> {
    let file = File::open(path)?;
    // ... process ...
    Ok(())
}
```

## Step-by-Step Conversion Process

### For Each Module

1. **Read the C++ header file** to understand the public API
2. **Identify the data structures** and translate to Rust
3. **Translate function signatures** to Rust idioms
4. **Implement core logic**, translating C++ algorithms to Rust
5. **Add tests** - copy C++ tests or write Rust equivalents
6. **Validate equivalence** - test with same inputs as C++ version

### Example: Converting `build_db.cc`

```rust
// 1. Read the header (build_db.h) - understand what it does
// 2. Create src/build_db.rs with module structure
// 3. Translate main function signature:
pub fn build_database(
    opts: &IndexOptions,
    sequences: &[Sequence],
    taxonomy: &Taxonomy,
) -> Result<CompactHashTable>

// 4. Translate the algorithm:
fn extract_and_hash_minimizers(...) { ... }
fn populate_hash_table(...) { ... }

// 5. Add tests:
#[test]
fn test_build_from_fasta() { ... }

#[test]
fn test_hash_table_population() { ... }
```

## Common Pitfalls

### 1. Mutable Borrowing

**Problem:**
```rust
fn modify_vector(v: &Vec<i32>) {
    v.push(5);  // ERROR: cannot borrow as mutable
}
```

**Solution:**
```rust
fn modify_vector(v: &mut Vec<i32>) {
    v.push(5);  // OK
}
```

### 2. String vs &str

**Problem:**
```rust
fn get_name(s: String) -> &str {
    &s  // ERROR: value dropped at end of scope
}
```

**Solution:**
```rust
fn get_name(s: &str) -> &str {
    s  // OK - borrowed reference
}
```

### 3. Unsafe Code

Minimize unsafe blocks. Use them only for:
- Memory-mapped file access
- FFI to C libraries
- Performance-critical pointer operations

Example:
```rust
unsafe {
    let ptr = data.as_ptr().add(offset) as *const T;
    Some(*ptr)
}
```

### 4. Lifetime Annotations

When returning references from function arguments:
```rust
fn get_node<'a>(nodes: &'a HashMap<Id, Node>, id: Id) -> Option<&'a Node> {
    nodes.get(&id)
}
```

## Testing Strategy

### Unit Tests (Per Module)

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_functionality() {
        // arrange
        let input = "test data";

        // act
        let result = process(input);

        // assert
        assert_eq!(result, expected_output);
    }
}
```

### Integration Tests

Create `tests/` directory for full-pipeline tests:

```rust
// tests/integration_test.rs
#[test]
fn test_full_pipeline() {
    let seqs = load_test_sequences();
    let db = build_database(&seqs);
    let results = classify(&db, &test_queries);
    assert_eq!(results.len(), test_queries.len());
}
```

### Cargo Test Commands

```bash
# Run all tests
cargo test

# Run specific test
cargo test test_name

# Run with output
cargo test -- --nocapture

# Run tests in release mode (faster)
cargo test --release

# Show test discovery without running
cargo test --no-run
```

## Performance Optimization

### Before Optimizing

1. Ensure correctness first (pass all tests)
2. Profile with `cargo flamegraph` or `perf`
3. Identify actual bottlenecks

### Common Optimizations

| Optimization | Example |
|--------------|---------|
| Avoid cloning | Use references `&T` instead of `T` |
| Use iterators | `vec.iter().map()` instead of for loops |
| Parallel iteration | `vec.par_iter()` from rayon |
| Memory pooling | Reuse allocations between iterations |
| SIMD | `packed_simd` for bulk operations |

### Release Builds

```bash
# Optimized binary
cargo build --release

# With further optimizations
cargo build --release -C opt-level=3 -C lto=fat -C codegen-units=1
```

## Documentation

### Code Comments

```rust
/// Brief description of the function.
///
/// Longer description with examples and behavior details.
///
/// # Arguments
/// * `param1` - Description of first parameter
/// * `param2` - Description of second parameter
///
/// # Returns
/// Description of return value
///
/// # Examples
/// ```
/// let result = function(arg1, arg2);
/// assert_eq!(result, expected);
/// ```
pub fn function(param1: Type1, param2: Type2) -> ReturnType {
    // implementation
}
```

### Module-level Documentation

```rust
//! # Module Name
//!
//! Brief description of what this module does.
//!
//! ## Features
//! - Feature 1
//! - Feature 2
//!
//! ## Translated from
//! C++ files: filename.cc, filename.h
```

## Next Steps

1. **Start with `build_db.rs`** - Core algorithm, but depends on completed modules
2. **Then `classify.rs`** - Main engine, orchestrates other modules
3. **Then utility commands** - Smaller, fewer dependencies
4. **Finally CLI** - Once core algorithms are working

## Useful Rust Resources

- [Rust Book](https://doc.rust-lang.org/book/) - Official guide
- [Rayon Documentation](https://docs.rs/rayon/) - Parallel iteration
- [Memmap2](https://docs.rs/memmap2/) - Memory mapping
- [Cargo Book](https://doc.rust-lang.org/cargo/) - Build system

## Git Integration

```bash
# Track conversion progress
git add kraken2_rust/src/*.rs
git commit -m "Add module_name translation

Translates module_name.cc/h from C++ to Rust.

- Implements all public functions
- Adds unit tests
- Maintains API compatibility"
```

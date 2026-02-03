# Kraken 2 Rust - Usage Examples

## Overview

The Kraken 2 Rust translation provides a clean API for building and using taxonomic sequence classification databases. This document shows how to use the main components.

## Building a Database

### Basic Database Construction

```rust
use kraken2_rust::api::DatabaseBuilder;

fn main() -> anyhow::Result<()> {
    let builder = DatabaseBuilder::new()
        .with_k(35)
        .with_l(31)
        .with_id_map("seqid2taxid.txt".to_string())
        .with_taxonomy("/path/to/taxonomy".to_string())
        .with_output("kraken_db".to_string());

    builder.build("sequences.fasta")?;
    Ok(())
}
```

### DNA vs Protein Databases

```rust
use kraken2_rust::api::DatabaseBuilder;

// DNA database (default)
let dna_builder = DatabaseBuilder::new()
    .with_k(35)
    .with_l(31);

// Protein database
let protein_builder = DatabaseBuilder::new()
    .with_protein(true)
    .with_k(15)
    .with_l(12);
```

### Parallel Processing

```rust
use kraken2_rust::api::DatabaseBuilder;

let builder = DatabaseBuilder::new()
    .with_threads(4)  // Use 4 threads for processing
    .with_k(35)
    .with_l(31);
```

## Classifying Sequences

### Basic Classification

```rust
use kraken2_rust::api::Classifier;

fn main() -> anyhow::Result<()> {
    let classifier = Classifier::new()
        .with_database("kraken_db".to_string())
        .with_report_file("classifications.txt".to_string());

    // The classifier is ready to use
    Ok(())
}
```

### Classification with Confidence Threshold

```rust
use kraken2_rust::api::Classifier;

let classifier = Classifier::new()
    .with_database("kraken_db".to_string())
    .with_confidence_threshold(0.8);  // 80% confidence threshold
```

### Quick Mode (Faster but Less Accurate)

```rust
use kraken2_rust::api::Classifier;

let classifier = Classifier::new()
    .with_database("kraken_db".to_string())
    .with_quick_mode(true);  // Enable quick mode
```

### MPA-Style Output

```rust
use kraken2_rust::api::Classifier;

let classifier = Classifier::new()
    .with_database("kraken_db".to_string())
    .with_mpa_format(true);  // Use MetaPhlan-style format
```

## Core Modules

### Utilities Module

String manipulation and bit operations:

```rust
use kraken2_rust::utilities::{split_string, expand_spaced_seed_mask};

// Split strings by delimiter
let parts = split_string("part1\tpart2\tpart3", "\t", None);
assert_eq!(parts.len(), 3);

// Split with max fields
let parts = split_string("a\tb\tc\td", "\t", Some(2));
assert_eq!(parts.len(), 2);
```

### Sequence Reader

Reading FASTA/FASTQ files:

```rust
use kraken2_rust::seqreader::SequenceReader;

fn main() -> anyhow::Result<()> {
    let mut reader = SequenceReader::new("sequences.fasta")?;

    while let Some(seq) = reader.next_sequence()? {
        println!("Sequence: {}", seq.header);
        println!("Length: {}", seq.seq.len());
        if let Some(quality) = &seq.quality {
            println!("Quality: {}", quality);
        }
    }
    Ok(())
}
```

### Taxonomy Management

Working with NCBI taxonomy:

```rust
use kraken2_rust::taxonomy::Taxonomy;

fn main() -> anyhow::Result<()> {
    let taxonomy = Taxonomy::load_from_ncbi(
        "nodes.dmp",
        "names.dmp"
    )?;

    // Get node information
    if let Some(name) = taxonomy.get_name(562) {
        println!("Taxon 562: {}", name);
    }

    // Compute LCA
    let lca = taxonomy.lca(562, 561);
    println!("LCA of 562 and 561: {}", lca);

    println!("Total taxa: {}", taxonomy.size());
    Ok(())
}
```

### Minimizer Extraction

Extracting minimizers from sequences:

```rust
use kraken2_rust::mmscanner::MinimizerScanner;
use kraken2_rust::kraken2_data::IndexOptions;

fn main() {
    let opts = IndexOptions {
        k: 35,
        l: 31,
        dna_db: true,
        ..Default::default()
    };

    let scanner = MinimizerScanner::new(&opts);
    let minimizers = scanner.scan("ATCGATCGATCGATCGATCGATCGATCGATCG");

    println!("Found {} minimizers", minimizers.len());
}
```

### Hash Table Operations

Compact hash table for k-mer storage:

```rust
use kraken2_rust::compact_hash::CompactHashTable;

fn main() {
    let table = CompactHashTable::new(10000);

    // Insert a minimizer -> taxon mapping
    table.insert(0x123456789ABCDEF0, 562)?;

    // Lookup
    if let Some(taxid) = table.lookup(0x123456789ABCDEF0) {
        println!("Taxon ID: {}", taxid);
    }

    println!("Capacity: {}", table.capacity());
}
```

## Advanced Usage

### Custom Builder Configuration

```rust
use kraken2_rust::api::DatabaseBuilder;

let mut builder = DatabaseBuilder::new();

// Access and modify options directly for advanced usage
{
    let opts = builder.options_mut();
    opts.k = 35;
    opts.l = 31;
    opts.spaced_seed_mask = 0x1000000000000000;
    opts.minimum_hit_groups = 2;
}

builder.build("sequences.fasta")?;
```

### Custom Classifier Configuration

```rust
use kraken2_rust::api::Classifier;

let mut classifier = Classifier::new();

// Access and modify options directly
{
    let opts = classifier.options_mut();
    opts.minimum_quality_score = 20;
    opts.minimum_hit_groups = 3;
    opts.quick_mode = false;
}
```

## Command-Line Usage

### Building a Database

```bash
# Using the Rust Kraken 2 tools directly
cargo run --release --bin kraken2-build -- \
    --threads 4 \
    --db-type standard \
    --library-dir ~/bacteria \
    ~/bacteria/my_kraken2_db
```

### Classifying Sequences

```bash
# Classify sequences
cargo run --release --bin kraken2 -- \
    --db ~/bacteria/my_kraken2_db \
    --threads 4 \
    --confidence 0.8 \
    --output output.txt \
    sequences.fasta
```

## API Design Principles

1. **Builder Pattern**: Use `DatabaseBuilder` and `Classifier` for configuration
2. **Fluent Interface**: Chain method calls for readable configuration
3. **Type Safety**: Rust's type system prevents invalid configurations at compile time
4. **Error Handling**: All fallible operations return `Result<T, E>`
5. **Memory Safety**: No buffer overflows or use-after-free errors possible

## Performance Considerations

- **Parallelization**: Use `with_threads()` for multi-threaded processing
- **Quick Mode**: Enable for faster (but less accurate) classification
- **Memory Mapping**: Automatically used for efficient file I/O
- **Minimizer Caching**: Minimizers are extracted once and reused

## Error Handling

```rust
use kraken2_rust::api::DatabaseBuilder;

fn main() {
    match DatabaseBuilder::new().build("nonexistent.fasta") {
        Ok(_) => println!("Success"),
        Err(e) => eprintln!("Error: {}", e),
    }
}
```

## Testing

Run the comprehensive test suite:

```bash
# Run all tests
cargo test --lib

# Run integration tests
cargo test --test integration_test

# Run with output
cargo test -- --nocapture
```

## Feature Flags (Future)

```toml
[features]
default = []
parallel = ["rayon"]
memory-mapping = ["memmap2"]
```

## Contributing

To add new functionality:

1. Add the module to `src/`
2. Export it in `lib.rs`
3. Add tests in the module's `#[cfg(test)]` section
4. Add integration tests in `tests/`
5. Update this documentation

## References

- [Kraken 2 Paper](https://doi.org/10.1186/s13059-019-1891-0)
- [Rust Book](https://doc.rust-lang.org/book/)
- [Cargo Documentation](https://doc.rust-lang.org/cargo/)

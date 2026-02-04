# HyperLogLog++ Implementation Notes

## Overview

This module implements the HyperLogLog++ algorithm for cardinality estimation, originally developed by Florian Breitwieser for the KrakenUniq taxonomic classification system.

The Rust implementation is a complete, idiomatic translation of the C++ code from `hyperloglogplus.cc` and `hyperloglogplus.h`.

## Key Features

### 1. Sparse Representation
- For low cardinalities, uses a `BTreeSet<u32>` to store encoded hash values
- Automatically switches to normal representation when the sparse list grows beyond m/4 elements
- Provides better accuracy for small cardinalities (< 1000 unique elements)

### 2. Three Cardinality Estimators

#### Ertl's Improved Estimator (Recommended)
```rust
let estimate = hll.ertl_cardinality();
```
- No empirical bias correction needed
- Good accuracy across all cardinality ranges
- Default estimator used by `cardinality()`

#### Heule's HLL++ Estimator
```rust
let estimate = hll.heule_cardinality(true); // with bias correction
```
- Uses empirical bias correction tables
- Switches between linear counting and HLL estimation
- Limited to precision ≤ 18

#### Flajolet's Original HLL Estimator
```rust
let estimate = hll.flajolet_cardinality(true); // use sparse precision
```
- The original HyperLogLog algorithm
- Uses linear counting for small cardinalities

### 3. Hash Functions

**MurmurHash3 Finalizer** (default):
```rust
use kraken2_rust::hyperloglogplus::murmurhash3_finalizer;
let hash = murmurhash3_finalizer(value);
```

**Wang Mixer** (alternative):
```rust
use kraken2_rust::hyperloglogplus::wang_mixer;
let hll = HyperLogLogPlusMinus::with_hash(12, true, wang_mixer);
```

## Usage Examples

### Basic Usage
```rust
use kraken2_rust::hyperloglogplus::HyperLogLogPlusMinus;

// Create with precision 12 (4096 registers), sparse mode enabled
let mut hll = HyperLogLogPlusMinus::new(12, true);

// Insert values
for i in 0..10000 {
    hll.insert(i);
}

// Get cardinality estimate
let estimate = hll.cardinality();
println!("Estimated unique count: {}", estimate);
println!("Actual count: 10000");
println!("Error: {}%", ((estimate as f64 - 10000.0) / 10000.0).abs() * 100.0);
```

### Merging Sketches
```rust
let mut hll1 = HyperLogLogPlusMinus::new(12, true);
let mut hll2 = HyperLogLogPlusMinus::new(12, true);

// Add data to both sketches
for i in 0..1000 {
    hll1.insert(i);
}
for i in 500..1500 {
    hll2.insert(i);
}

// Merge hll2 into hll1
hll1.merge(&hll2);

// Or use += operator
hll1 += &hll2;

// Estimate should be ~1500 (union of both sets)
let estimate = hll1.cardinality();
```

### Batch Insertion
```rust
let mut hll = HyperLogLogPlusMinus::new(14, true);

let data: Vec<u64> = (0..100000).collect();
hll.insert_slice(&data);

println!("Cardinality: {}", hll.cardinality());
```

### Precision vs Memory Tradeoff
```rust
// Low precision (16 registers) - less memory, higher error
let hll_p4 = HyperLogLogPlusMinus::new(4, true);  // ~1% standard error

// Medium precision (4096 registers) - balanced
let hll_p12 = HyperLogLogPlusMinus::new(12, true); // ~1.625% standard error

// High precision (262144 registers) - more memory, lower error
let hll_p18 = HyperLogLogPlusMinus::new(18, true); // ~0.406% standard error
```

## Algorithm Details

### Encoding Scheme (Sparse Mode)
In sparse mode, 64-bit hash values are encoded as 32-bit integers:
- First `p'` bits (25 bits) are used as the index
- If bits after position `p` are all zero, additional rank information is stored
- LSB indicates whether rank information is present

### Register Update (Normal Mode)
Each register stores the maximum rank observed for its index:
```
rank = position of first 1-bit after the index bits + 1
```

### Cardinality Estimation (Ertl)
```
estimate = (m² × α∞) / (m×τ(...) + Σ(c_k × 2^(-k)) + m×σ(...))
```
where:
- `m` = number of registers
- `α∞` = m / (2 × ln(2))
- `c_k` = register histogram
- `σ` and `τ` = correction functions

## Differences from C++ Implementation

### Type Choices
- `BTreeSet<u32>` instead of `set<uint32_t>` for sparse list (provides ordering)
- `fn(u64) -> u64` instead of function pointer with templates
- Native Rust integer types (u8, u32, u64) with explicit wrapping arithmetic

### Memory Management
- No manual memory management (Rust handles allocation/deallocation)
- Clone trait for deep copying instead of copy constructors
- Move semantics built into the language

### Error Handling
- Assertions for precondition checking (panic on invalid input)
- No exceptions (Rust doesn't have exceptions)

### Traits Implemented
- `Default` - creates HLL with precision 12, sparse mode
- `Clone` - deep copy of the entire structure
- `AddAssign<&HyperLogLogPlusMinus>` - merge operator (+=)

## Performance Characteristics

### Time Complexity
- **Insert**: O(log m) for sparse mode (BTreeSet), O(1) for normal mode
- **Cardinality**: O(m) where m is the number of registers
- **Merge**: O(m) for normal mode, O(s₁ + s₂) for sparse mode

### Space Complexity
- **Sparse mode**: O(unique elements) up to m/4
- **Normal mode**: O(2^p) bytes (one byte per register)

### Accuracy
Standard error for different precisions:
- p=4:  ~26%
- p=8:  ~6.5%
- p=12: ~1.625%
- p=14: ~0.8125%
- p=18: ~0.406%

## Testing

The module includes comprehensive tests covering:
- Basic insertion and cardinality estimation
- Sparse to normal mode conversion
- Merge operations (sparse and normal)
- Different estimator algorithms
- Hash function determinism
- Edge cases (very small/large cardinalities)
- Duplicate handling
- Reset functionality

Run tests with:
```bash
cargo test hyperloglogplus --lib
```

## References

1. Flajolet, P., et al. "HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm." (2007)
2. Heule, S., Nunkesser, M., & Hall, A. "HyperLogLog in practice: algorithmic engineering of a state of the art cardinality estimation algorithm." (2013)
3. Ertl, O. "New cardinality estimation algorithms for HyperLogLog sketches." (2017)

## Notes on Bias Correction Data

The current implementation includes truncated bias correction tables. For production use with the Heule estimator at precisions 5-18, you should populate the full bias data arrays from the C++ `hyperloglogplus-bias.h` file.

The Ertl estimator (default) does not require these tables and is recommended for most use cases.

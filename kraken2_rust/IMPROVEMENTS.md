# Potential Improvements for Kraken 2 Rust Implementation

Based on recent research from arXiv and related publications, here are concrete suggestions for making the implementation more lightweight, faster, accurate, or robust.

---

## 1. Performance: IDL Hash Function (High Priority)

**Source**: [arXiv:2406.14901](https://arxiv.org/abs/2406.14901) - "IDL: A Locality-Sensitive Hash Function for Efficient Genome Search"

**Current state**: `compact_hash.rs` uses standard hash indexing with linear probing.

**Improvement**: The IDL (Interleaved Double Lookups) hash function provides:
- 5x reduction in cache misses
- 2x query speedup
- Drop-in replacement compatibility

**Implementation sketch**:
```rust
/// IDL hash function - interleaves bits from two lookups for better locality
pub fn idl_hash(kmer: u64, table_size: usize) -> (usize, usize) {
    let h1 = murmur_hash(kmer);
    let h2 = murmur_hash(kmer.rotate_left(32));

    // Interleave bits for locality
    let interleaved = interleave_bits(h1, h2);
    let primary = (interleaved as usize) % table_size;
    let secondary = ((interleaved >> 32) as usize) % table_size;

    (primary, secondary)
}

fn interleave_bits(a: u64, b: u64) -> u64 {
    let mut result = 0u64;
    for i in 0..32 {
        result |= ((a >> i) & 1) << (2 * i);
        result |= ((b >> i) & 1) << (2 * i + 1);
    }
    result
}
```

**Files to modify**: `compact_hash.rs`

---

## 2. Memory: Simplified HyperLogLog (Medium Priority)

**Source**: Recent analysis shows simpler estimators match HyperLogLog++ accuracy

**Current state**: `hyperloglogplus.rs` implements full HLL++ with:
- Sparse/normal mode switching
- Bias correction tables (723 entries)
- Three different estimators

**Improvement options**:

### Option A: Remove bias correction tables
The Ertl MLE estimator doesn't need bias tables and provides similar accuracy:
```rust
/// Simplified cardinality estimate using Ertl's MLE
pub fn estimate_cardinality_simple(&self) -> f64 {
    let m = self.num_registers as f64;
    let registers = self.get_registers();

    // Count register values
    let mut counts = [0u32; 65];
    for &reg in &registers {
        counts[reg as usize] += 1;
    }

    // Ertl's MLE estimate
    let sum: f64 = (0..=64).map(|i| counts[i] as f64 * 2.0_f64.powi(-(i as i32))).sum();
    m * m / sum * 0.7213 / (1.0 + 1.079 / m)
}
```

### Option B: Use HyperLogLogLog (~40% less memory)
For very large cardinalities, HyperLogLogLog stores register differences instead of absolute values.

**Files to modify**: `hyperloglogplus.rs`

---

## 3. Accuracy: KATKA Kernels for Better Specificity (Medium Priority)

**Source**: [arXiv:2402.06935](https://arxiv.org/abs/2402.06935) - "KATKA: Exploiting Maximal Exact Matches"

**Current state**: Classification uses minimizer-based k-mer matching.

**Improvement**: KATKA kernels use maximal exact matches (MEMs) instead of k-mers:
- Better handling of sequencing errors
- More specific at strain level
- Can be combined with existing minimizer index

**Implementation approach**:
```rust
/// Find maximal exact matches between query and reference
pub struct KatkaKernel {
    /// Suffix array for MEM finding
    suffix_array: Vec<u32>,
    /// LCP array
    lcp_array: Vec<u32>,
}

impl KatkaKernel {
    /// Find all MEMs of length >= min_len
    pub fn find_mems(&self, query: &[u8], min_len: usize) -> Vec<MEM> {
        // Use suffix array + LCP for O(n) MEM finding
        todo!()
    }

    /// Score based on MEM coverage rather than k-mer count
    pub fn score_classification(&self, mems: &[MEM], taxonomy: &Taxonomy) -> TaxId {
        // Weight by MEM length and coverage
        todo!()
    }
}
```

**Files to modify**: New module `katka.rs`, integrate with `classify.rs`

---

## 4. Accuracy: Minimizer Digests for Compression (Low Priority)

**Source**: Same KATKA paper describes "minimizer digests"

**Current state**: Full minimizer values stored in hash table.

**Improvement**: Store only a digest (fingerprint) of minimizers:
- 50-70% space reduction
- Slightly higher false positive rate
- Useful for very large databases

```rust
/// Compressed minimizer representation
pub struct MinimizerDigest {
    /// Bloom filter for quick membership test
    bloom: BloomFilter,
    /// Compact storage for top minimizers only
    top_minimizers: Vec<u32>,
}
```

---

## 5. Robustness: ProbMinHash for ANI Estimation (Low Priority)

**Source**: Recent MinHash variants provide tighter variance bounds

**Current state**: No direct ANI estimation in classification.

**Improvement**: Add optional ANI confidence scoring:
```rust
/// ProbMinHash sketch for ANI estimation
pub struct ProbMinHashSketch {
    registers: Vec<f64>,
    num_registers: usize,
}

impl ProbMinHashSketch {
    /// Estimate ANI between two sketches
    pub fn estimate_ani(&self, other: &Self) -> f64 {
        let intersection = self.estimate_intersection(other);
        let union = self.estimate_union(other);
        intersection / union
    }
}
```

This could provide confidence scores for classifications, flagging low-ANI matches.

---

## 6. Strain-Level: Pangenome Integration (Future Work)

**Source**: PanTax and similar tools

For true strain-level classification, consider:
1. Build pangenome graph per species
2. Map reads to graph variants
3. Score based on variant coverage

This is a significant architectural change but provides the best strain resolution.

---

## Implementation Priority (Ordered by Complexity)

Below are all improvements ordered from simplest to most complex, with detailed complexity estimates.

### Complexity Level 1: Simple (~1-2 hours, single file changes)

| # | Improvement | Impact | Files | Lines of Code | Risk |
|---|-------------|--------|-------|---------------|------|
| 1 | **IDL hash function** | 2x query speed | `compact_hash.rs` | ~50 | Low |
| 2 | **Remove HLL++ bias tables** | Simpler code, less memory | `hyperloglogplus.rs` | ~30 (remove ~100) | Low |

**IDL Hash Details:**
- Add `idl_hash()` function (~20 lines)
- Add `interleave_bits()` helper (~15 lines)
- Modify `insert()` and `lookup()` to use new hash (~10 lines)
- No API changes, no new dependencies

**Simplified HLL Details:**
- Remove `BIAS_DATA` and `RAW_ESTIMATE_DATA` arrays (~700 lines removed)
- Replace `estimate_bias()` with Ertl's formula (~10 lines)
- Keep sparse mode for small cardinalities

---

### Complexity Level 2: Moderate (~4-8 hours, multiple files)

| # | Improvement | Impact | Files | Lines of Code | Risk |
|---|-------------|--------|-------|---------------|------|
| 3 | **Minimizer digests** | 50-70% space reduction | `mmscanner.rs`, `compact_hash.rs` | ~150 | Medium |
| 4 | **ProbMinHash ANI** | Confidence scores | New `probminhash.rs`, `classify.rs` | ~200 | Low |

**Minimizer Digests Details:**
- New `MinimizerDigest` struct with Bloom filter
- Modify `MinimizerScanner::scan()` to produce digests
- Add digest-based lookup in `CompactHashTable`
- Optional feature flag for backwards compatibility

**ProbMinHash ANI Details:**
- New module `probminhash.rs` (~150 lines)
- Integration hook in classification output (~50 lines)
- Optional feature, no changes to core algorithm

---

### Complexity Level 3: Significant (~2-4 days, new module + integration)

| # | Improvement | Impact | Files | Lines of Code | Risk |
|---|-------------|--------|-------|---------------|------|
| 5 | **KATKA kernels (MEM-based)** | Better strain accuracy | New `katka.rs`, `classify.rs` | ~400-600 | Medium |

**KATKA Kernels Details:**
- New module with suffix array construction (~200 lines)
- MEM finding algorithm (~150 lines)
- Scoring integration with existing LCA logic (~100 lines)
- Database builder changes for MEM index (~100 lines)
- Requires understanding of FM-index or suffix arrays

---

### Complexity Level 4: Major (~1-2 weeks, architectural change)

| # | Improvement | Impact | Files | Lines of Code | Risk |
|---|-------------|--------|-------|---------------|------|
| 6 | **Pangenome graph integration** | True strain resolution | Many files, new modules | ~2000+ | High |

**Pangenome Details:**
- New graph data structure for pangenomes
- Modified database builder for graph construction
- New classification algorithm for graph traversal
- Significant memory model changes
- Would benefit from external library (e.g., `petgraph`)

---

## Summary Table (All Improvements by Complexity)

| Rank | Improvement | Complexity | Time Est. | Impact | Dependencies |
|------|-------------|------------|-----------|--------|--------------|
| 1 | IDL hash function | ⭐ Simple | 1-2 hrs | High | None |
| 2 | Remove HLL++ bias tables | ⭐ Simple | 1-2 hrs | Medium | None |
| 3 | Minimizer digests | ⭐⭐ Moderate | 4-8 hrs | Medium | Bloom filter crate |
| 4 | ProbMinHash ANI | ⭐⭐ Moderate | 4-8 hrs | Low | None |
| 5 | KATKA kernels | ⭐⭐⭐ Significant | 2-4 days | High | Suffix array knowledge |
| 6 | Pangenome graphs | ⭐⭐⭐⭐ Major | 1-2 weeks | High | petgraph, research |

---

## Recommended Implementation Order

1. **Start here → IDL hash function** (biggest ROI: 2x speed for 1-2 hours work)
2. **Then → Remove HLL++ bias tables** (simplifies codebase)
3. **If space matters → Minimizer digests** (database size reduction)
4. **If accuracy matters → KATKA kernels** (requires more investment)
5. **Future → Pangenome graphs** (major project, consider as separate milestone)

---

## Quick Wins

1. **IDL hash in compact_hash.rs** - Immediate 2x speedup potential
2. **Remove bias tables from HLL++** - Simpler code, same accuracy
3. **Add optional MEM-based scoring** - Better strain specificity

## References

- IDL Hash: https://arxiv.org/abs/2406.14901
- KATKA: https://arxiv.org/abs/2402.06935
- HyperLogLogLog: https://arxiv.org/abs/2205.00564
- ProbMinHash: https://arxiv.org/abs/1911.00675
- PanTax: https://doi.org/10.1093/bioinformatics/btae028

//! HyperLogLog++ demonstration
//! 
//! This example demonstrates the basic usage of the HyperLogLog++ 
//! probabilistic cardinality estimation algorithm.

use kraken2_rust::hyperloglogplus::HyperLogLogPlusMinus;

fn main() {
    println!("=== HyperLogLog++ Cardinality Estimation Demo ===\n");

    // Example 1: Basic usage
    println!("Example 1: Basic cardinality estimation");
    let mut hll = HyperLogLogPlusMinus::new(12, true);
    
    let n_items = 10000;
    for i in 0..n_items {
        hll.insert(i);
    }
    
    let estimate = hll.cardinality();
    let error = ((estimate as f64 - n_items as f64) / n_items as f64).abs() * 100.0;
    println!("  Actual unique items: {}", n_items);
    println!("  Estimated count: {}", estimate);
    println!("  Error: {:.2}%\n", error);

    // Example 2: Handling duplicates
    println!("Example 2: Duplicate detection");
    let mut hll2 = HyperLogLogPlusMinus::new(12, true);
    
    // Insert 100 unique values, each 10 times
    for i in 0..100 {
        for _ in 0..10 {
            hll2.insert(i);
        }
    }
    
    println!("  Total insertions: {}", hll2.n_observed());
    println!("  Unique items: ~100");
    println!("  Estimated unique: {}\n", hll2.cardinality());

    // Example 3: Merging sketches
    println!("Example 3: Merging HyperLogLog sketches");
    let mut hll3a = HyperLogLogPlusMinus::new(14, true);
    let mut hll3b = HyperLogLogPlusMinus::new(14, true);
    
    // Set A: 0..1000
    for i in 0..1000 {
        hll3a.insert(i);
    }
    
    // Set B: 500..1500
    for i in 500..1500 {
        hll3b.insert(i);
    }
    
    println!("  Set A cardinality: {}", hll3a.cardinality());
    println!("  Set B cardinality: {}", hll3b.cardinality());
    
    hll3a.merge(&hll3b);
    println!("  Union (A ∪ B) cardinality: {}", hll3a.cardinality());
    println!("  Expected: ~1500\n");

    // Example 4: Different precisions
    println!("Example 4: Precision vs. accuracy trade-off");
    for precision in [4, 8, 12, 14, 18] {
        let mut hll_p = HyperLogLogPlusMinus::new(precision, true);
        let n = 10000;
        
        for i in 0..n {
            hll_p.insert(i);
        }
        
        let est = hll_p.cardinality();
        let err = ((est as f64 - n as f64) / n as f64).abs() * 100.0;
        let mem = 1 << precision; // bytes for normal mode
        
        println!("  p={:2}: {} registers, estimate={}, error={:.2}%", 
                 precision, mem, est, err);
    }
    
    println!("\n=== Demo Complete ===");
}

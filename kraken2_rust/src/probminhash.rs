/// ProbMinHash for ANI (Average Nucleotide Identity) Estimation
///
/// This module provides probabilistic similarity estimation between sequences
/// using the ProbMinHash algorithm. It can be used to add confidence scores
/// to taxonomic classifications.
///
/// Reference: arXiv:1911.00675 - "ProbMinHash: A Class of Locality-Sensitive
/// Hash Algorithms for the (Probability) Jaccard Similarity"

use std::cmp::Ordering;

/// Number of registers in the sketch (determines accuracy vs memory tradeoff)
const DEFAULT_NUM_REGISTERS: usize = 128;

/// ProbMinHash sketch for ANI estimation
///
/// This sketch stores the minimum hash values across multiple hash functions,
/// allowing efficient estimation of Jaccard similarity between sets.
#[derive(Clone, Debug)]
pub struct ProbMinHashSketch {
    /// Minimum hash values for each register
    registers: Vec<f64>,
    /// Number of registers
    num_registers: usize,
    /// Number of items added
    count: u64,
}

impl Default for ProbMinHashSketch {
    fn default() -> Self {
        Self::new(DEFAULT_NUM_REGISTERS)
    }
}

impl ProbMinHashSketch {
    /// Create a new ProbMinHash sketch with the specified number of registers
    ///
    /// # Arguments
    /// * `num_registers` - Number of hash registers (more = better accuracy, more memory)
    pub fn new(num_registers: usize) -> Self {
        ProbMinHashSketch {
            registers: vec![f64::INFINITY; num_registers],
            num_registers,
            count: 0,
        }
    }

    /// Insert a k-mer hash into the sketch
    ///
    /// # Arguments
    /// * `hash` - The hash value to insert
    pub fn insert(&mut self, hash: u64) {
        self.count += 1;

        for i in 0..self.num_registers {
            // Generate register-specific hash using a simple mixing function
            let mixed = hash_with_seed(hash, i as u64);
            // Convert to uniform [0, 1) distribution
            let value = hash_to_uniform(mixed);

            if value < self.registers[i] {
                self.registers[i] = value;
            }
        }
    }

    /// Insert multiple k-mer hashes
    pub fn insert_slice(&mut self, hashes: &[u64]) {
        for &hash in hashes {
            self.insert(hash);
        }
    }

    /// Estimate Jaccard similarity between this sketch and another
    ///
    /// # Arguments
    /// * `other` - The other sketch to compare against
    ///
    /// # Returns
    /// Estimated Jaccard similarity in range [0, 1]
    pub fn estimate_jaccard(&self, other: &Self) -> f64 {
        assert_eq!(
            self.num_registers, other.num_registers,
            "Sketches must have same number of registers"
        );

        let mut matches = 0;

        for i in 0..self.num_registers {
            // Check if minimum values are approximately equal
            // (indicating they likely came from the same element)
            if (self.registers[i] - other.registers[i]).abs() < 1e-10 {
                matches += 1;
            }
        }

        matches as f64 / self.num_registers as f64
    }

    /// Estimate ANI (Average Nucleotide Identity) from Jaccard similarity
    ///
    /// Uses the Mash distance formula: D = -1/k * ln(2J / (1+J))
    /// Then converts to ANI: ANI = 1 - D
    ///
    /// # Arguments
    /// * `other` - The other sketch to compare against
    /// * `k` - The k-mer size used
    ///
    /// # Returns
    /// Estimated ANI in range [0, 1], or None if Jaccard is too low
    pub fn estimate_ani(&self, other: &Self, k: usize) -> Option<f64> {
        let jaccard = self.estimate_jaccard(other);

        if jaccard <= 0.0 {
            return None;
        }

        // Mash distance formula
        let two_j = 2.0 * jaccard;
        let one_plus_j = 1.0 + jaccard;

        if two_j / one_plus_j <= 0.0 {
            return None;
        }

        let distance = -1.0 / (k as f64) * (two_j / one_plus_j).ln();

        // Convert distance to identity
        let ani = 1.0 - distance;

        // Clamp to valid range
        Some(ani.clamp(0.0, 1.0))
    }

    /// Estimate containment of self in other (|A ∩ B| / |A|)
    ///
    /// This is useful for determining how much of a query sequence
    /// is covered by a reference.
    pub fn estimate_containment(&self, other: &Self) -> f64 {
        assert_eq!(
            self.num_registers, other.num_registers,
            "Sketches must have same number of registers"
        );

        let mut contained = 0;

        for i in 0..self.num_registers {
            // Element is contained if other's min is <= self's min
            if other.registers[i] <= self.registers[i] + 1e-10 {
                contained += 1;
            }
        }

        contained as f64 / self.num_registers as f64
    }

    /// Merge another sketch into this one (union operation)
    pub fn merge(&mut self, other: &Self) {
        assert_eq!(
            self.num_registers, other.num_registers,
            "Sketches must have same number of registers"
        );

        for i in 0..self.num_registers {
            if other.registers[i] < self.registers[i] {
                self.registers[i] = other.registers[i];
            }
        }

        self.count += other.count;
    }

    /// Get the number of items added to the sketch
    pub fn count(&self) -> u64 {
        self.count
    }

    /// Get the number of registers
    pub fn num_registers(&self) -> usize {
        self.num_registers
    }

    /// Check if the sketch is empty
    pub fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Reset the sketch to initial state
    pub fn reset(&mut self) {
        for reg in &mut self.registers {
            *reg = f64::INFINITY;
        }
        self.count = 0;
    }

    /// Get a confidence score for a classification based on ANI
    ///
    /// # Arguments
    /// * `reference_sketch` - Sketch of the reference sequence
    /// * `k` - K-mer size
    /// * `min_ani` - Minimum ANI threshold (typically 0.95 for species-level)
    ///
    /// # Returns
    /// Confidence score in range [0, 1], where 1 is high confidence
    pub fn classification_confidence(
        &self,
        reference_sketch: &Self,
        k: usize,
        min_ani: f64,
    ) -> f64 {
        match self.estimate_ani(reference_sketch, k) {
            Some(ani) => {
                if ani >= min_ani {
                    // Scale ANI to confidence score
                    // ANI of 1.0 -> confidence 1.0
                    // ANI of min_ani -> confidence 0.5
                    // ANI below min_ani -> confidence scales down
                    let scale = (ani - min_ani) / (1.0 - min_ani);
                    0.5 + 0.5 * scale.max(0.0)
                } else {
                    // Below threshold: low confidence
                    ani / min_ani * 0.5
                }
            }
            None => 0.0,
        }
    }
}

/// Hash a value with a seed for register-specific hashing
#[inline]
fn hash_with_seed(value: u64, seed: u64) -> u64 {
    let mut x = value ^ seed;
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51afd7ed558ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
    x ^= x >> 33;
    x
}

/// Convert a 64-bit hash to a uniform [0, 1) distribution
#[inline]
fn hash_to_uniform(hash: u64) -> f64 {
    // Use the upper 53 bits for double precision
    (hash >> 11) as f64 / (1u64 << 53) as f64
}

/// ANI-based classification confidence calculator
///
/// This struct holds configuration for computing confidence scores
/// based on ANI estimates.
#[derive(Clone, Debug)]
pub struct AniConfidenceCalculator {
    /// K-mer size
    k: usize,
    /// Minimum ANI for species-level confidence
    species_threshold: f64,
    /// Minimum ANI for genus-level confidence
    genus_threshold: f64,
}

impl Default for AniConfidenceCalculator {
    fn default() -> Self {
        Self {
            k: 35,
            species_threshold: 0.95,
            genus_threshold: 0.80,
        }
    }
}

impl AniConfidenceCalculator {
    /// Create a new confidence calculator
    pub fn new(k: usize, species_threshold: f64, genus_threshold: f64) -> Self {
        Self {
            k,
            species_threshold,
            genus_threshold,
        }
    }

    /// Calculate confidence level based on ANI
    ///
    /// # Returns
    /// * `ConfidenceLevel` enum indicating the confidence category
    pub fn confidence_level(&self, query: &ProbMinHashSketch, reference: &ProbMinHashSketch) -> ConfidenceLevel {
        match query.estimate_ani(reference, self.k) {
            Some(ani) if ani >= self.species_threshold => ConfidenceLevel::HighConfidence(ani),
            Some(ani) if ani >= self.genus_threshold => ConfidenceLevel::MediumConfidence(ani),
            Some(ani) if ani > 0.0 => ConfidenceLevel::LowConfidence(ani),
            _ => ConfidenceLevel::NoMatch,
        }
    }
}

/// Confidence level for a classification
#[derive(Clone, Debug, PartialEq)]
pub enum ConfidenceLevel {
    /// High confidence match (ANI >= species threshold)
    HighConfidence(f64),
    /// Medium confidence match (genus threshold <= ANI < species threshold)
    MediumConfidence(f64),
    /// Low confidence match (0 < ANI < genus threshold)
    LowConfidence(f64),
    /// No significant match
    NoMatch,
}

impl ConfidenceLevel {
    /// Get the ANI value if present
    pub fn ani(&self) -> Option<f64> {
        match self {
            ConfidenceLevel::HighConfidence(ani)
            | ConfidenceLevel::MediumConfidence(ani)
            | ConfidenceLevel::LowConfidence(ani) => Some(*ani),
            ConfidenceLevel::NoMatch => None,
        }
    }

    /// Check if this is a confident match
    pub fn is_confident(&self) -> bool {
        matches!(self, ConfidenceLevel::HighConfidence(_) | ConfidenceLevel::MediumConfidence(_))
    }
}

impl PartialOrd for ConfidenceLevel {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self.ani(), other.ani()) {
            (Some(a), Some(b)) => a.partial_cmp(&b),
            (Some(_), None) => Some(Ordering::Greater),
            (None, Some(_)) => Some(Ordering::Less),
            (None, None) => Some(Ordering::Equal),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sketch_creation() {
        let sketch = ProbMinHashSketch::new(64);
        assert_eq!(sketch.num_registers(), 64);
        assert_eq!(sketch.count(), 0);
        assert!(sketch.is_empty());
    }

    #[test]
    fn test_sketch_default() {
        let sketch = ProbMinHashSketch::default();
        assert_eq!(sketch.num_registers(), DEFAULT_NUM_REGISTERS);
    }

    #[test]
    fn test_insert_and_count() {
        let mut sketch = ProbMinHashSketch::new(64);
        sketch.insert(12345);
        assert_eq!(sketch.count(), 1);
        assert!(!sketch.is_empty());

        sketch.insert(67890);
        assert_eq!(sketch.count(), 2);
    }

    #[test]
    fn test_insert_slice() {
        let mut sketch = ProbMinHashSketch::new(64);
        sketch.insert_slice(&[1, 2, 3, 4, 5]);
        assert_eq!(sketch.count(), 5);
    }

    #[test]
    fn test_jaccard_identical() {
        let mut sketch1 = ProbMinHashSketch::new(64);
        let mut sketch2 = ProbMinHashSketch::new(64);

        // Insert same elements
        for i in 0..100 {
            sketch1.insert(i);
            sketch2.insert(i);
        }

        let jaccard = sketch1.estimate_jaccard(&sketch2);
        assert!(jaccard > 0.9, "Identical sets should have Jaccard ~1.0, got {}", jaccard);
    }

    #[test]
    fn test_jaccard_disjoint() {
        let mut sketch1 = ProbMinHashSketch::new(64);
        let mut sketch2 = ProbMinHashSketch::new(64);

        // Insert different elements
        for i in 0..100 {
            sketch1.insert(i);
            sketch2.insert(i + 1000);
        }

        let jaccard = sketch1.estimate_jaccard(&sketch2);
        assert!(jaccard < 0.2, "Disjoint sets should have low Jaccard, got {}", jaccard);
    }

    #[test]
    fn test_jaccard_partial_overlap() {
        let mut sketch1 = ProbMinHashSketch::new(128);
        let mut sketch2 = ProbMinHashSketch::new(128);

        // 50% overlap
        for i in 0..100 {
            sketch1.insert(i);
        }
        for i in 50..150 {
            sketch2.insert(i);
        }

        let jaccard = sketch1.estimate_jaccard(&sketch2);
        // Expected Jaccard: 50 / 150 ≈ 0.33
        assert!(jaccard > 0.1 && jaccard < 0.6, "Partial overlap Jaccard should be moderate, got {}", jaccard);
    }

    #[test]
    fn test_ani_estimation() {
        let mut sketch1 = ProbMinHashSketch::new(128);
        let mut sketch2 = ProbMinHashSketch::new(128);

        // High overlap for high ANI
        for i in 0..1000 {
            sketch1.insert(i);
            sketch2.insert(i);
        }

        let ani = sketch1.estimate_ani(&sketch2, 31);
        assert!(ani.is_some());
        let ani_val = ani.unwrap();
        assert!(ani_val > 0.9, "Similar sequences should have high ANI, got {}", ani_val);
    }

    #[test]
    fn test_ani_low_similarity() {
        let mut sketch1 = ProbMinHashSketch::new(128);
        let mut sketch2 = ProbMinHashSketch::new(128);

        // No overlap
        for i in 0..100 {
            sketch1.insert(i);
            sketch2.insert(i + 10000);
        }

        let ani = sketch1.estimate_ani(&sketch2, 31);
        // May be None or very low
        if let Some(ani_val) = ani {
            assert!(ani_val < 0.5, "Dissimilar sequences should have low ANI");
        }
    }

    #[test]
    fn test_containment() {
        let mut sketch1 = ProbMinHashSketch::new(64);
        let mut sketch2 = ProbMinHashSketch::new(64);

        // sketch1 is subset of sketch2
        for i in 0..50 {
            sketch1.insert(i);
        }
        for i in 0..100 {
            sketch2.insert(i);
        }

        let containment = sketch1.estimate_containment(&sketch2);
        assert!(containment > 0.7, "Subset should have high containment, got {}", containment);
    }

    #[test]
    fn test_merge() {
        let mut sketch1 = ProbMinHashSketch::new(64);
        let mut sketch2 = ProbMinHashSketch::new(64);

        for i in 0..50 {
            sketch1.insert(i);
        }
        for i in 50..100 {
            sketch2.insert(i);
        }

        let count_before = sketch1.count();
        sketch1.merge(&sketch2);
        assert_eq!(sketch1.count(), count_before + sketch2.count());
    }

    #[test]
    fn test_reset() {
        let mut sketch = ProbMinHashSketch::new(64);
        sketch.insert(12345);
        assert!(!sketch.is_empty());

        sketch.reset();
        assert!(sketch.is_empty());
        assert_eq!(sketch.count(), 0);
    }

    #[test]
    fn test_hash_with_seed_deterministic() {
        let h1 = hash_with_seed(12345, 0);
        let h2 = hash_with_seed(12345, 0);
        assert_eq!(h1, h2);
    }

    #[test]
    fn test_hash_with_seed_different_seeds() {
        let h1 = hash_with_seed(12345, 0);
        let h2 = hash_with_seed(12345, 1);
        assert_ne!(h1, h2);
    }

    #[test]
    fn test_hash_to_uniform_range() {
        for i in 0..1000 {
            let uniform = hash_to_uniform(i * 12345);
            assert!(uniform >= 0.0 && uniform < 1.0, "Uniform value {} out of range", uniform);
        }
    }

    #[test]
    fn test_classification_confidence_high() {
        let mut query = ProbMinHashSketch::new(128);
        let mut reference = ProbMinHashSketch::new(128);

        // Identical sequences
        for i in 0..1000 {
            query.insert(i);
            reference.insert(i);
        }

        let confidence = query.classification_confidence(&reference, 31, 0.95);
        assert!(confidence > 0.7, "High ANI should give high confidence, got {}", confidence);
    }

    #[test]
    fn test_confidence_level_ordering() {
        let high = ConfidenceLevel::HighConfidence(0.98);
        let medium = ConfidenceLevel::MediumConfidence(0.85);
        let low = ConfidenceLevel::LowConfidence(0.5);
        let none = ConfidenceLevel::NoMatch;

        assert!(high > medium);
        assert!(medium > low);
        assert!(low > none);
    }

    #[test]
    fn test_confidence_level_is_confident() {
        assert!(ConfidenceLevel::HighConfidence(0.98).is_confident());
        assert!(ConfidenceLevel::MediumConfidence(0.85).is_confident());
        assert!(!ConfidenceLevel::LowConfidence(0.5).is_confident());
        assert!(!ConfidenceLevel::NoMatch.is_confident());
    }

    #[test]
    fn test_ani_confidence_calculator() {
        let calc = AniConfidenceCalculator::default();
        assert_eq!(calc.k, 35);
        assert_eq!(calc.species_threshold, 0.95);
        assert_eq!(calc.genus_threshold, 0.80);
    }

    #[test]
    fn test_ani_confidence_calculator_custom() {
        let calc = AniConfidenceCalculator::new(31, 0.97, 0.85);
        assert_eq!(calc.k, 31);
        assert_eq!(calc.species_threshold, 0.97);
        assert_eq!(calc.genus_threshold, 0.85);
    }

    #[test]
    fn test_confidence_level_ani() {
        assert_eq!(ConfidenceLevel::HighConfidence(0.98).ani(), Some(0.98));
        assert_eq!(ConfidenceLevel::NoMatch.ani(), None);
    }

    #[test]
    fn test_sketch_clone() {
        let mut sketch = ProbMinHashSketch::new(64);
        sketch.insert(12345);

        let cloned = sketch.clone();
        assert_eq!(cloned.count(), sketch.count());
        assert_eq!(cloned.num_registers(), sketch.num_registers());
    }
}

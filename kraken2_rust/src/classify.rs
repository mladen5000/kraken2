/// Taxonomic sequence classification engine
///
/// Translated from classify.cc
///
/// Main classification algorithm that assigns sequences to taxa based on minimizer hits.
/// Implements full RTL (Right-to-Left) confidence scoring, paired-end processing,
/// translated search (6-frame translation), and quality score filtering.

use crate::kraken2_data::{IndexOptions, TaxId, TAXID_MAX};
use crate::taxonomy::Taxonomy;
use crate::mmscanner::MinimizerScanner;
use crate::seqreader::{Sequence, SequenceReader};
use crate::reports::{ReportGenerator, ReportFormat};
use crate::compact_hash::CompactHashTable;
use crate::aa_translate::translate_to_all_frames;
use anyhow::Result;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;

/// Number of sequence fragments to process per thread batch
#[allow(dead_code)]
const NUM_FRAGMENTS_PER_THREAD: usize = 10000;

/// Special taxon marker for paired-end read boundary
const MATE_PAIR_BORDER_TAXON: TaxId = TAXID_MAX;

/// Special taxon marker for reading frame boundary in translated search
const READING_FRAME_BORDER_TAXON: TaxId = TAXID_MAX - 1;

/// Special taxon marker for ambiguous (N-containing) spans
const AMBIGUOUS_SPAN_TAXON: TaxId = TAXID_MAX - 2;

/// Hit counts type - maps taxon ID to count of hits
pub type TaxonCounts = HashMap<TaxId, u32>;

/// Classification result with detailed hit information
#[derive(Clone, Debug)]
pub struct DetailedClassification {
    pub sequence_id: String,
    pub classified: bool,
    pub taxid: TaxId,
    pub external_taxid: TaxId,
    pub kmers_matched: usize,
    pub kmers_total: usize,
    pub hit_groups: usize,
    pub hitlist: Vec<TaxId>,
    pub hitlist_string: String,
}

impl Default for DetailedClassification {
    fn default() -> Self {
        DetailedClassification {
            sequence_id: String::new(),
            classified: false,
            taxid: 0,
            external_taxid: 0,
            kmers_matched: 0,
            kmers_total: 0,
            hit_groups: 0,
            hitlist: Vec::new(),
            hitlist_string: String::new(),
        }
    }
}

/// Options for sequence classification
pub struct ClassifyOptions {
    pub index_filename: String,
    pub taxonomy_filename: String,
    pub options_filename: String,
    pub report_filename: String,
    pub classified_output_filename: Option<String>,
    pub unclassified_output_filename: Option<String>,
    pub kraken_output_filename: Option<String>,
    pub mpa_style_report: bool,
    pub report_kmer_data: bool,
    pub quick_mode: bool,
    pub report_zero_counts: bool,
    pub use_translated_search: bool,
    pub print_scientific_name: bool,
    pub confidence_threshold: f64,
    pub num_threads: usize,
    pub paired_end_processing: bool,
    pub single_file_pairs: bool,
    pub minimum_quality_score: i32,
    pub minimum_hit_groups: i32,
    pub use_memory_mapping: bool,
}

impl Default for ClassifyOptions {
    fn default() -> Self {
        ClassifyOptions {
            index_filename: String::new(),
            taxonomy_filename: String::new(),
            options_filename: String::new(),
            report_filename: String::new(),
            classified_output_filename: None,
            unclassified_output_filename: None,
            kraken_output_filename: None,
            mpa_style_report: false,
            report_kmer_data: false,
            quick_mode: false,
            report_zero_counts: false,
            use_translated_search: false,
            print_scientific_name: false,
            confidence_threshold: 0.0,
            num_threads: 1,
            paired_end_processing: false,
            single_file_pairs: false,
            minimum_quality_score: 0,
            minimum_hit_groups: 0,
            use_memory_mapping: false,
        }
    }
}

/// Classification statistics
#[derive(Default, Debug, Clone)]
pub struct ClassificationStats {
    pub total_sequences: u64,
    pub total_bases: u64,
    pub total_classified: u64,
}

impl ClassificationStats {
    pub fn classification_rate(&self) -> f64 {
        if self.total_sequences == 0 {
            0.0
        } else {
            self.total_classified as f64 / self.total_sequences as f64
        }
    }
}

/// Resolve the taxonomy tree to find the best classification
///
/// Uses RTL (Right-to-Left) confidence scoring:
/// 1. For each taxon, compute sum of hits for all taxa in its subtree
/// 2. Find taxon with highest score
/// 3. If score doesn't meet confidence threshold, walk up tree until it does
///
/// This is a direct translation of the C++ ResolveTree function.
pub fn resolve_tree(
    hit_counts: &TaxonCounts,
    taxonomy: &Taxonomy,
    total_minimizers: usize,
    confidence_threshold: f64,
) -> TaxId {
    if hit_counts.is_empty() || total_minimizers == 0 {
        return 0;
    }

    let mut max_taxon: TaxId = 0;
    let mut max_score: u32 = 0;
    let required_score = (confidence_threshold * total_minimizers as f64).ceil() as u32;

    // Sum each taxon's LTR path, find taxon with highest LTR score
    for (&taxon, _) in hit_counts {
        let mut score: u32 = 0;

        // Sum hits from all taxa that are descendants of this taxon
        for (&taxon2, &count) in hit_counts {
            if taxonomy.is_a_ancestor_of_b(taxon, taxon2) {
                score += count;
            }
        }

        if score > max_score {
            max_score = score;
            max_taxon = taxon;
        } else if score == max_score && max_taxon != 0 {
            // Tie-breaker: use LCA
            max_taxon = taxonomy.lca(max_taxon, taxon);
        }
    }

    // Reset max score to be only hits at the called taxon
    max_score = *hit_counts.get(&max_taxon).unwrap_or(&0);

    // Walk up tree until confidence threshold is met
    while max_taxon != 0 && max_score < required_score {
        max_score = 0;
        for (&taxon, &count) in hit_counts {
            // Add to score if taxon is in max_taxon's clade
            if taxonomy.is_a_ancestor_of_b(max_taxon, taxon) {
                max_score += count;
            }
        }

        if max_score >= required_score {
            return max_taxon;
        } else {
            // Move up to parent
            max_taxon = taxonomy.get_parent(max_taxon);
        }
    }

    max_taxon
}

/// Build a hitlist string for Kraken-style output
///
/// Format: "taxid:count taxid:count ..." with special markers:
/// - A:count for ambiguous spans
/// - |:| for mate pair boundary
/// - -:- for reading frame boundary
pub fn add_hitlist_string(taxa: &[TaxId], taxonomy: &Taxonomy) -> String {
    if taxa.is_empty() {
        return "0:0".to_string();
    }

    let mut result = String::new();
    let mut last_code = taxa[0];
    let mut code_count = 1;

    for &code in taxa.iter().skip(1) {
        if code == last_code {
            code_count += 1;
        } else {
            // Output the previous run
            if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
                if last_code == AMBIGUOUS_SPAN_TAXON {
                    result.push_str(&format!("A:{} ", code_count));
                } else {
                    let ext_code = taxonomy.get_external_id(last_code);
                    result.push_str(&format!("{}:{} ", ext_code, code_count));
                }
            } else {
                // Mate pair or reading frame marker
                if last_code == MATE_PAIR_BORDER_TAXON {
                    result.push_str("|:| ");
                } else {
                    result.push_str("-:- ");
                }
            }
            code_count = 1;
            last_code = code;
        }
    }

    // Output final run
    if last_code != MATE_PAIR_BORDER_TAXON && last_code != READING_FRAME_BORDER_TAXON {
        if last_code == AMBIGUOUS_SPAN_TAXON {
            result.push_str(&format!("A:{}", code_count));
        } else {
            let ext_code = taxonomy.get_external_id(last_code);
            result.push_str(&format!("{}:{}", ext_code, code_count));
        }
    } else {
        if last_code == MATE_PAIR_BORDER_TAXON {
            result.push_str("|:|");
        } else {
            result.push_str("-:-");
        }
    }

    result
}

/// Mask low-quality bases in a FASTQ sequence
///
/// Replaces bases with quality score below threshold with 'x'
pub fn mask_low_quality_bases(seq: &mut Sequence, minimum_quality_score: i32) {
    if let Some(ref quals) = seq.quality {
        if seq.seq.len() != quals.len() {
            return; // Invalid quality string length
        }

        let mut seq_bytes = seq.seq.as_bytes().to_vec();
        for (i, &q) in quals.as_bytes().iter().enumerate() {
            let quality = (q as i32) - ('!' as i32);
            if quality < minimum_quality_score {
                seq_bytes[i] = b'x';
            }
        }
        seq.seq = String::from_utf8_lossy(&seq_bytes).to_string();
    }
}

/// Trim pair info from sequence header (e.g., "/1" or "/2")
pub fn trim_pair_info(id: &str) -> &str {
    let sz = id.len();
    if sz <= 2 {
        return id;
    }
    let bytes = id.as_bytes();
    if bytes[sz - 2] == b'/' && (bytes[sz - 1] == b'1' || bytes[sz - 1] == b'2') {
        return &id[..sz - 2];
    }
    id
}

/// Classify a single sequence (or pair of sequences)
///
/// Full implementation matching the C++ ClassifySequence function.
/// Supports:
/// - Single-end and paired-end reads
/// - Translated search (6-frame translation for protein databases)
/// - Quick mode (early termination)
/// - Minimum hit groups filtering
/// - Confidence threshold scoring
pub fn classify_sequence_full(
    dna1: &Sequence,
    dna2: Option<&Sequence>,
    hash: &CompactHashTable,
    taxonomy: &Taxonomy,
    idx_opts: &IndexOptions,
    opts: &ClassifyOptions,
    scanner: &mut MinimizerScanner,
) -> DetailedClassification {
    let mut taxa: Vec<TaxId> = Vec::new();
    let mut hit_counts: TaxonCounts = HashMap::new();
    let _frame_ct = if opts.use_translated_search { 6 } else { 1 };
    let mut minimizer_hit_groups: i32 = 0;
    let mut call: TaxId = 0;

    // Process each mate (1 for single-end, 2 for paired-end)
    let mate_count = if dna2.is_some() && opts.paired_end_processing { 2 } else { 1 };

    'mate_loop: for mate_num in 0..mate_count {
        let current_seq = if mate_num == 0 { dna1 } else { dna2.unwrap() };

        // Get frames for translated search or just the original sequence
        let frames: Vec<String> = if opts.use_translated_search {
            translate_to_all_frames(&current_seq.seq)
        } else {
            vec![current_seq.seq.clone()]
        };

        for (frame_idx, frame) in frames.iter().enumerate() {
            scanner.load_sequence(frame);

            let mut last_minimizer: u64 = u64::MAX;
            let mut last_taxon: TaxId = TAXID_MAX;

            while let Some(minimizer) = scanner.next_minimizer() {
                let taxon: TaxId;

                if scanner.is_ambiguous() {
                    taxon = AMBIGUOUS_SPAN_TAXON;
                } else {
                    if minimizer != last_minimizer {
                        // Check minimum acceptable hash value if set
                        let skip_lookup = if idx_opts.minimum_acceptable_hash_value > 0 {
                            murmur_hash3(minimizer) < idx_opts.minimum_acceptable_hash_value
                        } else {
                            false
                        };

                        let looked_up_taxon = if !skip_lookup {
                            hash.get(minimizer)
                        } else {
                            0
                        };

                        last_taxon = looked_up_taxon;
                        last_minimizer = minimizer;

                        if looked_up_taxon != 0 {
                            minimizer_hit_groups += 1;
                        }

                        taxon = looked_up_taxon;
                    } else {
                        taxon = last_taxon;
                    }

                    if taxon != 0 {
                        // Quick mode: return immediately if we have enough hit groups
                        if opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups {
                            call = taxon;
                            break 'mate_loop;
                        }
                        *hit_counts.entry(taxon).or_insert(0) += 1;
                    }
                }

                taxa.push(taxon);
            }

            // Add reading frame border between frames (except after last frame)
            if opts.use_translated_search && frame_idx != 5 {
                taxa.push(READING_FRAME_BORDER_TAXON);
            }
        }

        // Add mate pair border between mates
        if opts.paired_end_processing && mate_num == 0 && dna2.is_some() {
            taxa.push(MATE_PAIR_BORDER_TAXON);
        }
    }

    // Calculate total k-mers (excluding markers)
    let mut total_kmers = taxa.len();
    if opts.paired_end_processing && dna2.is_some() {
        total_kmers = total_kmers.saturating_sub(1); // Account for mate pair marker
    }
    if opts.use_translated_search {
        // Account for reading frame markers (5 per mate for 6 frames)
        let frame_markers = if opts.paired_end_processing && dna2.is_some() { 10 } else { 5 };
        total_kmers = total_kmers.saturating_sub(frame_markers);
    }

    // Resolve tree to get final classification
    if call == 0 {
        call = resolve_tree(&hit_counts, taxonomy, total_kmers, opts.confidence_threshold);
    }

    // Void a call made by too few minimizer hit groups
    if call != 0 && minimizer_hit_groups < opts.minimum_hit_groups {
        call = 0;
    }

    // Build result
    let external_taxid = if call != 0 {
        taxonomy.get_external_id(call)
    } else {
        0
    };

    let kmers_matched: usize = hit_counts.values().map(|&v| v as usize).sum();

    let hitlist_string = if opts.quick_mode && call != 0 {
        format!("{}:Q", external_taxid)
    } else if taxa.is_empty() {
        "0:0".to_string()
    } else {
        add_hitlist_string(&taxa, taxonomy)
    };

    DetailedClassification {
        sequence_id: dna1.header.clone(),
        classified: call != 0,
        taxid: call,
        external_taxid,
        kmers_matched,
        kmers_total: total_kmers,
        hit_groups: minimizer_hit_groups as usize,
        hitlist: taxa,
        hitlist_string,
    }
}

/// Simple MurmurHash3 finalizer for hash value checking
fn murmur_hash3(key: u64) -> u64 {
    let mut k = key;
    k ^= k >> 33;
    k = k.wrapping_mul(0xff51afd7ed558ccd);
    k ^= k >> 33;
    k = k.wrapping_mul(0xc4ceb9fe1a85ec53);
    k ^= k >> 33;
    k
}

/// Format Kraken-style output line for a classification result
pub fn format_kraken_output(
    result: &DetailedClassification,
    dna1: &Sequence,
    dna2: Option<&Sequence>,
    taxonomy: &Taxonomy,
    opts: &ClassifyOptions,
) -> String {
    let mut output = String::new();

    // Classification status
    if result.classified {
        output.push_str("C\t");
    } else {
        output.push_str("U\t");
    }

    // Sequence header (trimmed for paired-end)
    if opts.paired_end_processing {
        output.push_str(trim_pair_info(&dna1.header));
    } else {
        output.push_str(&dna1.header);
    }
    output.push('\t');

    // Taxon ID or name
    if opts.print_scientific_name {
        if result.classified {
            if let Some(name) = taxonomy.get_name(result.taxid) {
                output.push_str(&format!("{} (taxid {})", name, result.external_taxid));
            } else {
                output.push_str(&format!("taxid {}", result.external_taxid));
            }
        } else {
            output.push_str("unclassified (taxid 0)");
        }
    } else {
        output.push_str(&result.external_taxid.to_string());
    }
    output.push('\t');

    // Sequence length(s)
    if opts.paired_end_processing && dna2.is_some() {
        output.push_str(&format!("{}|{}", dna1.seq.len(), dna2.unwrap().seq.len()));
    } else {
        output.push_str(&dna1.seq.len().to_string());
    }
    output.push('\t');

    // Hitlist
    output.push_str(&result.hitlist_string);
    output.push('\n');

    output
}

/// Classify sequences using a pre-built database (simplified version)
pub fn classify_sequences(
    opts: &ClassifyOptions,
    mut sequence_reader: SequenceReader,
) -> Result<ClassificationStats> {
    // Load taxonomy
    let taxonomy = Taxonomy::load_from_ncbi(&opts.taxonomy_filename, "")?;

    // Load index options
    let idx_opts = load_index_options(&opts.options_filename)?;

    // Initialize report generator
    let report_format = if opts.mpa_style_report {
        ReportFormat::Mpa
    } else {
        ReportFormat::Kraken
    };
    let _report_gen = ReportGenerator::new(report_format);

    // Initialize output files
    let mut _report_writer: Option<BufWriter<File>> = None;
    if !opts.report_filename.is_empty() {
        let file = File::create(&opts.report_filename)?;
        _report_writer = Some(BufWriter::new(file));
    }

    // Statistics accumulator
    let mut stats = ClassificationStats::default();

    // Process sequences using simplified classification
    while let Some(seq) = sequence_reader.next_sequence()? {
        stats.total_sequences += 1;
        stats.total_bases += seq.seq.len() as u64;

        // Perform simplified classification
        let result = classify_sequence_simple(&seq.seq, &taxonomy, &idx_opts);

        if result.classified {
            stats.total_classified += 1;
        }
    }

    eprintln!("Classification complete");
    eprintln!("Processed {} sequences", stats.total_sequences);
    eprintln!("Classified {:.2}%", stats.classification_rate() * 100.0);

    Ok(stats)
}

/// Classify a single sequence (simplified version without hash table)
fn classify_sequence_simple(
    seq: &str,
    _taxonomy: &Taxonomy,
    idx_opts: &IndexOptions,
) -> DetailedClassification {
    let scanner = MinimizerScanner::new(idx_opts);
    let minimizers = scanner.scan(seq);
    let num_minimizers = minimizers.len();

    if num_minimizers == 0 {
        return DetailedClassification::default();
    }

    // In a full implementation, we would look up each minimizer in the hash table
    // For now, return unclassified since we don't have a hash table loaded
    DetailedClassification {
        sequence_id: String::new(),
        classified: false,
        taxid: 0,
        external_taxid: 0,
        kmers_matched: 0,
        kmers_total: num_minimizers,
        hit_groups: 0,
        hitlist: Vec::new(),
        hitlist_string: "0:0".to_string(),
    }
}

/// Load index options from binary file
fn load_index_options(filename: &str) -> Result<IndexOptions> {
    use std::fs;
    use std::io::Read;

    // Try to read the binary options file
    if let Ok(mut file) = fs::File::open(filename) {
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        // Parse binary format (matches C++ struct layout)
        if buffer.len() >= std::mem::size_of::<IndexOptions>() {
            // For safety, we create a default and update fields
            let mut opts = IndexOptions::default();

            // Read fields in order (assuming little-endian, matching C++ layout)
            if buffer.len() >= 8 {
                opts.k = u32::from_le_bytes([buffer[0], buffer[1], buffer[2], buffer[3]]) as usize;
                opts.l = u32::from_le_bytes([buffer[4], buffer[5], buffer[6], buffer[7]]) as usize;
            }
            if buffer.len() >= 16 {
                opts.spaced_seed_mask = u64::from_le_bytes([
                    buffer[8], buffer[9], buffer[10], buffer[11],
                    buffer[12], buffer[13], buffer[14], buffer[15],
                ]);
            }
            if buffer.len() >= 24 {
                opts.toggle_mask = u64::from_le_bytes([
                    buffer[16], buffer[17], buffer[18], buffer[19],
                    buffer[20], buffer[21], buffer[22], buffer[23],
                ]);
            }
            if buffer.len() >= 25 {
                opts.dna_db = buffer[24] != 0;
            }
            if buffer.len() >= 32 {
                opts.minimum_acceptable_hash_value = u64::from_le_bytes([
                    buffer[25], buffer[26], buffer[27], buffer[28],
                    buffer[29], buffer[30], buffer[31], buffer[32].min(0),
                ]);
            }

            eprintln!("Loaded index options from {} (k={}, l={})", filename, opts.k, opts.l);
            return Ok(opts);
        }
    }

    // Fall back to defaults
    let opts = IndexOptions::default();
    eprintln!("Using default index options (file not found or invalid: {})", filename);
    Ok(opts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classification_stats_rate() {
        let mut stats = ClassificationStats::default();
        stats.total_sequences = 100;
        stats.total_classified = 75;

        assert_eq!(stats.classification_rate(), 0.75);
    }

    #[test]
    fn test_classification_stats_rate_zero() {
        let stats = ClassificationStats::default();
        assert_eq!(stats.classification_rate(), 0.0);
    }

    #[test]
    fn test_classify_options_default() {
        let opts = ClassifyOptions::default();
        assert_eq!(opts.num_threads, 1);
        assert!(!opts.quick_mode);
        assert!(!opts.paired_end_processing);
        assert!(!opts.single_file_pairs);
    }

    #[test]
    fn test_detailed_classification_default() {
        let result = DetailedClassification::default();
        assert!(!result.classified);
        assert_eq!(result.taxid, 0);
        assert_eq!(result.external_taxid, 0);
        assert_eq!(result.hit_groups, 0);
    }

    #[test]
    fn test_detailed_classification_creation() {
        let result = DetailedClassification {
            sequence_id: "seq1".to_string(),
            classified: true,
            taxid: 562,
            external_taxid: 562,
            kmers_matched: 30,
            kmers_total: 35,
            hit_groups: 5,
            hitlist: vec![562, 561, 562],
            hitlist_string: "562:2 561:1".to_string(),
        };

        assert_eq!(result.sequence_id, "seq1");
        assert!(result.classified);
        assert_eq!(result.taxid, 562);
        assert_eq!(result.hit_groups, 5);
    }

    #[test]
    fn test_resolve_tree_empty() {
        let hit_counts = TaxonCounts::new();
        let taxonomy = Taxonomy::new();
        let result = resolve_tree(&hit_counts, &taxonomy, 100, 0.5);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_resolve_tree_single_taxon() {
        let mut hit_counts = TaxonCounts::new();
        hit_counts.insert(562, 10);

        let taxonomy = Taxonomy::new();
        let result = resolve_tree(&hit_counts, &taxonomy, 10, 0.0);
        // With no confidence threshold, should return the taxon
        assert_eq!(result, 562);
    }

    #[test]
    fn test_trim_pair_info() {
        assert_eq!(trim_pair_info("seq1/1"), "seq1");
        assert_eq!(trim_pair_info("seq1/2"), "seq1");
        assert_eq!(trim_pair_info("seq1"), "seq1");
        assert_eq!(trim_pair_info("a"), "a");
        assert_eq!(trim_pair_info(""), "");
    }

    #[test]
    fn test_mask_low_quality_bases() {
        let mut seq = Sequence {
            header: "test".to_string(),
            seq: "ATCG".to_string(),
            quality: Some("!5II".to_string()), // Quality scores: 0, 20, 40, 40
        };

        mask_low_quality_bases(&mut seq, 10);
        assert_eq!(seq.seq, "xTCG"); // First base masked (quality 0 < 10)
    }

    #[test]
    fn test_add_hitlist_string_simple() {
        let taxa = vec![562, 562, 562, 561, 561];
        let taxonomy = Taxonomy::new();
        let result = add_hitlist_string(&taxa, &taxonomy);
        // Should show runs of same taxon
        assert!(result.contains("562:3"));
        assert!(result.contains("561:2"));
    }

    #[test]
    fn test_add_hitlist_string_with_ambiguous() {
        let taxa = vec![562, AMBIGUOUS_SPAN_TAXON, AMBIGUOUS_SPAN_TAXON, 561];
        let taxonomy = Taxonomy::new();
        let result = add_hitlist_string(&taxa, &taxonomy);
        assert!(result.contains("A:2")); // Ambiguous span
    }

    #[test]
    fn test_murmur_hash3() {
        // Test that hash produces deterministic results
        let h1 = murmur_hash3(12345);
        let h2 = murmur_hash3(12345);
        assert_eq!(h1, h2);

        // Different inputs should produce different outputs
        let h3 = murmur_hash3(12346);
        assert_ne!(h1, h3);
    }

    #[test]
    fn test_mate_pair_border_constants() {
        // Verify constants are set correctly
        assert_eq!(MATE_PAIR_BORDER_TAXON, TAXID_MAX);
        assert_eq!(READING_FRAME_BORDER_TAXON, TAXID_MAX - 1);
        assert_eq!(AMBIGUOUS_SPAN_TAXON, TAXID_MAX - 2);
    }

    #[test]
    fn test_classification_stats_clone() {
        let mut stats = ClassificationStats::default();
        stats.total_sequences = 100;
        stats.total_classified = 75;
        stats.total_bases = 5000;

        let cloned = stats.clone();
        assert_eq!(cloned.total_sequences, 100);
        assert_eq!(cloned.total_classified, 75);
        assert_eq!(cloned.total_bases, 5000);
    }

    #[test]
    fn test_classification_stats_debug() {
        let stats = ClassificationStats::default();
        let debug = format!("{:?}", stats);
        assert!(debug.contains("ClassificationStats"));
    }

    #[test]
    fn test_detailed_classification_clone() {
        let result = DetailedClassification {
            sequence_id: "seq1".to_string(),
            classified: true,
            taxid: 562,
            external_taxid: 562,
            kmers_matched: 30,
            kmers_total: 35,
            hit_groups: 5,
            hitlist: vec![562, 561],
            hitlist_string: "562:1 561:1".to_string(),
        };

        let cloned = result.clone();
        assert_eq!(cloned.sequence_id, "seq1");
        assert_eq!(cloned.taxid, 562);
        assert_eq!(cloned.hitlist.len(), 2);
    }

    #[test]
    fn test_detailed_classification_debug() {
        let result = DetailedClassification::default();
        let debug = format!("{:?}", result);
        assert!(debug.contains("DetailedClassification"));
    }

    #[test]
    fn test_mask_low_quality_bases_no_quality() {
        let mut seq = Sequence {
            header: "test".to_string(),
            seq: "ATCG".to_string(),
            quality: None,
        };
        mask_low_quality_bases(&mut seq, 10);
        // Should remain unchanged when no quality string
        assert_eq!(seq.seq, "ATCG");
    }

    #[test]
    fn test_mask_low_quality_bases_length_mismatch() {
        let mut seq = Sequence {
            header: "test".to_string(),
            seq: "ATCGATCG".to_string(),
            quality: Some("IIII".to_string()), // Wrong length
        };
        mask_low_quality_bases(&mut seq, 10);
        // Should remain unchanged with mismatched lengths
        assert_eq!(seq.seq, "ATCGATCG");
    }

    #[test]
    fn test_mask_low_quality_bases_all_high_quality() {
        let mut seq = Sequence {
            header: "test".to_string(),
            seq: "ATCG".to_string(),
            quality: Some("IIII".to_string()), // Quality 40 for all
        };
        mask_low_quality_bases(&mut seq, 10);
        // All high quality, should remain unchanged
        assert_eq!(seq.seq, "ATCG");
    }

    #[test]
    fn test_mask_low_quality_bases_all_low_quality() {
        let mut seq = Sequence {
            header: "test".to_string(),
            seq: "ATCG".to_string(),
            quality: Some("!!!!".to_string()), // Quality 0 for all
        };
        mask_low_quality_bases(&mut seq, 10);
        // All low quality, should all be masked
        assert_eq!(seq.seq, "xxxx");
    }

    #[test]
    fn test_trim_pair_info_variations() {
        assert_eq!(trim_pair_info("read1/1"), "read1");
        assert_eq!(trim_pair_info("read1/2"), "read1");
        assert_eq!(trim_pair_info("read1/3"), "read1/3"); // Not /1 or /2
        assert_eq!(trim_pair_info("ab"), "ab"); // Too short to have /1
        assert_eq!(trim_pair_info("/1"), "/1"); // Edge case
    }

    #[test]
    fn test_add_hitlist_string_empty() {
        let taxa: Vec<TaxId> = vec![];
        let taxonomy = Taxonomy::new();
        let result = add_hitlist_string(&taxa, &taxonomy);
        assert_eq!(result, "0:0");
    }

    #[test]
    fn test_add_hitlist_string_single_taxon() {
        let taxa = vec![562];
        let taxonomy = Taxonomy::new();
        let result = add_hitlist_string(&taxa, &taxonomy);
        assert_eq!(result, "562:1");
    }

    #[test]
    fn test_add_hitlist_string_mate_pair_border() {
        let taxa = vec![562, MATE_PAIR_BORDER_TAXON, 561];
        let taxonomy = Taxonomy::new();
        let result = add_hitlist_string(&taxa, &taxonomy);
        assert!(result.contains("|:|"));
    }

    #[test]
    fn test_add_hitlist_string_reading_frame_border() {
        let taxa = vec![562, READING_FRAME_BORDER_TAXON, 561];
        let taxonomy = Taxonomy::new();
        let result = add_hitlist_string(&taxa, &taxonomy);
        assert!(result.contains("-:-"));
    }

    #[test]
    fn test_murmur_hash3_known_values() {
        // Test with 0
        let h0 = murmur_hash3(0);
        assert_ne!(h0, 0); // Should mix even 0

        // Test with max value
        let hmax = murmur_hash3(u64::MAX);
        assert_ne!(hmax, u64::MAX);
    }

    #[test]
    fn test_resolve_tree_zero_minimizers() {
        let mut hit_counts = TaxonCounts::new();
        hit_counts.insert(562, 10);

        let taxonomy = Taxonomy::new();
        let result = resolve_tree(&hit_counts, &taxonomy, 0, 0.5);
        assert_eq!(result, 0);
    }

    #[test]
    fn test_classify_options_custom() {
        let opts = ClassifyOptions {
            num_threads: 8,
            quick_mode: true,
            paired_end_processing: true,
            confidence_threshold: 0.5,
            minimum_hit_groups: 3,
            ..Default::default()
        };

        assert_eq!(opts.num_threads, 8);
        assert!(opts.quick_mode);
        assert!(opts.paired_end_processing);
        assert_eq!(opts.confidence_threshold, 0.5);
        assert_eq!(opts.minimum_hit_groups, 3);
    }
}

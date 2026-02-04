/// Taxonomy management for NCBI taxonomy tree
///
/// Translated from taxonomy.cc/taxonomy.h
///
/// Manages the NCBI taxonomy hierarchy, including parent-child relationships,
/// scientific names, and Lowest Common Ancestor (LCA) computation

use crate::kraken2_data::TaxId;
use std::collections::HashMap;

/// Represents a taxonomy node
#[derive(Clone, Debug)]
pub struct TaxonomyNode {
    pub taxid: TaxId,
    pub parent: TaxId,
    pub rank: String,
    pub name: String,
}

/// Manages the taxonomy hierarchy
pub struct Taxonomy {
    nodes: HashMap<TaxId, TaxonomyNode>,
    root_id: TaxId,
}

impl Taxonomy {
    /// Create a new empty taxonomy
    pub fn new() -> Self {
        Taxonomy {
            nodes: HashMap::new(),
            root_id: 1, // NCBI root is always TaxID 1
        }
    }

    /// Load taxonomy from NCBI nodes.dmp and names.dmp files
    pub fn load_from_ncbi(_nodes_file: &str, _names_file: &str) -> std::io::Result<Self> {
        // TODO: Implement loading from actual NCBI files
        Ok(Taxonomy::new())
    }

    /// Add a node to the taxonomy
    pub fn add_node(&mut self, taxid: TaxId, parent: TaxId, rank: &str, name: &str) {
        let node = TaxonomyNode {
            taxid,
            parent,
            rank: rank.to_string(),
            name: name.to_string(),
        };
        self.nodes.insert(taxid, node);
    }

    /// Get a taxonomy node by ID
    pub fn get_node(&self, taxid: TaxId) -> Option<&TaxonomyNode> {
        self.nodes.get(&taxid)
    }

    /// Get the parent of a taxon (returns 0 if not found)
    pub fn get_parent(&self, taxid: TaxId) -> TaxId {
        self.nodes.get(&taxid).map(|node| node.parent).unwrap_or(0)
    }

    /// Get the parent of a taxon as Option
    pub fn get_parent_opt(&self, taxid: TaxId) -> Option<TaxId> {
        self.nodes.get(&taxid).map(|node| node.parent)
    }

    /// Find the Lowest Common Ancestor (LCA) of two taxa
    pub fn lca(&self, taxid1: TaxId, taxid2: TaxId) -> TaxId {
        if taxid1 == taxid2 {
            return taxid1;
        }

        // Build path to root for first taxon
        let mut path1 = std::collections::HashSet::new();
        let mut current = taxid1;
        while let Some(node) = self.nodes.get(&current) {
            path1.insert(current);
            if current == node.parent {
                break; // Reached root
            }
            current = node.parent;
        }

        // Walk from second taxon until we find a common ancestor
        let mut current = taxid2;
        loop {
            if path1.contains(&current) {
                return current;
            }

            if let Some(node) = self.nodes.get(&current) {
                if current == node.parent {
                    break; // Reached root without finding LCA
                }
                current = node.parent;
            } else {
                break;
            }
        }

        self.root_id
    }

    /// Find the LCA of multiple taxa
    pub fn lca_many(&self, taxa: &[TaxId]) -> TaxId {
        if taxa.is_empty() {
            return self.root_id;
        }

        let mut result = taxa[0];
        for &taxid in &taxa[1..] {
            result = self.lca(result, taxid);
        }
        result
    }

    /// Get the scientific name of a taxon
    pub fn get_name(&self, taxid: TaxId) -> Option<&str> {
        self.nodes.get(&taxid).map(|node| node.name.as_str())
    }

    /// Get the rank of a taxon
    pub fn get_rank(&self, taxid: TaxId) -> Option<&str> {
        self.nodes.get(&taxid).map(|node| node.rank.as_str())
    }

    /// Get the root taxon ID (usually 1 for NCBI taxonomy)
    pub fn root(&self) -> TaxId {
        self.root_id
    }

    /// Get the number of taxa in the taxonomy
    pub fn size(&self) -> usize {
        self.nodes.len()
    }

    /// Check if taxon A is an ancestor of taxon B (including A == B)
    ///
    /// This implements IsAAncestorOfB from the C++ code.
    /// Returns true if A is on the path from B to the root.
    pub fn is_a_ancestor_of_b(&self, a: TaxId, b: TaxId) -> bool {
        if a == b {
            return true;
        }

        // Walk from B towards root, checking for A
        let mut current = b;
        while let Some(node) = self.nodes.get(&current) {
            if current == a {
                return true;
            }
            if current == node.parent {
                break; // Reached root
            }
            current = node.parent;
        }

        false
    }

    /// Get the external ID for a taxon (for now, same as internal ID)
    ///
    /// In the full C++ implementation, nodes have separate internal and external IDs.
    /// External IDs are what's shown in output (NCBI taxonomy IDs).
    pub fn get_external_id(&self, taxid: TaxId) -> TaxId {
        // For now, return the taxid itself as external ID
        // In full implementation, would lookup node's external_id field
        taxid
    }
}

impl Default for Taxonomy {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_taxonomy() -> Taxonomy {
        let mut tax = Taxonomy::new();
        // Build a simple taxonomy tree:
        // 1 (root) -> 2 (Bacteria) -> 3 (Firmicutes) -> 4 (Bacilli)
        //          -> 5 (Archaea) -> 6 (Euryarchaeota)
        tax.add_node(1, 1, "root", "Root");
        tax.add_node(2, 1, "superkingdom", "Bacteria");
        tax.add_node(3, 2, "phylum", "Firmicutes");
        tax.add_node(4, 3, "class", "Bacilli");
        tax.add_node(5, 1, "superkingdom", "Archaea");
        tax.add_node(6, 5, "phylum", "Euryarchaeota");
        tax
    }

    #[test]
    fn test_taxonomy_creation() {
        let tax = create_test_taxonomy();
        assert_eq!(tax.size(), 6);
        assert_eq!(tax.get_name(2), Some("Bacteria"));
    }

    #[test]
    fn test_taxonomy_default() {
        let tax = Taxonomy::default();
        assert_eq!(tax.size(), 0);
        assert_eq!(tax.root(), 1);
    }

    #[test]
    fn test_get_node() {
        let tax = create_test_taxonomy();

        let node = tax.get_node(3).unwrap();
        assert_eq!(node.taxid, 3);
        assert_eq!(node.parent, 2);
        assert_eq!(node.rank, "phylum");
        assert_eq!(node.name, "Firmicutes");
    }

    #[test]
    fn test_get_node_not_found() {
        let tax = create_test_taxonomy();
        assert!(tax.get_node(999).is_none());
    }

    #[test]
    fn test_get_parent() {
        let tax = create_test_taxonomy();

        assert_eq!(tax.get_parent(4), 3);
        assert_eq!(tax.get_parent(3), 2);
        assert_eq!(tax.get_parent(2), 1);
        assert_eq!(tax.get_parent(1), 1); // Root's parent is itself
    }

    #[test]
    fn test_get_parent_not_found() {
        let tax = create_test_taxonomy();
        assert_eq!(tax.get_parent(999), 0);
    }

    #[test]
    fn test_get_parent_opt() {
        let tax = create_test_taxonomy();

        assert_eq!(tax.get_parent_opt(4), Some(3));
        assert_eq!(tax.get_parent_opt(999), None);
    }

    #[test]
    fn test_lca_same_taxon() {
        let tax = create_test_taxonomy();
        assert_eq!(tax.lca(3, 3), 3);
    }

    #[test]
    fn test_lca_parent_child() {
        let tax = create_test_taxonomy();

        // LCA of parent and child is parent
        assert_eq!(tax.lca(2, 3), 2);
        assert_eq!(tax.lca(3, 2), 2);
    }

    #[test]
    fn test_lca_siblings() {
        let tax = create_test_taxonomy();

        // LCA of siblings is their common parent
        assert_eq!(tax.lca(2, 5), 1); // Bacteria and Archaea -> Root
    }

    #[test]
    fn test_lca_distant_relatives() {
        let tax = create_test_taxonomy();

        // LCA of Bacilli (4) and Euryarchaeota (6) is Root (1)
        assert_eq!(tax.lca(4, 6), 1);
    }

    #[test]
    fn test_lca_many_empty() {
        let tax = create_test_taxonomy();
        assert_eq!(tax.lca_many(&[]), 1); // Returns root for empty list
    }

    #[test]
    fn test_lca_many_single() {
        let tax = create_test_taxonomy();
        assert_eq!(tax.lca_many(&[4]), 4);
    }

    #[test]
    fn test_lca_many_multiple() {
        let tax = create_test_taxonomy();

        // LCA of multiple taxa in same lineage
        assert_eq!(tax.lca_many(&[2, 3, 4]), 2);

        // LCA of taxa from different branches
        assert_eq!(tax.lca_many(&[4, 6]), 1);
    }

    #[test]
    fn test_get_rank() {
        let tax = create_test_taxonomy();

        assert_eq!(tax.get_rank(1), Some("root"));
        assert_eq!(tax.get_rank(2), Some("superkingdom"));
        assert_eq!(tax.get_rank(3), Some("phylum"));
        assert_eq!(tax.get_rank(999), None);
    }

    #[test]
    fn test_get_name_not_found() {
        let tax = create_test_taxonomy();
        assert_eq!(tax.get_name(999), None);
    }

    #[test]
    fn test_is_ancestor_same() {
        let tax = create_test_taxonomy();
        assert!(tax.is_a_ancestor_of_b(3, 3)); // A is ancestor of itself
    }

    #[test]
    fn test_is_ancestor_direct_parent() {
        let tax = create_test_taxonomy();
        assert!(tax.is_a_ancestor_of_b(2, 3)); // Bacteria is ancestor of Firmicutes
        assert!(!tax.is_a_ancestor_of_b(3, 2)); // Firmicutes is NOT ancestor of Bacteria
    }

    #[test]
    fn test_is_ancestor_distant() {
        let tax = create_test_taxonomy();
        assert!(tax.is_a_ancestor_of_b(1, 4)); // Root is ancestor of Bacilli
        assert!(tax.is_a_ancestor_of_b(2, 4)); // Bacteria is ancestor of Bacilli
        assert!(!tax.is_a_ancestor_of_b(4, 2)); // Bacilli is NOT ancestor of Bacteria
    }

    #[test]
    fn test_is_ancestor_different_branches() {
        let tax = create_test_taxonomy();
        assert!(!tax.is_a_ancestor_of_b(2, 6)); // Bacteria is NOT ancestor of Euryarchaeota
        assert!(!tax.is_a_ancestor_of_b(5, 4)); // Archaea is NOT ancestor of Bacilli
    }

    #[test]
    fn test_get_external_id() {
        let tax = create_test_taxonomy();
        // Currently returns same as internal ID
        assert_eq!(tax.get_external_id(4), 4);
    }

    #[test]
    fn test_taxonomy_node_clone() {
        let node = TaxonomyNode {
            taxid: 42,
            parent: 1,
            rank: "species".to_string(),
            name: "E. coli".to_string(),
        };
        let cloned = node.clone();
        assert_eq!(node.taxid, cloned.taxid);
        assert_eq!(node.name, cloned.name);
    }

    #[test]
    fn test_taxonomy_node_debug() {
        let node = TaxonomyNode {
            taxid: 42,
            parent: 1,
            rank: "species".to_string(),
            name: "E. coli".to_string(),
        };
        let debug = format!("{:?}", node);
        assert!(debug.contains("42"));
        assert!(debug.contains("E. coli"));
    }
}

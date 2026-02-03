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

    #[test]
    fn test_taxonomy_creation() {
        let mut tax = Taxonomy::new();
        tax.add_node(1, 1, "root", "Root");
        tax.add_node(2, 1, "superkingdom", "Bacteria");
        tax.add_node(3, 2, "phylum", "Firmicutes");

        assert_eq!(tax.size(), 3);
        assert_eq!(tax.get_name(2), Some("Bacteria"));
    }

    #[test]
    fn test_lca_same_taxon() {
        let tax = Taxonomy::new();
        let result = tax.lca(1, 1);
        assert_eq!(result, 1);
    }

    #[test]
    fn test_lca_parent_child() {
        let mut tax = Taxonomy::new();
        tax.add_node(1, 1, "root", "Root");
        tax.add_node(2, 1, "superkingdom", "Bacteria");

        let result = tax.lca(1, 2);
        assert_eq!(result, 1);
    }
}

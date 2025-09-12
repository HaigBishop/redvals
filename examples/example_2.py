"""
Example 2 - Loading Decorated Trees

Copyright (c) 2025 Haig Bishop

MIT License (see LICENSE file for details)

This example demonstrates working with previously decorated trees:
1. Load the decorated tree files (.pkl) into a RedTree object
2. Convert between GTDB and redvals node IDs
3. Access other node information
4. Compute RED distances between any two nodes
5. Get the RED distance for a taxon

"""

# Import the RedTree class
from redvals import RedTree


# GTDB Release -------------
GTDB_RELEASE = "r220"

# Decorated Tree Files -------------
# A Bio.Phylo.Newick.Tree object holding the decorated archeal GTDB phylogenetic tree
ARC_DECORATED_TREE_PATH = f"decorated_trees/ar53_{GTDB_RELEASE}_decorated.pkl"
# A Bio.Phylo.Newick.Tree object holding the decorated bacterial GTDB phylogenetic tree
BAC_DECORATED_TREE_PATH = f"decorated_trees/bac120_{GTDB_RELEASE}_decorated.pkl"

# Pre-Computed Taxon Name Mappings -------------
PRECOMPUTED_TAXON_MAPPING = f"./taxon_mappings/taxon_to_node_mapping_{GTDB_RELEASE}.pkl"


# 1. Initialise (already decorated) RedTree object -----------
# We are using the decorated trees (.pkl files) as input
red_trees = RedTree(BAC_DECORATED_TREE_PATH, ARC_DECORATED_TREE_PATH, verbose=False)

# The trees are already decorated with RED values
print("Are the trees decorated:", red_trees.is_decorated())


# 2. Convert between different node ID formats -----------
# Get the redvals ID for this GTDB ID
gtdb_id = "GB_GCA_002687935.1"
bacterial_redvals_id = red_trees.get_redvals_id(gtdb_id)
print(f"\nGTDB ID '{gtdb_id}' corresponds to redvals ID '{bacterial_redvals_id}'")

# Get the redvals ID for this GTDB ID
gtdb_id = "GB_GCA_000230955.3"
bacterial_redvals_id = red_trees.get_redvals_id(gtdb_id)
print(f"GTDB ID '{gtdb_id}' corresponds to redvals ID '{bacterial_redvals_id}'")

# Get the GTDB ID for this redvals ID
redvals_id = "bac00000001"
gtdb_id = red_trees.get_gtdb_id(redvals_id)
print(f"Redvals ID '{redvals_id}' corresponds to GTDB ID '{gtdb_id}'")

# Get the GTDB ID for this redvals ID
redvals_id = "arc00002281"
gtdb_id = red_trees.get_gtdb_id(redvals_id)
print(f"Redvals ID '{redvals_id}' corresponds to GTDB ID '{gtdb_id}'")


# 3. Get node information -----------
# Get information for a specific node
redvals_id = "bac00000001"
node_info = red_trees.get_node_info(redvals_id)
print('\n', node_info)

# Get information for a specific node
redvals_id = "GB_GCA_028725265.1"
node_info = red_trees.get_node_info(redvals_id)
print('\n', node_info)

# Get information for a specific node
redvals_id = "bac00152898"
node_info = red_trees.get_node_info(redvals_id)
print('\n', node_info)

# Get information for a specific node
redvals_id = "arc00002281"
node_info = red_trees.get_node_info(redvals_id)
print('\n', node_info)

# Get information for a specific node
redvals_id = "arc00007699"
node_info = red_trees.get_node_info(redvals_id)
print('\n', node_info)



# 4. Calculate RED distances between pairs of nodes -----------
# Two leaf nodes (close by each other)
leaf_node_1 = "GB_GCA_947502505.1"
leaf_node_2 = "RS_GCF_001186155.3"
red_distance, mrca_node_id = red_trees.dist_between_nodes(leaf_node_1, leaf_node_2)
print("\nTwo Leaf Nodes (close by each other):")
print(f"{leaf_node_1} and {leaf_node_2}:")
print(f"The RED distance between them: {red_distance:.3f}")
# Get the MRCA node information
mrca_node_info = red_trees.get_node_info(mrca_node_id)
print(f"Their MRCA node is: {mrca_node_id} with RED value: {mrca_node_info.red_value:.3f}")

# Two leaf nodes (their MRCA is the root)
leaf_node_1 = "bac00000001"
leaf_node_2 = "RS_GCF_001186155.3"
red_distance, mrca_node_id = red_trees.dist_between_nodes(leaf_node_1, leaf_node_2)
print("\nTwo Leaf Nodes (their MRCA is the root):")
print(f"{leaf_node_1} and {leaf_node_2}:")
print(f"The RED distance between them: {red_distance:.3f}")
# Get the MRCA node information
mrca_node_info = red_trees.get_node_info(mrca_node_id)
print(f"Their MRCA node is: {mrca_node_id} with RED value: {mrca_node_info.red_value:.3f}")

# A leaf and an internal node
leaf_node = "GB_GCA_018970085.1"
internal_node = "bac00149750"
red_distance, mrca_node_id = red_trees.dist_between_nodes(leaf_node, internal_node)
print(f"\nLeaf and Internal Node:")
print(f"{leaf_node} and {internal_node}:")
print(f"The RED distance between them: {red_distance:.3f}")
# Get the MRCA node information
mrca_node_info = red_trees.get_node_info(mrca_node_id)
print(f"Their MRCA node is: {mrca_node_id} with RED value: {mrca_node_info.red_value:.3f}")

# Two internal nodes (close by each other)
internal_node_1 = "arc00007710"
internal_node_2 = "arc00007700"
red_distance, mrca_node_id = red_trees.dist_between_nodes(internal_node_1, internal_node_2)
print("\nTwo Internal Nodes:")
print(f"{internal_node_1} and {internal_node_2}:")
print(f"The RED distance between them: {red_distance:.3f}")
# Get the MRCA node information
mrca_node_info = red_trees.get_node_info(mrca_node_id)
print(f"Their MRCA node is: {mrca_node_id} with RED value: {mrca_node_info.red_value:.3f}")

# Get distance between identical nodes
node_id = "bac00000420"
red_distance, mrca_node_id = red_trees.dist_between_nodes(node_id, node_id)
print("\nIdentical Nodes:")
print(f"{node_id} and {node_id}:")
print(f"The RED distance between them: {red_distance:.3f}")
# Get the MRCA node information
mrca_node_info = red_trees.get_node_info(mrca_node_id)
print(f"Their MRCA node is: {mrca_node_id} with RED value: {mrca_node_info.red_value:.3f}")



# 5. Get the RED distances in taxa -----------

# You can choose to either compute the mapping of taxa to nodes (~30 minutes) or load a pre-existing mapping
load_precomputed_mapping = True

if load_precomputed_mapping:
    red_trees.load_taxa_to_node_mapping(PRECOMPUTED_TAXON_MAPPING)
else:
    # For computation of mapping, you require this file available on GTDB website
    seqs_fasta_path = f"ssu_all_{GTDB_RELEASE}.fna" 
    red_trees.map_taxa_to_nodes(seqs_fasta_path, save_result_path=PRECOMPUTED_TAXON_MAPPING)


# Get the RED distance for a taxon
taxon_name = "d__Bacteria"
red_distance = red_trees.get_distance_in_taxon(taxon_name)
print(f"The RED distance for {taxon_name}: \t\t{red_distance:.6f}")

# Get the RED distance for a taxon
taxon_name = "p__Nitrospirota"
red_distance = red_trees.get_distance_in_taxon(taxon_name)
print(f"The RED distance for {taxon_name}: \t\t{red_distance:.6f}")

# Get the RED distance for a taxon
taxon_name = "g__Escherichia"
red_distance = red_trees.get_distance_in_taxon(taxon_name)
print(f"The RED distance for {taxon_name}: \t\t{red_distance:.6f}")

# Get the RED distance for a taxon
taxon_name = "s__Spirillospora terrae"
red_distance = red_trees.get_distance_in_taxon(taxon_name)
print(f"The RED distance for {taxon_name}: \t{red_distance:.6f}")

"""
Example 2 - Loading Decorated Trees

This example demonstrates working with previously decorated trees:
1. Load the decorated tree files (.pkl) into a RedTree object
2. Convert between GTDB and redvals node IDs
3. Access other node information
4. Compute RED distances between any two nodes

"""

# Import the RedTree class
from redvals import RedTree

# Decorated Tree Files -------------
# A Bio.Phylo.Newick.Tree object holding the decorated archeal GTDB phylogenetic tree
ARC_DECORATED_TREE_PATH = "decorated_trees/ar53_r220_decorated.pkl"
# A Bio.Phylo.Newick.Tree object holding the decorated bacterial GTDB phylogenetic tree
BAC_DECORATED_TREE_PATH = "decorated_trees/bac120_r220_decorated.pkl"


# 1. Initialise (already decorated) RedTree object -----------
# We are using the decorated trees (.pkl files) as input
red_trees = RedTree(BAC_DECORATED_TREE_PATH, ARC_DECORATED_TREE_PATH)

# The trees are already decorated with RED values
print("\nAre the trees decorated:", red_trees.is_decorated())


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
print(f"\nNode {redvals_id}:")
print(f"  redvals ID: {node_info.redvals_id}")
print(f"  GTDB ID: {node_info.gtdb_id}")
print(f"  Domain: {node_info.domain}")
print(f"  RED value: {node_info.red_value}")
print(f"  RED distance: {node_info.red_distance}")
print(f"  Is terminal node: {node_info.is_terminal}")

# Get information for a specific node
redvals_id = "bac00002898"
node_info = red_trees.get_node_info(redvals_id)
print(f"\nNode {redvals_id}:")
print(f"  redvals ID: {node_info.redvals_id}")
print(f"  GTDB ID: {node_info.gtdb_id}")
print(f"  Domain: {node_info.domain}")
print(f"  RED value: {node_info.red_value}")
print(f"  RED distance: {node_info.red_distance}")
print(f"  Is terminal node: {node_info.is_terminal}")

# Get information for a specific node
redvals_id = "GB_GCA_028725265.1"
node_info = red_trees.get_node_info(redvals_id)
print(f"\nNode {redvals_id}:")
print(f"  redvals ID: {node_info.redvals_id}")
print(f"  GTDB ID: {node_info.gtdb_id}")
print(f"  Domain: {node_info.domain}")
print(f"  RED value: {node_info.red_value}")
print(f"  RED distance: {node_info.red_distance}")
print(f"  Is terminal node: {node_info.is_terminal}")

# Get information for a specific node
redvals_id = "arc00002281"
node_info = red_trees.get_node_info(redvals_id)
print(f"\nNode {redvals_id}:")
print(f"  redvals ID: {node_info.redvals_id}")
print(f"  GTDB ID: {node_info.gtdb_id}")
print(f"  Domain: {node_info.domain}")
print(f"  RED value: {node_info.red_value}")
print(f"  RED distance: {node_info.red_distance}")
print(f"  Is terminal node: {node_info.is_terminal}")

# Get information for a specific node
redvals_id = "RS_GCF_000337735.1"
node_info = red_trees.get_node_info(redvals_id)
print(f"\nNode {redvals_id}:")
print(f"  redvals ID: {node_info.redvals_id}")
print(f"  GTDB ID: {node_info.gtdb_id}")
print(f"  Domain: {node_info.domain}")
print(f"  RED value: {node_info.red_value}")
print(f"  RED distance: {node_info.red_distance}")
print(f"  Is terminal node: {node_info.is_terminal}")


# 4. Calculate RED distances between pairs of nodes -----------
# Example node pairs
node_pairs = [
    ("bac00000001", "bac00000002"),  # Two XXX nodes
    ("bac00090342", "bac00054321"),  # Two XXX nodes
    ("arc00000147", "arc00002342"),  # Two XXX nodes
]

for node1_id, node2_id in node_pairs:
    # Calculate RED distance between nodes
    red_distance, mrca_node_id = red_trees.dist_between_nodes(node1_id, node2_id)

    # Get the MRCA node information
    mrca_node_info = red_trees.get_node_info(mrca_node_id)
    
    print(f"\nThe RED distance between {node1_id} and {node2_id} is {red_distance:.6f}")
    print(f"Their MRCA node is: {mrca_node_id} with RED value: {mrca_node_info.red_value:.6f}")


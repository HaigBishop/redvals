"""
Preprocesses GTDB phylogenetic trees. Writes them as .pkl files
 - bac120 and ar53 are done seperately
 - adds RED values to every node
 - adds ((1 - RED) * 2) values to every node (called distance_between_ancestors)
 - gives every node a clade_id attribute (a unique integer)

Decorates a tree with RED values and writes as a .pkl file.
This is preprocessing for before dataset generation.
"""

from Bio import Phylo
import pandas as pd
import time
import pickle


def decorate_tree(red_values_df, tree, start_clade_id=1, verbose=True):
    """
    This function added RED and ((1 - RED) * 2) values to all nodes in a tree...
    Amoung other things!

    1. sort rows of red values descending order
    2. loop through each row of RED values
        a. if nodes is one node, find that leaf node clade
        b. else if nodes is two nodes, find the common ancestor clade
        c. assign a red attribute to the clade with the red value
        d. assign a "distance between ancestors" attribute to the clade equal to ((1 - RED) * 2)
        e. assigns a clade ID
    3. loop through all clades in the tree and count how many have no .red attribute
        a. if more than zero, print how many have no .red attribute
    4. return tree
    """
    # Sort rows of the DataFrame by RED values in descending order
    red_values_df = red_values_df.sort_values(by="RED", ascending=False)
    
    # Initialize clade_id counter
    clade_id_counter = start_clade_id

    # Assign RED values to the tree
    total_rows = len(red_values_df)
    checkpoint = max(1, total_rows // 2000)  # Number of checkpoints
    start_time = time.time()  # Record the start time

    for i, (_, row) in enumerate(red_values_df.iterrows(), 1):
        node = row["nodes"]
        red_value = row["RED"]

        # If the node is a single leaf node
        if '|' not in node:
            clade = tree.find_any(name=node)
            if clade:
                # Leaf clade !
                clade.red = red_value
                clade.distance_between_ancestors = (1 - red_value) * 2
                clade.clade_id = clade_id_counter
                clade_id_counter += 1
        
        # If the node represents 2 nodes (MRCA of two nodes)
        else:
            two_nodes = node.split('|')
            common_ancestor = tree.common_ancestor(two_nodes)
            if common_ancestor:
                # Proper clade !
                if not hasattr(common_ancestor, "clade_id"):
                    common_ancestor.red = red_value
                    common_ancestor.distance_between_ancestors = (1 - red_value) * 2
                    common_ancestor.clade_id = clade_id_counter
                    clade_id_counter += 1

        # Print progress every X%
        if verbose and i % checkpoint == 0:
            progress = (i / total_rows) * 100
            elapsed_time = time.time() - start_time
            estimated_total_time = (elapsed_time / i) * total_rows
            estimated_remaining_time = (estimated_total_time - elapsed_time) / 60
            print(f"Progress: {progress:.2f}%, Estimated time remaining: {estimated_remaining_time:.2f} minutes")


    # Count and report clades without the RED attribute
    unannotated_clades = [clade for clade in tree.find_clades() if not hasattr(clade, 'red')]
    if unannotated_clades:
        print(unannotated_clades)
        print(f"WARNING! {len(unannotated_clades)} clades have no RED attribute.")

    return tree

def add_lineage_to_leaves(tree):
    """
    Adds a .lineage attribute to every leaf node, which is a tuple of clade_ids
    from the root to that leaf node (including the leaf node).
    """
    def traverse(clade, lineage):
        # Append the current clade's clade_id to the lineage
        lineage = lineage + (clade.clade_id,)
        if clade.is_terminal():
            # Assign the lineage tuple to the leaf node
            clade.lineage = lineage
        else:
            # Recurse on child clades
            for child in clade.clades:
                traverse(child, lineage)
    # Start traversal from the root of the tree
    traverse(tree.root, ())


def get_most_recent_common_ancestor(tree, leaf_1, leaf_2, clade_dictionary):
    """
    Returns the clade object of the MRCA of the two leaf nodes in the tree.
    Relies on the leaf nodes having lineage attributes 
    And relies on all nodes having clade_id attributes which are found in the clade_dictionary
    """
    lineage1 = leaf_1.lineage
    lineage2 = leaf_2.lineage

    # Find the minimum length of the two lineages
    min_length = min(len(lineage1), len(lineage2))

    # Initialize MRCA clade_id
    mrca_clade_id = None

    # Iterate over both lineages to find the deepest common clade_id
    for i in range(min_length):
        if lineage1[i] == lineage2[i]:
            mrca_clade_id = lineage1[i]
        else:
            break

    if mrca_clade_id is not None:
        return clade_dictionary[mrca_clade_id]
    else:
        return None  # No common ancestor found (should not happen in a connected tree)


def get_distance_between_leaves(tree, leaf_1, leaf_2, clade_dictionary):
    """
    Returns the distance_between_ancestors of the most recent common ancestor (MRCA)
    of the two leaf nodes in the tree.

    Parameters:
        tree: Phylo tree object with nodes that have clade_id attributes.
        leaf_1, leaf_2: The two leaf nodes with .lineage attributes.
        clade_dictionary: Dictionary mapping clade_id to node objects.

    Returns:
        float: The distance_between_ancestors of the MRCA, or None if no MRCA is found.
    """
    mrca = get_most_recent_common_ancestor(tree, leaf_1, leaf_2, clade_dictionary)
    return mrca.distance_between_ancestors if mrca else None



def make_clade_id_dictionary(tree):
    """
    Makes a dictionary which maps all clade_id of all nodes to the node object.
    WARNING: This function must be used to make the dictionary in the python session that it is used in... pretty sure.

    Parameters:
        tree: Phylo tree object with nodes that have a clade_id attribute.

    Returns:
        dict: A dictionary mapping clade_id (int) to the corresponding node object.
    """
    return {clade.clade_id: clade for clade in tree.find_clades() if hasattr(clade, 'clade_id')}



if __name__ == "__main__":
    # Pick a domain foo
    DOMAIN_NAME = 'bac120'      # bac120 or ar53

    # Start clade IDs at either 0000000 or 100000-
    start_clade_id = 0 if DOMAIN_NAME == 'bac120' else 100000

    # Load the Newick tree file of the OG tree
    tree_file = "./trees/" + DOMAIN_NAME + "_r220.tree"
    og_tree = Phylo.read(tree_file, "newick")

    print(type(og_tree))
    input()


    # Load the TSV file of the red values
    red_values_tsv_file = "D:/16S_databases/GTDBtk/mrca_red/gtdbtk_r220_" + DOMAIN_NAME + ".tsv"
    red_values_df = pd.read_csv(red_values_tsv_file, sep="\t", header=None, names=["nodes", "RED"])


    # Add RED values to the tree (and distance_between_ancestors and clade_id)
    redded_tree = decorate_tree(red_values_df, og_tree, start_clade_id=start_clade_id)

    # Add lineages to the tree
    add_lineage_to_leaves(redded_tree)


    # Write object to file
    tree_object_file = "D:/16S_databases/GTDBtk/" + DOMAIN_NAME + "_r220_decorated.pkl"
    with open(tree_object_file, "wb") as f:
        pickle.dump(redded_tree, f)
    print(f"Decorated tree object saved to {tree_object_file}")

    # (LOAD WITH):
    # with open(tree_object_file, "rb") as f:
    #     loaded_tree = pickle.load(f)


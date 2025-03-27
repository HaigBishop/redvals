"""redvals: A tool for obtaining RED (Relative Evolutionary Divergence) values from GTDB phylogenetic trees

This module provides functionality for handling phylogenetic trees from the Genome Taxonomy Database (GTDB) 
and calculating Relative Evolutionary Divergence (RED) values. It supports both bacterial and archaeal trees
in a single RedTree object.

The module contains two main classes:
    - RedTree: For handling phylogenetic trees with RED values
    - NodeInfo: For storing information about nodes in the trees

Features:
    - Load and parse GTDB phylogenetic trees in Newick format
    - Decorate trees with RED values from TSV files
    - Save and load decorated trees using pickle format
    - Convert between different node identifier formats (GTDB IDs and redvals IDs)
    - Calculate RED distances between any two nodes
    - Retrieve detailed information about nodes

"""

import os
from Bio import Phylo
import pickle
from tqdm import tqdm
import pandas as pd


# DEFAULT FILE PATHS ========================================================
# Original Newick Tree Files -------------
# The original Newick format bacterial GTDB phylogenetic tree file
BAC_TREE_PATH = "trees/bac120_r220.tree"
# The original Newick format archeal GTDB phylogenetic tree file
ARC_TREE_PATH = "trees/ar53_r220.tree"

# RED Values TSV Files -------------
# The TSV file containing RED values for nodes in the bacterial GTDB phylogenetic tree
BAC_RED_VALUES_PATH = "red_values/gtdbtk_r220_bac120.tsv"
# The TSV file containing RED values for nodes in the archeal GTDB phylogenetic tree
ARC_RED_VALUES_PATH = "red_values/gtdbtk_r220_ar53.tsv"

# Decorated Tree Files -------------
# A Bio.Phylo.Newick.Tree object holding the decorated bacterial GTDB phylogenetic tree
BAC_DECORATED_TREE_PATH = "decorated_trees/bac120_r220_decorated.pkl"
# A Bio.Phylo.Newick.Tree object holding the decorated archeal GTDB phylogenetic tree
ARC_DECORATED_TREE_PATH = "decorated_trees/ar53_r220_decorated.pkl"


class RedTree:
    def __init__(self, bac_tree_path=None, arc_tree_path=None, verbose=False):
        """
        Initialise the RedTree object.
        - bac_tree_path and arc_tree_path can be either:
          - paths to Newick format .tree files
          - paths to .pkl files containing decorated Bio.Phylo.Newick.Tree objects derived from RedTree.write_trees()
        """
        # Set the verbose attribute
        self.verbose = verbose
        
        # If no paths are provided, use the default paths
        if bac_tree_path is None:
            bac_tree_path = BAC_TREE_PATH
        if arc_tree_path is None:
            arc_tree_path = ARC_TREE_PATH

        # Enforce that both are the same format
        if bac_tree_path.endswith(".tree") != arc_tree_path.endswith(".tree") or bac_tree_path.endswith(".pkl") != arc_tree_path.endswith(".pkl"):
            raise ValueError("Both trees must be in the same format")
        
        # If the trees are .tree files
        if bac_tree_path.endswith(".tree"):
            self.bac_tree = Phylo.read(bac_tree_path, "newick")
            self.arc_tree = Phylo.read(arc_tree_path, "newick")
            print(f"\nLoaded Newick format trees from {bac_tree_path} and {arc_tree_path}")

            # The trees are not decorated
            self.bac_tree.is_decorated = False
            self.arc_tree.is_decorated = False
            # Assign redvals IDs to all the nodes in the trees
            self.assign_redvals_ids()
            # Make a dictionaries to do with node IDs
            self.make_node_id_dicts()
            
            print("Successfully initialised RedTree object using Newick format trees...")
            print("Note: The the nodes have IDs, but they are not yet decorated with RED values.\n")

        # If the trees are .pkl files
        elif bac_tree_path.endswith(".pkl"):
            self.bac_tree = pickle.load(open(bac_tree_path, "rb"))
            self.arc_tree = pickle.load(open(arc_tree_path, "rb"))
            print(f"\nLoaded decorated .pkl files from {bac_tree_path} and {arc_tree_path}")
            # Check if the trees are valid Bio.Phylo.Newick.Tree objects
            if not isinstance(self.bac_tree, Phylo.Newick.Tree) or not isinstance(self.arc_tree, Phylo.Newick.Tree):
                raise ValueError("The trees must be valid Bio.Phylo.Newick.Tree objects")
            # Check if the trees are decorated
            decorated = self.check_decorated()
            # Set the is_decorated attribute of the trees
            self.bac_tree.is_decorated = decorated
            self.arc_tree.is_decorated = decorated
            # If the trees are not decorated
            if not decorated:
                raise ValueError("The provided Bio.Phylo.Newick.Tree .pkl files are not decorated")
            # Make a dictionaries to do with node IDs
            self.make_node_id_dicts()
            print("Successfully initialised RedTree object using decorated .pkl files...")
            print("The the nodes are decorated with IDs and RED values.")
            # Make a dictionary for NodeInfo objects
            self.make_node_info_dict()
            print("Successfully created NodeInfo mappings for the IDs of all nodes.\n")

    def assign_redvals_ids(self):
        """
        Assign redvals IDs to all the nodes in the trees.
         - The redvals IDs are assigned deterministically according to the Newick format tree file
         - IDs follow the format "bac00000001", "bac00000002", "arc00000001", "arc00000002", etc.
         - all nodes are given a .redvals_id attribute
        """
        # Helper function to assign IDs to nodes in a tree
        def assign_ids(tree, prefix):
            counter = 1
            # Assign IDs to terminal nodes first
            for node in tree.get_terminals():
                node.redvals_id = f"{prefix}{counter:08d}"
                counter += 1
            # Then assign IDs to internal nodes
            for node in tree.get_nonterminals():
                node.redvals_id = f"{prefix}{counter:08d}"
                counter += 1
            return counter - 1

        # Assign IDs to bacterial tree nodes
        bac_count = assign_ids(self.bac_tree, "bac")
        
        # Assign IDs to archaeal tree nodes
        arc_count = assign_ids(self.arc_tree, "arc")
        print(f"Assigned IDs to all {bac_count} bacterial and all {arc_count} archaeal nodes")

        return bac_count, arc_count

    def make_node_id_dicts(self):
        """
        Make dictionaries to do with node IDs.
         - self.get_node_dict: a dictionary mapping node IDs to nodes
         - self.get_redvals_id_dict: a dictionary mapping GTDB IDs to redvals IDs
         - self.get_gtdb_id_dict: a dictionary mapping redvals IDs to GTDB IDs
        """
        # Create a set to check for duplicate and invalid names
        name_set = set(['', 'None', None])
        # Get all nodes in both trees
        all_nodes = list(self.bac_tree.get_terminals()) + list(self.bac_tree.get_nonterminals()) + list(self.arc_tree.get_terminals()) + list(self.arc_tree.get_nonterminals())

        # Populate the dictionaries
        self.get_node_dict = {}
        self.get_redvals_id_dict = {}
        self.get_gtdb_id_dict = {}

        # For every node in both trees
        for node in all_nodes:
            # Have GTDB ID? e.g. 'GB_GCA_018399855.1'
            if hasattr(node, 'name'):
                if node.name not in name_set:
                    name_set.add(node.name)
                    self.get_node_dict[node.name] = node
                    self.get_redvals_id_dict[node.name] = None # (this may be replaced in following code)
            # Have redvals ID? e.g. 'bac00000001'
            if hasattr(node, 'redvals_id'):
                self.get_node_dict[node.redvals_id] = node
                self.get_gtdb_id_dict[node.redvals_id] = None # (this may be replaced in following code)
            else:
                raise ValueError(f"Node '{node.name}' does not have a redvals ID. This should never happen.")
            # Have both GTDB ID and redvals ID?
            if hasattr(node, 'name') and hasattr(node, 'redvals_id'):
                self.get_gtdb_id_dict[node.redvals_id] = node.name
                self.get_redvals_id_dict[node.name] = node.redvals_id
        
        print(f"Created mappings for the IDs of {len(self.get_node_dict)} nodes")

    def make_node_info_dict(self):
        """
        Make a dictionary to do with node IDs.
         - self.get_node_info_dict: a dictionary mapping node IDs to NodeInfo objects
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated to make the node info dictionary")
        
        # Populate the dictionary
        self.get_node_info_dict = {}

        def add_to_node_info_dict(nodes, domain, is_terminal):
            # For every node
            for node in nodes:
                # Check if the node is in the gtdb_id_dict
                if node.redvals_id not in self.get_gtdb_id_dict:
                    raise ValueError(f"Node '{node.redvals_id}' not found in the tree")
                # Get the GTDB ID (may be None)
                gtdb_id = self.get_gtdb_id_dict[node.redvals_id]
                node_info = NodeInfo(domain, node.red_value, gtdb_id, node.redvals_id, is_terminal, node.red_distance)
                # Store the node info in the dictionary
                self.get_node_info_dict[node.redvals_id] = node_info
                if gtdb_id is not None:
                    self.get_node_info_dict[gtdb_id] = node_info

        # Add all archaeal leaf nodes
        add_to_node_info_dict(self.arc_tree.get_terminals(), 'arc', True)
        # Add all bacterial leaf nodes
        add_to_node_info_dict(self.bac_tree.get_terminals(), 'bac', True)
        # Add all archaeal internal nodes
        add_to_node_info_dict(self.arc_tree.get_nonterminals(), 'arc', False)
        # Add all bacterial internal nodes
        add_to_node_info_dict(self.bac_tree.get_nonterminals(), 'bac', False)
        
        print(f"Created NodeInfo mappings for the IDs of {len(self.get_node_dict)} nodes")

    def get_node(self, node_id):
        """
        Given a node ID, return the node.
        e.g. "bac00000001" -> Bio.Phylo.Newick.Clade object
        e.g. "GB_GCA_018399855.1" -> Bio.Phylo.Newick.Clade object
        """
        
        if node_id not in self.get_node_dict:
            raise KeyError(f"Node ID '{node_id}' not found in the tree")
            
        return self.get_node_dict[node_id]
    
    def get_gtdb_id(self, redvals_id):
        """
        Given a redvals ID, return the GTDB ID.
        e.g. "bac00000001" -> "GB_GCA_018399855.1"
        """
            
        if redvals_id not in self.get_gtdb_id_dict:
            raise KeyError(f"Redvals ID '{redvals_id}' not found or does not have a corresponding GTDB ID")
        
        gtdb_id = self.get_gtdb_id_dict[redvals_id]

        if gtdb_id is None:
            print(f"Redvals ID '{redvals_id}' does not have a corresponding GTDB ID. Likely an unlablled internal node.")
            
        return gtdb_id

    def get_redvals_id(self, gtdb_id):
        """
        Given a GTDB ID, return the redvals ID.
        e.g. "GB_GCA_018399855.1" -> "bac00000001"
        """
            
        if gtdb_id not in self.get_redvals_id_dict:
            raise KeyError(f"GTDB ID '{gtdb_id}' not found in the tree")
        
        redvals_id = self.get_redvals_id_dict[gtdb_id]

        if redvals_id is None:
            raise ValueError(f"GTDB ID '{gtdb_id}' does not have a corresponding redvals ID. This should never happen.")
            
        return redvals_id

    def get_node_info(self, node_id):
        """
        Given a node ID, return the node information.
        - The node IDs can be either:
            - IDs assigned by redvals - e.g. "bac00001342" or "arc00000281"
            - IDs derived from the GTDB phylogenetic trees - e.g. "GB_GCA_002687935.1" or "GB_GCA_000230955.3"
        """
        if node_id not in self.get_node_info_dict:
            raise KeyError(f"Node ID '{node_id}' not found in the tree")
        
        return self.get_node_info_dict[node_id]
    
    def is_decorated(self):
        """
        Return True if the tree is decorated (i.e. has RED values), False otherwise.
        """
        # Check if the trees have the attribute is_decorated
        have_attrs = hasattr(self.bac_tree, "is_decorated") and hasattr(self.arc_tree, "is_decorated")
        if not have_attrs:
            # Check if the trees are decorated
            decorated = self.check_decorated()
            # Set the is_decorated attribute of the trees
            self.bac_tree.is_decorated = decorated
            self.arc_tree.is_decorated = decorated
        # Check if the trees are decorated
        is_decorated = self.bac_tree.is_decorated and self.arc_tree.is_decorated
        return is_decorated
    
    def check_decorated(self):
        """
        Check if the trees are decorated. Unlike is_decorated(), this function will actually check the attributes of all nodes in the trees.

        To be decorated, a tree must:
        - have a .red_value attribute for every node
        - have a .red_distance attribute for every node
        - have a .redvals_id attribute for every node

        Returns True if the trees are decorated.
        """
        is_decorated = True
        
        # Get all nodes from both trees
        all_nodes = list(self.bac_tree.get_terminals()) + list(self.bac_tree.get_nonterminals()) + \
                    list(self.arc_tree.get_terminals()) + list(self.arc_tree.get_nonterminals())
        
        # Check all required attributes for each node
        for node in all_nodes:
            if not hasattr(node, 'red_value') or not hasattr(node, 'red_distance') or not hasattr(node, 'redvals_id'):
                is_decorated = False
                break

        return is_decorated
    
    def decorate_from_tsv(self, bac_red_values_path=None, arc_red_values_path=None, progress_bar=True):
        """
        Decorate the trees (with RED values) from TSV files containing RED values.
        """
        # If no paths are provided, use the default paths
        if bac_red_values_path is None:
            bac_red_values_path = BAC_RED_VALUES_PATH
        if arc_red_values_path is None:
            arc_red_values_path = ARC_RED_VALUES_PATH
        # Check if the trees are decorated
        if self.is_decorated():
            raise ValueError("The trees are already decorated")
        # Check the files exist
        if not os.path.exists(bac_red_values_path):
            raise FileNotFoundError(f"Bacterial RED values file not found at {bac_red_values_path}")
        if not os.path.exists(arc_red_values_path):
            raise FileNotFoundError(f"Archaeal RED values file not found at {arc_red_values_path}")

        def decorate_tree(red_values_path, tree, progress_bar=False):
            # Read the TSV file
            red_values_df = pd.read_csv(red_values_path, sep="\t", header=None, names=["nodes", "RED"])
            # Randomly shuffle rows of the DataFrame
            red_values_df = red_values_df.sample(frac=1).reset_index(drop=True)
            
            # Setup progress bar if requested
            if progress_bar:
                total_rows = len(red_values_df)
                pbar = tqdm(total=total_rows, desc="Decorating tree", unit="nodes")
            
            # For every RED value in the DataFrame
            for i, (_, row) in enumerate(red_values_df.iterrows(), 1):
                node_column = row["nodes"]
                red_value = row["RED"]

                # If the node column has only one leaf node
                if '|' not in node_column:
                    # Check that the RED is 1.0
                    if red_value != 1.0:
                        raise ValueError(f"RED value for leaf node '{node_column}' is not 1.0")
                    # Get the node
                    leaf_node = self.get_node(node_column)
                    # Assign the RED value
                    leaf_node.red_value = red_value
                    # Assign the RED distance
                    leaf_node.red_distance = (1 - red_value) * 2
                
                # If the node column has two leaf nodes
                else:
                    # Get the two leaf nodes
                    leaf_nodes = node_column.split('|')
                    # Get the nodes
                    leaf_node_1 = self.get_node(leaf_nodes[0])
                    leaf_node_2 = self.get_node(leaf_nodes[1])
                    # Get the MRCA
                    mrca_node = tree.common_ancestor(leaf_node_1, leaf_node_2)
                    # Check that the MRCA is not None
                    if mrca_node is None:
                        raise ValueError(f"MRCA for pair of leaf nodes '{node_column}' is not found")
                    # Assign the RED value
                    mrca_node.red_value = red_value
                    # Assign the RED distance
                    mrca_node.red_distance = (1 - red_value) * 2
                
                # Update progress bar if enabled
                if progress_bar:
                    pbar.update(1)
                    # Add percentage complete to description
                    if i % 100 == 0 or i == total_rows:
                        percent = (i / total_rows) * 100
                        pbar.set_description(f"Decorating... ({percent:.1f}%)")
            
            # Close progress bar
            if progress_bar:
                pbar.close()
        
        # Decorate both trees
        print(f"Decorating the archaeal tree with RED values from {arc_red_values_path}...")
        decorate_tree(arc_red_values_path, self.arc_tree, progress_bar=progress_bar)
        print(f"Decorating the bacterial tree with RED values from {bac_red_values_path}...")
        decorate_tree(bac_red_values_path, self.bac_tree, progress_bar=progress_bar)

        # Check that ALL nodes are decorated
        if not self.check_decorated():
            raise ValueError("The trees were failed to be decorated properly")
        else:
            print("Trees successfully decorated with RED values")
            # Set the is_decorated attribute of the trees
            self.bac_tree.is_decorated = True
            self.arc_tree.is_decorated = True
            # Make a dictionary for NodeInfo objects
            self.make_node_info_dict()
            print()
    
    def decorate_from_calc(self):
        """
        Decorate the trees (with RED values) by calculating RED values.
        """
        raise NotImplementedError("This function is not yet implemented")
        # # Check if the trees are decorated
        # if self.is_decorated():
        #     raise ValueError("The trees are already decorated")

    def get_mrca_node_from_ids(self, node_1_id, node_2_id):
        """
        Given two node IDs, return their MRCA.
        - The two nodes could be leaf nodes, internal nodes or one of each
        - The node IDs can be either:
            - IDs assigned by redvals - e.g. "bac00001342" or "arc00000281"
            - IDs derived from the GTDB phylogenetic trees - e.g. "GB_GCA_002687935.1" or "GB_GCA_000230955.3"
        """
        # Get the nodes
        node_1 = self.get_node(node_1_id)
        node_2 = self.get_node(node_2_id)
        # Get the MRCA
        return self.get_mrca_node(node_1, node_2)

    def get_mrca_node(self, node_1, node_2):
        """
        Given two nodes, return their MRCA.
        - The two nodes could be leaf nodes, internal nodes or one of each
        """
        # Get the domain of the nodes
        node_1_domain = node_1.redvals_id[0:3]
        node_2_domain = node_2.redvals_id[0:3]
        # Ensure the domains are the same
        if node_1_domain != node_2_domain:
            raise ValueError(f"The two nodes '{node_1.redvals_id}' and '{node_2.redvals_id}' are not from the same domain")
        # Get the correct tree
        tree = self.bac_tree if node_1_domain == "bac" else self.arc_tree
        # Get the MRCA
        mrca_node = tree.common_ancestor(node_1, node_2)
        # Check that the MRCA is not None
        if mrca_node is None:
            print(f"MRCA for nodes '{node_1.redvals_id}' and '{node_2.redvals_id}' is not found. Returning None.")
        return mrca_node

    def dist_between_nodes(self, node_1_id, node_2_id):
        """
        Given two node IDs, return their RED distance and the node ID of their MRCA, taking into account internal nodes.
        - The two nodes could be leaf nodes, internal nodes or one of each
        - The node IDs can be either:
            - IDs assigned by redvals - e.g. "bac00001342" or "arc00000281"
            - IDs derived from the GTDB phylogenetic trees - e.g. "GB_GCA_002687935.1" or "GB_GCA_000230955.3"
        - MRCA node ID is returned as a redval-assigned ID (e.g. "bac00001342" or "arc00000281")
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated")
        # Get the nodes
        node_1 = self.get_node(node_1_id)
        node_2 = self.get_node(node_2_id)
        # Get the MRCA
        mrca_node = self.get_mrca_node(node_1, node_2)
        # Get the RED distance (between ancestor nodes)
        mrca_red_distance = mrca_node.red_distance
        # Subtract the RED values of the two nodes (incase they are internal nodes)
        red_distance_between_nodes = mrca_red_distance - node_1.red_distance - node_2.red_distance
        return red_distance_between_nodes, mrca_node.redvals_id
    
    def write_decorated_trees(self, bac_tree_path=None, arc_tree_path=None):
        """
        Write the trees to file.
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated")
        # If the paths are not provided, use the default paths
        if bac_tree_path is None:
            bac_tree_path = BAC_DECORATED_TREE_PATH
        if arc_tree_path is None:
            arc_tree_path = ARC_DECORATED_TREE_PATH
        # Ensure the files don't already exist
        if os.path.exists(bac_tree_path):
            raise FileExistsError(f"Bacterial decorated tree file already exists at {bac_tree_path}")
        if os.path.exists(arc_tree_path):
            raise FileExistsError(f"Archaeal decorated tree file already exists at {arc_tree_path}")
        # If the directory doesn't exist, create it
        if not os.path.exists(os.path.dirname(bac_tree_path)):
            os.makedirs(os.path.dirname(bac_tree_path))
        if not os.path.exists(os.path.dirname(arc_tree_path)):
            os.makedirs(os.path.dirname(arc_tree_path))
        # Write the trees to file (pickle format)
        pickle.dump(self.bac_tree, open(bac_tree_path, "wb"))
        pickle.dump(self.arc_tree, open(arc_tree_path, "wb"))
        print(f"Wrote decorated trees (as .pkl) to {bac_tree_path} and {arc_tree_path}")
    
    def get_node_domain(self, node_id):
        """
        Given a node ID, return the node domain.
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated")
        # Check if the node ID is valid
        if node_id not in self.get_node_info_dict:
            raise KeyError(f"Node ID '{node_id}' not found in the tree")
        # Return the node domain
        return self.get_node_info_dict[node_id].domain
    
    def get_node_gtdb_id(self, node_id):
        """
        Given a node ID, return the GTDB ID.
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated")
        # Check if the node ID is valid
        if node_id not in self.get_node_info_dict:
            raise KeyError(f"Node ID '{node_id}' not found in the tree")
        return self.get_node_info_dict[node_id].gtdb_id
    
    def get_node_redvals_id(self, node_id):
        """
        Given a node ID, return the redvals ID.
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated")
        # Check if the node ID is valid
        if node_id not in self.get_node_info_dict:
            raise KeyError(f"Node ID '{node_id}' not found in the tree")
        return self.get_node_info_dict[node_id].redvals_id
    
    def get_node_red_value(self, node_id):
        """
        Given a node ID, return the RED value.
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated")
        # Check if the node ID is valid
        if node_id not in self.get_node_info_dict:
            raise KeyError(f"Node ID '{node_id}' not found in the tree")
        return self.get_node_info_dict[node_id].red_value
    
    def is_node_terminal(self, node_id):
        """
        Given a node ID, return True if the node is a terminal node, False otherwise.
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated")
        # Check if the node ID is valid
        if node_id not in self.get_node_info_dict:
            raise KeyError(f"Node ID '{node_id}' not found in the tree")
        return self.get_node_info_dict[node_id].is_terminal
    
    def get_node_red_distance(self, node_id):
        """
        Given a node ID, return the RED distance.
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated")
        # Check if the node ID is valid
        if node_id not in self.get_node_info_dict:
            raise KeyError(f"Node ID '{node_id}' not found in the tree")
        return self.get_node_info_dict[node_id].red_distance

class NodeInfo:
    def __init__(self, domain, red_value, gtdb_id, redvals_id, is_terminal, red_distance=None):
        """
        Initialise the NodeInfo object.
        """
        self.domain = domain            # "bac" or "arc"
        self.gtdb_id = gtdb_id          # e.g. "GB_GCA_002687935.1" or None
        self.redvals_id = redvals_id    # e.g. "bac00001342"
        self.red_value = red_value      # e.g. 0.95
        self.is_terminal = is_terminal  # True or False
        self.red_distance = (1 - red_value) * 2 if red_distance is None else red_distance

    def __str__(self):
        return f"NodeInfo(domain={self.domain}, gtdb_id={self.gtdb_id}, redvals_id={self.redvals_id}, red_value={self.red_value}, is_terminal={self.is_terminal}, red_distance={self.red_distance})"
    
    def __repr__(self):
        return self.__str__()

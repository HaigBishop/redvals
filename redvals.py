"""redvals: A tool for obtaining RED (Relative Evolutionary Divergence) values from GTDB phylogenetic trees

Copyright (c) 2025 Haig Bishop

MIT License (see LICENSE file for details)

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
from collections import defaultdict
import random


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
    """
    A class representing both the bacterial and archaeal GTDB phylogenetic trees, with RED values.

    Decoration of the trees is done using the decorate_from_tsv() method. In this process, every node in both trees is assigned the following attributes:
      RED Value (red_value):
        - For any nodeX, nodeX.red_value is the RED value of the nodeX, directly derived from the TSV file.
      RED Distance (red_distance):
        - For any nodeX, nodeX.red_distance is the RED distance between any two descendants of nodeX for which their MRCA is nodeX, calculated as 2 * (1 - nodeX.red_value).
    """
    def __init__(self, bac_tree_path=None, arc_tree_path=None, verbose=False):
        """
        Initialise the RedTree object.
        - bac_tree_path and arc_tree_path can be either:
          - paths to Newick format .tree files
          - paths to .pkl files containing decorated Bio.Phylo.Newick.Tree objects derived from RedTree.write_decorated_trees()
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
            if self.verbose:
                print(f"\nLoaded Newick format trees from {bac_tree_path} and {arc_tree_path}")

            # The trees are not decorated
            self.bac_tree.is_decorated = False
            self.arc_tree.is_decorated = False
            # Assign redvals IDs to all the nodes in the trees
            self.assign_redvals_ids()
            # Make a dictionaries to do with node IDs
            self.make_node_id_dicts()
            
            if self.verbose:
                print("Successfully initialised RedTree object using Newick format trees...")
                print("Note: The the nodes have IDs, but they are not yet decorated with RED values.\n")

        # If the trees are .pkl files
        elif bac_tree_path.endswith(".pkl"):
            self.bac_tree = pickle.load(open(bac_tree_path, "rb"))
            self.arc_tree = pickle.load(open(arc_tree_path, "rb"))
            if self.verbose:
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
            if self.verbose:
                print("Successfully initialised RedTree object using decorated .pkl files...")
                print("The the nodes are decorated with IDs and RED values.")
            # Make a dictionary for NodeInfo objects
            self.make_node_info_dict()
            if self.verbose:
                print("Successfully created NodeInfo mappings for the IDs of all nodes.\n")
        
        # Initialise node_from_taxon_name_dict as None
        self.node_from_taxon_name_dict = None

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
        
        if self.verbose:
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
        
        if self.verbose:
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
        
        if self.verbose:
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
    
    def get_nodes(self, domain='both', node_type='both'):
        """
        Return all nodes in the tree (with optional filters).

        Args:
            domain (str): Can be 'both', 'bac' or 'arc'.
            node_type (str): Can be 'both', 'terminal'/'leaf' or 'nonterminal'/'internal'.

        Returns:
            list: A list of nodes.
        """
        # Check for valid args
        if domain not in ['both', 'bac', 'arc']:
            raise ValueError(f"Invalid domain: '{domain}'")
        if node_type not in ['both', 'terminal', 'nonterminal', 'leaf', 'internal']:
            raise ValueError(f"Invalid node type: '{node_type}'")

        # Get the correct tree
        node_list = []
        
        # Handle terminal/leaf nodes
        if node_type in ['both', 'terminal', 'leaf']:
            if domain in ['both', 'bac']:
                node_list.extend(self.bac_tree.get_terminals())
            if domain in ['both', 'arc']:
                node_list.extend(self.arc_tree.get_terminals())
                
        # Handle nonterminal/internal nodes
        if node_type in ['both', 'nonterminal', 'internal']:
            if domain in ['both', 'bac']:
                node_list.extend(self.bac_tree.get_nonterminals())
            if domain in ['both', 'arc']:
                node_list.extend(self.arc_tree.get_nonterminals())

        # Return the list of nodes
        return node_list
    
    def get_node_ids(self, id_type='redvals', domain='both', node_type='both', keep_none=False):
        """
        Return all node IDs in the tree (with optional filters).

        Args:
            id_type (str): Can be 'redvals' or 'gtdb'.
            domain (str): Can be 'both', 'bac' or 'arc'.
            node_type (str): Can be 'both', 'terminal'/'leaf' or 'nonterminal'/'internal'.
            keep_none (bool): If True, keep None IDs.

        Returns:
            list: A list of node IDs.
        """
        # Get the nodes
        nodes = self.get_nodes(domain, node_type)
        # Get the node IDs
        if id_type == 'redvals':
            return [node.redvals_id for node in nodes]
        elif id_type == 'gtdb':
            if keep_none:
                return [node.name for node in nodes]
            else:
                return [node.name for node in nodes if node.name not in ['', 'None', None]]
        else:
            raise ValueError(f"Invalid ID type: '{id_type}'")
    
    def get_gtdb_id(self, redvals_id):
        """
        Given a redvals ID, return the GTDB ID.
        e.g. "bac00000001" -> "GB_GCA_018399855.1"
        """
            
        if redvals_id not in self.get_gtdb_id_dict:
            raise KeyError(f"Redvals ID '{redvals_id}' not found or does not have a corresponding GTDB ID")
        
        gtdb_id = self.get_gtdb_id_dict[redvals_id]

        if gtdb_id is None and self.verbose:
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
        if mrca_node is None and self.verbose:
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
        # Calculate the RED distance
        red_distance_between_nodes = (node_1.red_value - mrca_node.red_value) + (node_2.red_value - mrca_node.red_value)
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
        if self.verbose:
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

    def map_taxa_to_nodes(self, seqs_fasta_path, save_result_path=None):
        """
        Assigns taxon names to nodes and creates a mapping from taxon names to nodes.

        Reads a FASTA file where headers contain GTDB IDs and taxonomic lineages.
        It assigns a 'taxon_names' attribute (a list of strings, e.g., ['g__Escherichia', 'f__Enterobacteriaceae']) 
        to the node representing the Most Recent Common Ancestor (MRCA) of all leaves belonging 
        to that taxon found in the FASTA.
        It also populates the `self.node_from_taxon_name_dict` dictionary, mapping
        taxon names directly to their corresponding MRCA node objects.
        
        Args:
            seqs_fasta_path (str): Path to the FASTA file with sequence information and taxonomy.
                                   Headers should be formatted like: >{GTDB_ID}~{lineage} [...]
            save_result_path (str): Path to save the mapping to. (e.g. "./taxon_mapping/taxon_to_node_mapping.pkl")

        Example seqs_fasta_path (e.g. "D:/16S_databases/ssu_all_r220.fna"):
            >RS_GCF_018344175.1~NZ_JAAMUS010000074.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli [location=16..1459] [ssu_len=1444] [contig_len=1463]
            ACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATC
            >RS_GCF_001246675.1~NZ_CXGB01000177.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli [location=5822..7324] [ssu_len=1503] [contig_len=7713]
            ATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGG
            >RS_GCF_000335255.2~NZ_AOEB01000154.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli [location=2..1012] [ssu_len=1011] [contig_len=27380]
            TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGA
        (These files are available on the GTDB website - e.g. https://gtdb.ecogenomic.org/downloads  ->  /public/gtdb/data/releases/release220/220.0/genomic_files_all)

        This method is long and complex, due to optimisations which drastically speed up calls to common_ancestor.
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated before mapping taxa to nodes.")
        
        # Check if the FASTA file exists
        if not os.path.exists(seqs_fasta_path):
            raise FileNotFoundError(f"Sequence FASTA file not found at {seqs_fasta_path}")

        def parse_fasta_line(line):
            if not line.startswith('>'):
                return None
            # Get the header (e.g. "RS_GCF_000335255.2~NZ_AOEB01000154.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteoba...")
            header = line[1:].strip()
            # Split the header into GTDB ID and lineage part
            parts = header.split('~', 1)
            if len(parts) != 2:
                print(f"Warning: Skipping malformed header: {line.strip()}")
                return None
            # Get the GTDB ID (e.g. RS_GCF_905219285.2)
            gtdb_id = parts[0]
            # Extract the lineage_part (after the first space, and before the first ' [')
            # e.g. "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enter..."
            lineage_full_string = parts[1]
            lineage_part = lineage_full_string.split(' ', 1)[1].split(' [', 1)[0]
            return gtdb_id, lineage_part
        
        # Get all leaf node GTDB IDs from the decorated trees
        all_gtdb_leaf_ids = self.get_node_ids(id_type="gtdb", node_type="leaf")
        all_gtdb_leaf_ids = set(all_gtdb_leaf_ids)
        
        # Create the dictionary mapping GTDB IDs to taxonomy lists
        # e.g. "RS_GCF_018344175.1" -> ["d__Bacteria", "p__Pseudomonadota", "c__Gammaproteobacteria", "o__Enterobacterales", "f__Enterobacteriaceae", "g__Escherichia", "s__Escherichia coli"]
        gtdb_id_to_taxa = {}
        all_taxa_set = set()
        with open(seqs_fasta_path, 'r') as f:
            for line in f:
                parsed_result = parse_fasta_line(line)
                if parsed_result:
                    gtdb_id, lineage_part = parsed_result
                    # Check if the GTDB ID is in the set of all leaf node GTDB IDs
                    if gtdb_id in all_gtdb_leaf_ids:
                        taxa_list = lineage_part.split(';')
                        if gtdb_id in gtdb_id_to_taxa and taxa_list != gtdb_id_to_taxa[gtdb_id]:
                            print(f"Warning: Duplicate GTDB ID '{gtdb_id}' found in the FASTA file with different taxonomies.")
                            print(f"  First taxonomy: {gtdb_id_to_taxa[gtdb_id]}")
                            print(f"  Second taxonomy: {taxa_list}")
                        gtdb_id_to_taxa[gtdb_id] = taxa_list
                        all_taxa_set.update(taxa_list)
        if self.verbose:
            print(f"Successfully read the taxonomies of {len(gtdb_id_to_taxa)} sequences from {seqs_fasta_path}")
        
        # Create a dictionary mapping taxon names to lists of GTDB IDs
        # e.g. "g__Escherichia" -> ["RS_GCF_018344175.1", "RS_GCF_001246675.1", ...]
        taxon_to_gtdb_ids = defaultdict(list)

        if self.verbose:
            print("Grouping leaf nodes by taxon...")
        # Iterate through the parsed taxonomies
        for gtdb_id, taxa_list in tqdm(gtdb_id_to_taxa.items(), desc="Grouping taxa"):
            # For each taxon name in the lineage of the current GTDB ID
            for taxon_name in taxa_list:
                # Append the GTDB ID to the list for that taxon name
                taxon_to_gtdb_ids[taxon_name].append(gtdb_id)
        
        if self.verbose:
             print(f"Grouped leaf nodes into {len(taxon_to_gtdb_ids)} unique taxa.")

        # Now we have:
        # taxon_to_gtdb_ids: 
        # "g__Escherichia" -> ["RS_GCF_018344175.1", "RS_GCF_001246675.1", ...]
        # gtdb_id_to_taxa:
        # "RS_GCF_018344175.1" -> ["d__Bacteria", "p__Pseudomonadota", "c__Gammaproteobacteria", "o__Enterobacterales", "f__Enterobacteriaceae", "g__Escherichia", "s__Escherichia coli"]

        # Construct reduced_taxon_to_gtdb_ids, where the lists of IDs are subsets of the original lists in taxon_to_gtdb_ids
        # More specifically, for each taxon name, we identify the highest rank at which the sequences do not all have the same taxonomic label, then select a single sequence from each label
        reduced_taxon_to_gtdb_ids = {}

        # Construct reduced_taxon_to_gtdb_ids
        if self.verbose:
            print("Selecting representative leaf nodes for each taxon...")
        reduced_taxon_to_gtdb_ids = {}
        # Define rank order and mapping
        ranks = ['d', 'p', 'c', 'o', 'f', 'g', 's']
        rank_map = {rank: i for i, rank in enumerate(ranks)}

        # Iterate through each taxon and its associated GTDB IDs
        for taxon_name, gtdb_ids_list in tqdm(taxon_to_gtdb_ids.items(), desc="Reducing taxa lists"):
            # Handle taxa with less than 3 members - MRCA is trivial or undefined in this context
            if len(gtdb_ids_list) <= 2:
                reduced_taxon_to_gtdb_ids[taxon_name] = gtdb_ids_list # Keep original list (0, 1, or 2 elements)
                continue

            # Try to determine the rank of the current taxon_name
            current_rank_prefix = taxon_name[0]
            if current_rank_prefix not in rank_map:
                if self.verbose:
                    print(f"Warning: Could not determine rank for taxon '{taxon_name}'. Using first/last IDs as fallback.")
                reduced_taxon_to_gtdb_ids[taxon_name] = [gtdb_ids_list[0], gtdb_ids_list[-1]]
                continue
                
            current_rank_index = rank_map[current_rank_prefix]

            # Species level taxa have no lower ranks to check diversity against
            if current_rank_index == rank_map['s']:
                 reduced_taxon_to_gtdb_ids[taxon_name] = [gtdb_ids_list[0], gtdb_ids_list[-1]]
                 continue

            representatives_found = False
            # Check ranks below the current taxon's rank for diversity
            for check_rank_index in range(current_rank_index + 1, len(ranks)):
                groups_at_rank = defaultdict(list)
                unknown_rank_label = f"{ranks[check_rank_index]}__unknown_or_missing"

                # Group IDs based on their label at the check_rank_index
                for gtdb_id in gtdb_ids_list:
                    taxa_list = gtdb_id_to_taxa.get(gtdb_id) # Get the full lineage for this ID
                    if taxa_list and len(taxa_list) > check_rank_index:
                         # Check if the label at this rank index has the correct prefix
                        label_at_check_rank = taxa_list[check_rank_index]
                        if label_at_check_rank.startswith(ranks[check_rank_index] + '__'):
                             groups_at_rank[label_at_check_rank].append(gtdb_id)
                        else:
                             # Lineage might be inconsistent or malformed at this rank
                             groups_at_rank[unknown_rank_label].append(gtdb_id)
                    else:
                        # Lineage doesn't extend to this rank
                        groups_at_rank[unknown_rank_label].append(gtdb_id)
                
                # If we found more than one group at this rank, diversity exists
                if len(groups_at_rank) > 1:
                    # Select the first ID from each group as a representative
                    reduced_list = [ids[0] for ids in groups_at_rank.values() if ids] # Ensure list is not empty
                    # Need at least two representatives to find an MRCA
                    if len(reduced_list) >= 2:
                        reduced_taxon_to_gtdb_ids[taxon_name] = reduced_list
                        representatives_found = True
                        break # Stop checking lower ranks
                    else: 
                        # This case (finding diversity but ending up with < 2 representatives) 
                        # is unlikely but possible if groups had empty lists (shouldn't happen with current logic)
                        # or if only one group had valid members plus an 'unknown' group.
                        # Continue to potentially find better diversity lower down.
                        pass

            # Fallback: If no diversity was found at any lower rank
            if not representatives_found:
                # Use the first and last ID from the original list
                reduced_taxon_to_gtdb_ids[taxon_name] = [gtdb_ids_list[0], gtdb_ids_list[-1]]

        if self.verbose:
            print(f"Reduced lists of representative leaf nodes.")

        # Now we have:
        # reduced_taxon_to_gtdb_ids: 
        # "g__Escherichia" -> ["RS_GCF_018344175.1", "RS_GCF_001246675.1", ...]

        # Initialize the dictionary mapping taxon names to MRCA nodes
        self.node_from_taxon_name_dict = {}
        
        # Shuffle the order of reduced_taxon_to_gtdb_ids so that we can estimate time of completion better
        taxon_items = list(reduced_taxon_to_gtdb_ids.items())
        random.shuffle(taxon_items)

        # Iterate through shuffled reduced_taxon_to_gtdb_ids, find MRCA for each using the selected representatives,
        # and populate self.node_from_taxon_name_dict and the '.taxon_names' attribute on the MRCA nodes
        print("Finding MRCA nodes for each taxon...")
        nodes_with_taxa_assigned = 0
        failed_taxa_count = 0

        # Iterate through the shuffled taxa and their representative GTDB IDs
        for taxon_name, representative_gtdb_ids in tqdm(taxon_items, desc="Mapping taxa to MRCA nodes"):
            # Skip if there are fewer than 1 representative IDs, as MRCA is trivial or meaningless
            if len(representative_gtdb_ids) < 1:
                print(f"Warning: Skipping MRCA node for taxon '{taxon_name}': Needs at least 1 representative leaf nodes to find a meaningful MRCA.")
                failed_taxa_count += 1
                continue

            try:
                # Get the node objects for the representative GTDB IDs
                representative_nodes = [self.get_node(gtdb_id) for gtdb_id in representative_gtdb_ids]
                
                # Determine the correct tree (bacterial or archaeal)
                the_domain = representative_nodes[0].redvals_id[:3]
                tree = self.bac_tree if the_domain == "bac" else self.arc_tree

                # Find the MRCA of the representative nodes
                mrca_node = tree.common_ancestor(*representative_nodes) # Use * to unpack the list

                if mrca_node:
                    # Add the mapping from taxon name to the MRCA node
                    self.node_from_taxon_name_dict[taxon_name] = mrca_node
                    nodes_with_taxa_assigned += 1

                    # Add the taxon name to the MRCA node's taxon_names attribute (initialize if needed)
                    if not hasattr(mrca_node, 'taxon_names'):
                        mrca_node.taxon_names = set()
                    mrca_node.taxon_names.add(taxon_name)
                else:
                    # This case should be rare if input data is correct
                    print(f"Warning: Could not find MRCA for taxon '{taxon_name}' with representatives: {representative_gtdb_ids}. Skipping.")
                    failed_taxa_count += 1

            except KeyError as e:
                print(f"Warning: Skipping taxon '{taxon_name}'. Could not find node for GTDB ID: {e}. Ensure all IDs in the FASTA exist in the tree.")
                failed_taxa_count += 1
            except Exception as e:
                print(f"Warning: An unexpected error occurred while processing taxon '{taxon_name}': {e}. Skipping.")
                failed_taxa_count += 1

        print(f"\nFinished mapping taxa to nodes.")
        if self.verbose:
            print(f"  The 'node_from_taxon_name_dict' mapping is now populated.")
            total_taxa = len(reduced_taxon_to_gtdb_ids)
            assigned_count = len(self.node_from_taxon_name_dict)
            print(f"  Successfully mapped {assigned_count} taxa to MRCA nodes.")
            # Note: nodes_with_taxa_assigned might be higher than assigned_count if multiple taxa map to the same node.
            # print(f"  Assigned taxon names to {nodes_with_taxa_assigned} unique nodes.") # This count might be misleading if multiple taxa map to the same node. Use len(self.node_from_taxon_name_dict) instead.
            print(f"  Failed to map or skipped {failed_taxa_count} taxa (due to insufficient representatives, cross-domain issues, missing IDs, or MRCA errors).")
            print(f"  Total unique taxa considered: {total_taxa}.")
            print(f"  Total unique taxa in the input FASTA: {len(all_taxa_set)}.")
        
        if save_result_path:
            self.save_taxa_to_node_mapping(save_result_path)

    def save_taxa_to_node_mapping(self, save_result_path):
        """
        Save the taxon names to MRCA nodes mapping to a pickle file.
        save_result_path (str): Path to save the mapping to. (e.g. "./taxon_mapping/taxon_to_node_mapping.pkl")
        """
        if self.node_from_taxon_name_dict is None:
            raise ValueError("node_from_taxon_name_dict has not been created. Populate it first.")
        is_valid_type = save_result_path.endswith('.pkl')
        is_valid_dir = os.path.isdir(os.path.dirname(save_result_path))
        if not is_valid_type or not is_valid_dir:
            raise ValueError(f"Invalid save_result_path: {save_result_path}. Call save_taxa_to_node_mapping with a valid path.")
        # First, convert mappings from taxon names -> nodes to taxon names -> redvals IDs
        taxon_name_to_redvals_id_dict = {taxon_name: node.redvals_id for taxon_name, node in self.node_from_taxon_name_dict.items()}
        # Then, save the dictionary to a pickle file
        with open(save_result_path, 'wb') as f:
            pickle.dump(taxon_name_to_redvals_id_dict, f)

    def load_taxa_to_node_mapping(self, load_result_path):
        """
        Load the taxon names to MRCA nodes mapping from a pickle file.
        load_result_path (str): Path to load the mapping from. (e.g. "./taxon_mapping/taxon_to_node_mapping.pkl")
        """
        if not os.path.exists(load_result_path):
            raise FileNotFoundError(f"File not found: {load_result_path}")
        with open(load_result_path, 'rb') as f:
            redvals_id_from_taxon_name_dict = pickle.load(f)
        # Now use redvals_id_from_taxon_name_dict to create the node_from_taxon_name_dict
        self.node_from_taxon_name_dict = {taxon_name: self.get_node(redvals_id) for taxon_name, redvals_id in redvals_id_from_taxon_name_dict.items()}

    def get_distance_in_taxon(self, taxon_name):
        """
        Given a taxon name, return the RED distance between any pair of leaf nodes (sequences) that have the node representing that taxon as their MRCA.

        taxon_name (str): The name of the taxon (e.g., "g__Escherichia" or "d__Bacteria" or "s__Spirillospora terrae").
        """
        # Check if the trees are decorated
        if not self.is_decorated():
            raise ValueError("The trees must be decorated to get distance from taxon")
        # Check if the node_from_taxon_name_dict exists
        if self.node_from_taxon_name_dict is None:
            raise ValueError("node_from_taxon_name_dict has not been created. Run a method to populate it first.")
        # Look up the node
        if taxon_name not in self.node_from_taxon_name_dict:
            raise KeyError(f"Taxon name '{taxon_name}' not found in node_from_taxon_name_dict")
        
        node = self.node_from_taxon_name_dict[taxon_name]
        
        # Check if the node has the red_distance attribute
        if not hasattr(node, 'red_distance'):
             # This case should ideally not happen if the tree is decorated, but check for robustness
            raise AttributeError(f"Node for taxon '{taxon_name}' does not have a red_distance attribute.")
            
        return node.red_distance
    
    def get_node_from_taxon_name(self, taxon_name):
        """
        Given a taxon name, return the node representing that taxon.
        """
        return self.node_from_taxon_name_dict[taxon_name]
    


class NodeInfo:
    """Container for information about a node in a decorated GTDB phylogenetic tree.
    
    Attributes:
        domain (str): Domain of the node ("bac" or "arc")
        gtdb_id (str): GTDB identifier for the node (e.g., "GB_GCA_002687935.1")
        redvals_id (str): Internal ID assigned by redvals (e.g., "bac00001342")
        red_value (float): RED value for the node
        is_terminal (bool): Whether the node is a terminal/leaf node (True) or an internal node (False)
        red_distance (float): RED distance between two leaf nodes which share this node as their MRCA, calculated as 2 * (1 - red_value)
    """
    
    def __init__(self, domain, red_value, gtdb_id, redvals_id, is_terminal, red_distance=None):
        """Initialize a NodeInfo object with information about a tree node."""
        self.domain = domain
        self.gtdb_id = gtdb_id
        self.redvals_id = redvals_id
        self.red_value = red_value
        self.is_terminal = is_terminal
        self.red_distance = None if red_value is None else (1 - red_value) * 2 if red_distance is None else red_distance

    def __str__(self):
        """Return a human-readable string representation of the NodeInfo object."""
        return (f"NodeInfo:\n"
                f"  redvals ID: {self.redvals_id}\n"
                f"  GTDB ID: {self.gtdb_id}\n"
                f"  Domain: {self.domain}\n"
                f"  RED value: {self.red_value}\n"
                f"  RED distance: {self.red_distance}\n"
                f"  Is terminal node: {self.is_terminal}")
    
    def __repr__(self):
        """Return a string representation that can be used to recreate the NodeInfo object."""
        return f"NodeInfo(domain='{self.domain}', gtdb_id='{self.gtdb_id}', redvals_id='{self.redvals_id}', red_value={self.red_value}, is_terminal={self.is_terminal}, red_distance={self.red_distance})"

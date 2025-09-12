"""
Example 1 - Decorating Trees from Scratch

Copyright (c) 2025 Haig Bishop

MIT License (see LICENSE file for details)

This example demonstrates decorating a tree with RED values:
1. Load the original Newick tree files (.tree) into a RedTree object
2. Decorate the trees with calculated RED values
3. Write the decorated trees to files (.pkl)

"""

# Import the RedTree class
from redvals import RedTree

import os


# Original Newick Tree Files -------------
# GTDB Release -------------
GTDB_RELEASE = "r226"

# The original Newick format archeal GTDB phylogenetic tree file
ARC_TREE_PATH = f"trees/ar53_{GTDB_RELEASE}.tree"
# The original Newick format bacterial GTDB phylogenetic tree file
BAC_TREE_PATH = f"trees/bac120_{GTDB_RELEASE}.tree"

# Decorated Tree Files (OUTPUT) -------------
# A Bio.Phylo.Newick.Tree object holding the decorated archeal GTDB phylogenetic tree
ARC_DECORATED_TREE_PATH = f"decorated_trees/ar53_{GTDB_RELEASE}_decorated.pkl"
# A Bio.Phylo.Newick.Tree object holding the decorated bacterial GTDB phylogenetic tree
BAC_DECORATED_TREE_PATH = f"decorated_trees/bac120_{GTDB_RELEASE}_decorated.pkl"


# 1. Initialise RedTree object -----------
# We are using the original undecorated trees (.tree files) as input
red_trees = RedTree(BAC_TREE_PATH, ARC_TREE_PATH)


# 2. Decorate the trees -----------
# The trees will not be decorated with RED values because we loaded them from .tree files
print("Are the trees are decorated:", red_trees.is_decorated())
# Decorate the trees with RED values
red_trees.decorate_from_calc()

# Get some random nodes and print their RED values
node = red_trees.get_node("bac00090342")
print("The RED value of bac00090342 is", node.red_value)

node = red_trees.get_node("bac00000002")
print("The RED value of bac00000002 is", node.red_value)

node = red_trees.get_node("arc00002342")
print("The RED value of arc00002342 is", node.red_value)

node = red_trees.get_node("arc00000147")
print("The RED value of arc00000147 is", node.red_value)


# 3. Write the decorated trees to files -----------
# The trees are now decorated with RED values
print("Are the trees are decorated:", red_trees.is_decorated())

# If the tree is decorated (with RED values)
if red_trees.is_decorated():
    # If the decorated trees don't already exist
    if not os.path.exists(BAC_DECORATED_TREE_PATH) and not os.path.exists(ARC_DECORATED_TREE_PATH):
        # Write the decorated trees to file
        red_trees.write_decorated_trees(BAC_DECORATED_TREE_PATH, ARC_DECORATED_TREE_PATH)




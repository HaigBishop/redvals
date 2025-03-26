"""
Example usage of the output of main.py
"""


# Load Dictionary ----------------------------
# First, load X
# TODO


# Use Dictionary ----------------------------
# Print the RED value for 
# TODO


from Bio import Phylo
# Load the Newick tree file of bacterial GTDB tree
bac120_tree = Phylo.read("./trees/bac120_r220.tree", "newick")
# Print some info about the tree
terminal_nodes = bac120_tree.get_terminals()
nonterminal_nodes = bac120_tree.get_nonterminals()
print("10 leaf node names:", [node.name for node in terminal_nodes[:10]])
print("10 internal node names:", [node.name for node in nonterminal_nodes[:10]])
print("Number of leaf nodes:", len(terminal_nodes))
print("Number of internal nodes:", len(nonterminal_nodes))


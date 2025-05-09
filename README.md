# redvals
Tool for obtaining and accessing Relative Evolutionary Divergence (RED) values from GTDB phylogenetic trees.


## Set Up
1. Clone repository
```
git clone https://github.com/HaigBishop/redvals.git
cd redvals
```
2. Install dependencies (biopython, pandas and tqdm)
```
conda create -n redvals_env python=3.12 biopython pandas tqdm
conda activate redvals_env
```
3. Use package
 - See usage below
 - Or see in depth examples in `example_1.py` and `example_2.py`


## Usage
### Import redvals.RedTree
The `RedTree` object represents the bacterial and archaeal phylogenetic trees.
```python
from redvals import RedTree
```

### Load Undecorated Trees (Init Method 1)
This loads original Newick trees with no RED values. Also assigns 'redvals IDs' to all nodes.
```python
red_trees = RedTree("trees/bac120_r220.tree", "trees/ar53_r220.tree")
```

### Decorate Trees from TSVs (takes 30-90 minutes)
This adds RED values and RED distances to every node in both trees.
```python
red_trees.decorate_from_tsv("red_values/gtdbtk_r220_bac120.tsv", "red_values/gtdbtk_r220_ar53.tsv")
```

### Write Decorated Trees
This saves the trees with the RED values and distances as .pkl files, avoiding repetition of the time-consuming decoration process.
```python
red_trees.write_decorated_trees("decorated_trees/bac120_r220_decorated.pkl", "decorated_trees/ar53_r220_decorated.pkl")
```

### Load Decorated Trees (Init Method 2)
Load the already decorated trees originating from write_decorated_trees.
```python
red_trees = RedTree("decorated_trees/bac120_r220_decorated.pkl", "decorated_trees/ar53_r220_decorated.pkl")
```

### Convert Node IDs
The IDs found in the orginal GTDB `.tree` files are maintained (called `gtdb_ids`). But for most internal nodes there are no IDs. Therefore new IDs are assigned to *all* nodes called `redvals_ids`. You can use either type at any time, but if you wish to convert between them you can use the code here.
```python
bacterial_redvals_id = red_trees.get_redvals_id("GB_GCA_002687935.1")
gtdb_id = red_trees.get_gtdb_id("bac00000001")
```

### Get All Node IDs
You can retrieve all node IDs with optional filtering.
```python
all_gtdb_ids = red_trees.get_node_ids(domain='both', node_type='both', id_type='gtdb') # (many internal nodes do not have GTDB IDs)
internal_archaeal_redvals_ids = red_trees.get_node_ids(domain='arc', node_type='internal', id_type='redvals')
```

### Get All Nodes
You can retrieve all nodes with optional filtering.
```python
all_leaf_nodes = red_trees.get_nodes(domain='both', node_type='leaf')
all_bacterial_nodes = red_trees.get_nodes(domain='bac', node_type='both')
all_bacterial_leaf_nodes = red_trees.get_nodes(domain='bac', node_type='leaf')
```

### Access Node Info
Many attributes of nodes can easily be accessed using get_node_info.
```python
node_info = red_trees.get_node_info("bac00000001")
print(f"GTDB ID: {node_info.gtdb_id}, redvals ID: {node_info.redvals_id}, Domain: {node_info.domain}")
print(f"RED value: {node_info.red_value}, RED distance: {node_info.red_distance}, Is terminal?: {node_info.is_terminal}")
```

### Compute RED Distances
For any two nodes dist_between_nodes gives the RED distance between them and the `redvals_ids` of their MRCA. This works for two leaf nodes, two internal nodes, or an internal and a leaf node.
```python
red_distance, mrca_node_id = red_trees.dist_between_nodes("bac00000001", "RS_GCF_001186155.3")
```

### Map Taxon Names to Nodes
Given a taxon name (e.g. g__Escherichia) we can use get_distance_in_taxon(taxon_name) to get the RED distance between any pair of leaf nodes who's MRCA is the node representing that taxon. Before using get_distance_in_taxon(taxon_name) however, map_nodes_to_taxa(seqs_fasta) must be run, where seqs_fasta is a GTDB database FASTA file such as ssu_all_r220.fna.
```python
# Running map_taxa_to_nodes may take 20-40 minutes
red_trees.map_taxa_to_nodes("D:/16S_databases/ssu_all_r220.fna", 
    save_result_path="./taxon_mapping/taxon_to_node_mapping.pkl")
# Then these are fast
red_distance = red_trees.get_distance_in_taxon("p__Nitrospirota")
s_terrae_node = red_trees.get_node_from_taxon_name("s__Spirillospora terrae")
```

### Load Pre-Computed Taxon Name Mappings
If you previously ran map_taxa_to_nodes and saved the result, you can use load_taxa_to_node_mapping to quicky load taxon -> node mappings
```python
red_trees.load_taxa_to_node_mapping("./taxon_mapping/taxon_to_node_mapping.pkl")
red_distance = red_trees.get_distance_in_taxon("p__Nitrospirota")
s_terrae_node = red_trees.get_node_from_taxon_name("s__Spirillospora terrae")
```


## Files

### Original Newick Tree Files
 - `./trees/bac120_r220.tree`
 - `./trees/ar53_r220.tree`

These are the GTDB phylogenetic trees (release 220) in Newick format. Obtained from: https://gtdb.ecogenomic.org/downloads

Example Usage:
```python
from Bio import Phylo

# Load the Newick tree file of bacterial GTDB tree
bac120_tree = Phylo.read("./trees/bac120_r220.tree", "newick")

# Print some info about the tree
terminal_nodes = bac120_tree.get_terminals()
nonterminal_nodes = arc53_tree.get_nonterminals()
print("5 leaf node names:", [node.name for node in terminal_nodes[:5]])
print("5 internal node names:", [node.name for node in nonterminal_nodes[:5]])
print("Number of leaf nodes:", len(terminal_nodes))
print("Number of internal nodes:", len(nonterminal_nodes))
```

### RED Values TSV Files
 - `./red_values/gtdbtk_r220_bac120.tsv`
 - `./red_values/gtdbtk_r220_ar53.tsv`

These TSV files contain the RED values for all nodes in each tree. They originate from the release 220 gtdbtk_package directory from: https://gtdb.ecogenomic.org/downloads
i.e. https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/

There is one row for each node (terminal and nonterminal). The second column holds the RED value assigned to the given node. The first column holds one or two leaf IDs (AKA genome IDs). For example in the first column we might find "GB_GCA_020721905.1" or "GB_GCA_026414805.1|GB_GCA_020721905.1". When there is a single leaf ID, this is a leaf node, and therefore the RED value in column two is 1.0. But if there are two leaf IDs seperated by a "|" symbol, this row represents an internal node. Specifically, it represents the node that is the most recent common ancestor (MRCA) of the two leaf IDs in column one.

### Decorated Tree Files
 - `./decorated_trees/bac120_r220_decorated.pkl`
 - `./decorated_trees/ar53_r220_decorated.pkl`

These files are Python pickle files, each containing a Bio.Phylo.Newick.Tree object. Loading these objects (see below) results in the same object as loading the .tree files (as above), but the difference is that they are decorated with RED values, and RED distances.

Example Usage:
```python
from Bio import Phylo
import pickle

# Load the pickle file of decorated bacterial GTDB tree
with open("./out/bac120_r220_decorated.pkl", "rb") as f:
    decorated_bac120_tree = pickle.load(f)

# Access a RED value using the tree
node_520 = decorated_bac120_tree.get_node("bac00000520")
print("RED value of 'bac00000520':", node_520.red_value)
```

### Taxon to Node Mapping
 - `./taxon_mapping/taxon_to_node_mapping.pkl`

This files is a Python pickle file containing a dictionary which simply maps `taxon_name` to `redvals_id` where `redvals_id` is the ID of the node representing that clade. For example, "g__Escherichia" -> "bac00002631" or "p__Nitrospirota" -> "bac00079003". This is used for get_node_from_taxon_name() and get_distance_in_taxon().




## How to Contribute
Contributions are very welcome.

1. **Fork the repository** and clone it locally
2. **Create a new branch** for your feature or bug fix
3. **Make your changes** and write tests if applicable
4. **Run existing tests** to ensure nothing was broken
5. **Submit a pull request** with a clear description of your changes


## To Do
 - add good docstrings to RedTree
 - add comprehensive tests to all_tests.py
 - add to PyPI


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Data Licensing

This tool processes data from the Genome Taxonomy Database (GTDB). The GTDB data is subject to its own licensing terms:
- GTDB data (including phylogenetic trees and taxonomy) is available under a [CC BY-SA 4.0 license](https://creativecommons.org/licenses/by-sa/4.0/)


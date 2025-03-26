# redvals
Tool for obtaining RED values from GTDB phylogenetic trees


## Usage
1. Install dependencies (biopython and pandas)
```
conda install -c conda-forge biopython
conda install -c conda-forge pandas
```
2. Unless the `out/` directory is already populated, run `main.py`
3. Make use of the output files (as in `example.py`)


## Input Files
####trees/{DOMAIN}_r220.tree
These are the GTDB phylogenetic trees in Newick format.
- where from
- example loading 

####red_values/gtdbtk_r220_{DOMAIN}.tsv
These TSV files contain the RED values for all internal nodes for each tree.
- explain format


## Output Files
####out/{DOMAIN}_r220_decorated.pkl
These files are Python pickle files, each containing a Bio.Phylo.Newick.Tree object.
- decoration (compare to .tree and say they combine the .tsv and .tree files)
- example use

####out/dict.pkl
TODO


![Phylogenetic tree visualization](res/tree.png)

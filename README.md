# PhyloCharMod (publication data)
Scripts and data for publication

## Phylogenetic prediction of functional sequence modules

An original approach to characterize functional motifs based on:
1. Detection of conserved sequence modules (using Partial Local Multiple Alignment)
2. Phylogenetic inference of species/genes/modules/functions evolutionary histories
3. Identification of co-appearances of modules and functions

![Data recuperation](img/method.png)

## Usage

### To replicate publication results :

```cd data/min5_human_214_t10m1M20/```

```python3 ../../phylocharmod/integrate_3phylo.py seadogMD_214.output gene_tree_214/214.tree --pastml_tab acs_dir_seadogMD_214_gene/pastml_seadogMD_214_gene_leaf_Manual_214_combined_ancestral_states.tab --domains_csv domains_214.csv --itol``` 

### To use on you own data :

```python3 phylocharmod/phylocharmod.py```

(see ```--help``` for details)

```usage: phylocharmod.py [-h] [--output_directory OUTPUT_DIRECTORY] [--species_tree SPECIES_TREE] [--gene_tree GENE_TREE] [--plma_file PLMA_FILE] [--reconc_domains] multi_fasta_file leaf_functions_csv

positional arguments:
  multi_fasta_file      Multi fasta file, with specific formated header >RefSeq_taxid (ex : >XP_012810820.2_8364)
  leaf_functions_csv    csv file containing for each of our sequence, the list of his functions (ex : XP_012810820.2, P59509 | P999999)

optional arguments:
  -h, --help            show this help message and exit
  --output_directory OUTPUT_DIRECTORY
                        output directory name
  --species_tree SPECIES_TREE
                        Species tree to use as a support for the reconciliations (WARNING, must correspond to the taxid use in the other files !)
  --gene_tree GENE_TREE
                        Gene tree to use as a support for the pastML and DGS reconciliation inference (WARNING, must correspond to the sequences in the multi fasta file !)
  --plma_file PLMA_FILE
                        Paloma-2 output file (.agraph format, .dot, or .oplma format)
  --reconc_domains      Do a DGS reconciliation with known modules (pfam / prosite)
```

## Dependencies

All these programs are necessary to run *PhyloCharMod*, and must be in ;
```/usr/local/bin/```
If not, their path must be specified in the config file ;
```phylocharmod/config.txt```

#### Python3 packages

[ETE Toolkit](http://etetoolkit.org/), a python3 framework for the analysis and visualisation of trees.

[Six]()
```pip install six```

[Bio](https://github.com/biopython/biopython), the Biopython Project is an international association of developers of freely available Python tools for computational molecular biology.

#### Softwares

[Muscle](http://www.drive5.com/muscle/), one of the best-performing multiple alignment programs, [Conda package](https://anaconda.org/bioconda/muscle)

[PhyML](https://github.com/stephaneguindon/phyml), maximum likelihood phylogenetic inference for the gene and the module trees, [Conda package](https://anaconda.org/bioconda/phyml)

[TreeFix](https://www.cs.hmc.edu/~yjw/software/treefix/), Statistically Informed Gene Tree Error Correction Using Species Trees, [Conda package](https://anaconda.org/OcMalde/treefix)

[SEADOG-MD](https://compbio.engr.uconn.edu/software/seadog/), for DGS-reconciliation.

[PastML](https://pastml.pasteur.fr/), for ancestral characters inference, [Pip package](https://pypi.org/project/pastml/)




# PhyloCharMod (publication data)
Scripts and data for publication

## Phylogenetic prediction of functional sequence modules

An original approach to characterize functional motifs based on:
1. Detection of conserved sequenc modules (using Partial Local Multiple Alignment)
2. Phylogenetic inference of species/genes/modules/functions evolutionary histories
3. Identification of co-appearances of modules and functions

![Data recuperation](img/method.png)

## Usage

### To replicate publication results :

```cd data/min5_human_214_t10m1M20/```

```python3 ../../phylocharmod/integrate_3phylo.py seadogMD_214.output gene_tree_214/214.tree --pastml_tab acs_dir_seadogMD_214_gene/pastml_seadogMD_214_gene_leaf_Manual_214_combined_ancestral_states.tab --domains_csv domains_214.csv --itol``` 

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

[Bioservices](https://pypi.org/project/bioservices/), (only for annotations recuperation, if annotations are provided an other way, its not necessary.

#### Muscle

[Muscle](http://www.drive5.com/muscle/), one of the best-performing multiple alignment programs.

[Conda package](https://anaconda.org/bioconda/muscle)

#### Trimal

[Trimal](http://trimal.cgenomics.org/), a tool for the automated removal of spurious sequences or poorly aligned regions from a multiple sequence alignment. 

[Conda package](https://anaconda.org/bioconda/trimal)

#### PRANK (no required at the moment)

[PRANK](http://wasabiapp.org/software/prank/), a probabilistic multiple alignment program used to align the gene sequences.

[Conda package](https://anaconda.org/bioconda/prank)

#### PhyML

[PhyML](https://github.com/stephaneguindon/phyml), maximum likelihood phylogenetic inference for the gene and the module trees.

[Conda package](https://anaconda.org/bioconda/phyml)

#### TreeFix

[TreeFix](https://www.cs.hmc.edu/~yjw/software/treefix/), Statistically Informed Gene Tree Error Correction Using Species Trees.

[Conda package](https://anaconda.org/OcMalde/treefix)

#### Paloma

[Paloma](http://tools.genouest.org/tools/protomata/learn/), Partial Local Multiple alignment, for the module segmentation of our sequences.

???

#### SEADOG-MD

[SEADOG-MD](https://compbio.engr.uconn.edu/software/seadog/), for DGS-reconciliation.

???

#### PastML

[PastML](https://pastml.pasteur.fr/), for ancestral characters inference.

[Pip package](https://pypi.org/project/pastml/)

## Usage

```python3 phylocharmod myfasta.fasta gene_functions.csv```

How it work, figure, schemas, results




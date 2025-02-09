# BetaMining
A Python 3 package to mine the AlphaFold2 monomer database for secondary structure predictions. Written for the [Halfmann Lab](https://research.stowers.org/halfmannlab/) at the [Stowers Institute for Medical Research](https://www.stowers.org/). This project was started by Paula Berry, and resumed by myself.

## Installation
This package can be installed via `pip install git+https://github.com/jacobnjensen1/BetaMining.git`.

### Dependencies
This package requires Python 3.5 or greater. The main packages used to perform calculations are:

```
biopython>=1.79
biopandas>=0.4.1
numpy>=1.21
pandas>=1.3
ProDy>=2.0.2
scikit-learn>=1.0.2
```

`DSSP` is optional based on the desired analysis.

## Usage
When installed, `BetaMining` can be run from the command line through the wrapper `run_beta_mining.py`.

`run_beta_mining.py -d` or `run_beta_mining.py --defaults` generates a default config YAML file in the current working directory, `default_config.yaml`. It contains the options for configuring a `BetaMining` session.

`run_beta_mining.py -c <CONFIG>` or `run_beta_mining.py --config <CONFIG>` executes the script using the parameters and filepaths indicated in the YAML formatted config file provided by the user.

`BetaMining` is intended to be used on `.pdb` files generated by AlphaFold, or downloaded from the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/). The bulk download `.tar` archives must be unpacked before use, but `BetaMining` can handle both `.pdb` and `.pdb.gz` files natively.

#### `structure_dictionary.json`

Structures are defined in `structure_dictionary.json`, which contains some useful defaults. If you want to change the defaults, make a copy of the original and specify the path to the new file in `<config>.yaml`. The useful characteristics of a target region are:
  * `"name"`: the name of the structure
  * `"regex"`: a regular expression used to identify that structure in a string composed of `"mask_symbol"`s
  * `"regex_flank"`: the number of residues on each side of the structure to include in the output files
  * `"include"`: characteristics that a potentially interesting region needs to have. Possibilities are:
    * `"contacts"`: defines a class of intramolecular contacts.
    * `"twist"`: requires the twist of residues in the region to fall between two bounds, options are `"mean": [min,max]`, `"min": [min,max]`, and `"max": [min,max]`
    * `"plddt"`: requires the confidence scores of residues in the region to fall between two bounds, the options are the same as `"twist"`
    * `"flank_plddt"`: calculates confidence of the `"size"` flanking residues of the region. `"mean_ratio"` requires the average flank confidence / average in region confidence to fall between two bounds. `"mean"` requires the average flank confidence to fall between two bounds.
  * `"exclude"`: characteristics that a potentially interesting region cannot have. This section has the same options as `"include"` with the exception of `"flank_plddt"`.
    * `"out_of_plane_sheets"`: defines restrictions on the orientation of adjacent beta sheets. This is optional, but highly recommended, and requires `DSSP`.

#### `config.yaml`

This file specifies the locations of the input and output and some options for the analysis.
Detailed explanations can be found in the file.

### Outputs
The default behavior of `BetaMining` is to produce a `.csv` containing the sequences it identified as a hit, as well as metadata pulled from the `.pdb` file. It will also produce these sequences as a multiple-sequence `.fasta` file, with limited metadata in the header, also pulled from the `.pdb` file.

Complete calculations and secondary structure assignments can be generated for every `.pdb` file analyzed. This is to aid in the identification of characteristics in the protein structure predictions that might explain why or why not `BetaMining` is identifying them as hits.
This is performed when `"per_residue_output"` is `True` in the config file, but should not be performed when running `BetaMining` on a large number of models.

## Concept Overview
The `BetaMining` algorithm looks for occult (hidden) secondary structure patterns in AlphaFold monomer protein structure predictions. Many proteins that are involved in amyloid formation have aggregation regions that are unstructured in monomer form, and are therefore classified as "low confidence" regions in the AlphaFold prediction framework. However, the predicted dihedral angles of these regions are often within the Ramachandran plot area where they appear in complex.

### beta-run
To classify an amino acid sequence as a `beta-run`, the sequence must:

* contain two or more beta-sheet assigned residues in a row
* contain three or more of these "beta-sheet" runs
* separated by no more than 6 non-beta-sheet assigned residues

Within the `beta-run` each residue must:

* not be within 6 angstroms of another residue more than 3 residues away
* have another beta-sheet assigned residue within 14 angstroms, but not less than 6 angstroms away, and more than 2 residues away

### sextuple-helix
An experimental mode for identifying characteristic death domain protein structures. To classify an amino acid sequence as a `sextuple-helix`, the sequence must:

* be between 80 and 120 residues long
* contain 6 instances of 4 or more alpha-helix assigned residues in a row
* contain multiple alpha-helix assigned residues within 10 angstroms of another alpha-helix assigned residue that is more than 4 residues away

## Future directions

Do you have a secondary structure you think might be hiding in monomeric predictions? Let me know!

### To do
* detailed logging to indicate why different proteins were rejected or not
* create example notebooks


## References

### Modules

#### ProDy
Project Link: http://prody.csb.pitt.edu/index.html

* PDB parsing
* dihedral angle calculation

Zhang S, Krieger JM, Zhang Y, Kaya C, Kaynak B, Mikulska-Ruminska K, Doruker P, Li H, Bahar I [ProDy 2.0: Increased scale and scope after 10 years of protein dynamics modelling with Python](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab187/6211036) **2021** *Bioinformatics*, btab187.

Bakan A, Meireles LM, Bahar I [ProDy: Protein Dynamics Inferred from Theory and Experiments](https://academic.oup.com/bioinformatics/article/27/11/1575/217006?login=true) **2011** *Bioinformatics* 27(11):1575-1577

Bakan A, Dutta A, Mao W, Liu Y, Chennubhotla C, Lezon TR, Bahar I [Evol and ProDy for Bridging Protein Sequence Evolution and Structural Dynamics](http://bioinformatics.oxfordjournals.org/content/30/18/2681.long) **2014** *Bioinformatics* 30(18):2681-2683

#### BioPython
Project Link: https://biopython.org/

* PDB header parsing

Cock PA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B and de Hoon MJL (2009) [Biopython: freely available Python tools for computational molecular biology and bioinformatics.](https://academic.oup.com/bioinformatics/article/25/11/1422/330687) *Bioinformatics*, 25, 1422-1423

Hamelryck T and Manderick B (2003) [PDB file parser and structure class implemented in Python.](https://academic.oup.com/bioinformatics/article/19/17/2308/205793) *Bioinformatics*, 22, 2308-2310


#### BioPandas
Project Link: http://rasbt.github.io/biopandas/

* PDB header parsing
* PDB conversion to Pandas Dataframe
* 3-letter amino acid codes to single letter amino acid codes conversion

Sebastian Raschka. Biopandas: Working with molecular structures in pandas dataframes. *The Journal of Open Source Software*, 2(14), jun 2017. doi: 10.21105/joss.00279. URL http://dx.doi.org/10.21105/joss.00279.

#### NumPy
Project Link: https://numpy.org/

* triangular matrix calculations
* matrix addition
* vector operations

Harris, C.R., Millman, K.J., van der Walt, S.J. et al. [Array programming with NumPy.](https://www.nature.com/articles/s41586-020-2649-2) *Nature* 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.

#### Sci-kit Learn
Project Link: https://scikit-learn.org/stable/

* pair-wise distance calculation

[Scikit-learn: Machine Learning in Python](https://jmlr.csail.mit.edu/papers/v12/pedregosa11a.html), Pedregosa *et al.,* JMLR 12, pp. 2825-2830, 2011.

#### pandas
Project Link: https://pandas.pydata.org/docs/index.html

* dataframe functions
* data import and export functions

[Data structures for statistical computing in python](https://conference.scipy.org/proceedings/scipy2010/pdfs/mckinney.pdf), McKinney, Proceedings of the 9th Python in Science Conference, Volume 445, 2010.

### Protein Structure

#### protein structure prediction
AlphaFold Database Link: https://alphafold.ebi.ac.uk/

Jumper, J *et al.* [Highly accurate protein structure prediction with AlphaFold.](https://www.nature.com/articles/s41586-021-03819-2) Nature (2021).

Varadi, M *et al.* [AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models.](https://academic.oup.com/nar/article/50/D1/D439/6430488) Nucleic Acids Research (2021).

#### DSSP
Touw, W.G. *et al.* [A series of PDB-related databanks for everyday needs.](https://academic.oup.com/nar/article/43/D1/D364/2435537) Nucleic Acids Research (2014).

Kabsch, Wolfgang and Christian Sander [Dictionary of protein secondary structure: pattern recognition of hydrogen‐bonded and geometrical features](https://onlinelibrary.wiley.com/doi/abs/10.1002/bip.360221211) Biopolymers (1983).

#### secondary structure assignment from dihedral angles
Daniel Ting, Guoli Wang, Maxim Shapovalov, Rajib Mitra, Michael I. Jordan, Roland L. Dunbrack, Jr. [Neighbor-dependent Ramachandran probability distributions of amino acids developed from a hierarchical Dirichlet process model.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000763) *PLOS Comp. Biol.* (April 2010).

#### backbone twist calculations and values

Note: twist is no longer a default parameter in `structure_dictionary.json`

Cyrus Chothia, Conformation of twisted β-pleated sheets in proteins, *Journal of Molecular Biology,* Volume 75, Issue 2, **1973**, Pages 295-302, ISSN 0022-2836,
https://doi.org/10.1016/0022-2836(73)90022-3.

Shamovsky, I.L., Ross, G.M., Riopelle, R.J. [Theoretical Studies on the Origin of β-sheet Twisting](https://pubs.acs.org/doi/pdf/10.1021/jp002590t). *J. Phys. Chem. B*, **2000**. 104, 47, 11296-11307. https://doi.org/10/1021/jp002590t

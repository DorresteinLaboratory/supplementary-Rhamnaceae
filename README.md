# supplementary-Rhamnaceae
This repository contains supplementary materials relating to the manuscript "Comprehensive mass spectrometry-guided plant specialized metabolite phenotyping reveals metabolic diversity in the cosmopolitan plant family Rhamnaceae".

# Citation

Kang KB#, Ernst M#, van der Hooft JJJ#, da Silva RR, Park J, Medema M, Sung SH, Dorrestein PC. Comprehensive mass spectrometry-guided phenotyping of plant specialized metabolite reveals metabolic diversity in the cosmopolitan plant family Rhamnaceae. *The Plant Journal* (in press) (2019) (#Denotes equal contribution) [DOI:10.1111/tpj.14292](https://doi.org/10.1111/tpj.14292)

# Jupyter notebooks

### MZmine_GroupMapping.ipynb

Jupyter notebook used to map metadata onto the mass spectral molecular network in Cytoscape version 3.4.0 ([Shannon et al., 2003](https://genome.cshlp.org/content/13/11/2498.full)) for the MZmine 2 ([Pluskal et al., 2010](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-395)) preprocessed dataset.

# Folders

## ClassyFire

The folder called **ClassyFire** contains R scripts and data files used for performing automated chemical classification of the <i>in silico</i> annotated structures using ClassyFire ([Djoumbou Feunang et al., 2016](https://jcheminf.springeropen.com/articles/10.1186/s13321-016-0174-y)).

## Mass2Motifs_2_MolecularNetwork

The folder called **Mass2Motifs_2_MolecularNetwork** contains R scripts and data files used for mapping Mass2Motifs ([Van der Hooft et al., 2016](http://www.pnas.org/content/113/48/13738.full); [Wandy et al., 2018](https://academic.oup.com/bioinformatics/article/34/2/317/4158166)) on the mass spectral molecular networks ([Wang et al., 2016](https://www.nature.com/articles/nbt.3597); [Watrous et al., 2012](http://www.pnas.org/content/109/26/E1743)). 

## Heatmap_Subclasses

The folder called **Heatmap_subclasses** contains Jupyter notebook and data files for chemical subclass distribution heatmap (Figure 3(c)).

## Heatmap_Mass2Motifs

The folder called **Heatmap_Mass2Motifs** contains Jupyter notebook and data files for Mass2Motifs distribution heatmap (Figure 4).
 
### Motifs_featuretable.ipynb

Jupyter notebook used for creating Mass2Motifs-related feature table with ion count cutoff 1000.

### Heatmap_Motif_Drawing.ipynb

Jupyter notebook used for drawing the heatmap.

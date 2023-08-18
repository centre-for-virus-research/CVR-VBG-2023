---
layout: default
title: Phylogenetics
rank: 4
---
 
 # Phylogenetics Practical

David L Robertson, MRC-University of Glasgow Centre for Virus Research

[**david.l.robertson@glasgow.ac.uk**](mailto:david.l.robertson@glasgow.ac.uk)


**Aim**

To introduce multiple sequence alignment and the inference of evolutionary history. You will learn how to align homologous virus sequence data and construct a phylogenetic tree, use some different methods and how to test the reliability of clades in your phylogeny. 


**Task**

To generate the phylogenetic tree in Figure 1A of the paper [Iyer et al. 2017, Resistance to type 1 interferons is a major determinant of HIV-1 transmission fitness. PNAS 114(4):E590-E599](https://www.pnas.org/doi/10.1073/pnas.1620144114#fig01) using the same methods as the authors: *Nucleotide sequences were aligned using CLUSTALW, with ambiguous regions removed. Maximum likelihood trees with bootstrap support (1,000 replicates) were constructed using PhyML.*

To complete your analysis the key stages to consider are 1/ sequence alignment, 2/ the tree inference method, 3/ tree visualization and 4/ checking the reliability of clustering with bootstrapping.


**Software**

You can use the alignment software CLUSTALW by typing ‘clustalw2’ on the bioinformatics server Alpha2 <alpha2.cvr.gla.ac.uk>. Alternative alignment software you can try includes Muscle or Mafft (type ‘muscle’ or ‘mafft’ on the command line to see options). 

To use PhyML for tree inference, type ‘phyml’ on the command line (or available online at http://www.atgc-montpellier.fr/phyml/). See PhyML’s online helpfile for further guidance on options: http://www.atgc-montpellier.fr/phyml/usersguide.php. Note, PhyML takes PHYLIP formatted alignments which can be generated with CLUSTALW. Alternative phylogenetic software to try includes RaXML and IQ-TREE (commands: 'raxml-ng-mpi'and 'iqtree2'). 

FigTree (‘figtree’) is useful for visualizing phylogenetic trees and highlighting specific variants. 

Alignments, tree methods and visualization can also be carried out with graphical user interface software, for example, SeaView (command: ‘/software/seaview-v5.0.5/seaview’) or UGENE (‘ugene’). 



**Data**

The data set from the Iyer paper is quite large (available at /home4/VBG_data/Phylogenetics on Alpha) so will take some time to align. You can use the fasta file with fewer sequences (95) from the linked patients CH595 and CH455 that were presented in their figure 1. **Copy files to your own directory.**


Once you’ve generated some trees answer the questions below:

**Question 1**. What can you infer from your evolutionary tree about the relationship of virus from the two individuals: CH596 and CH455? What two properties of the phylogenetic tree support this relationship?


**Question 2**. What commands have you used to replicate the Figure 1A tree in Iyer et al, 2017? Why is PhyML, despite being a maximum likelihood method, relatively quick? What improvement would using the software RAxML bring to the analysis?


**Question 3**. Try a different method of tree inference such as the distance-method neighbor joining (available in CLUSTALW or SeaView). Does the result change in any meaningful way? What are the main differences between this method and maximum likelihood? What about the alignment method? 


**Question 4**. Does a different alignment method, e.g., MUSCLE or MAFFT, give the same result? Are there potential issues with using CLUSTALW (hint, what version of CLUSTAL does SeaView implement)?


**Question 5**. Have a go implementing bootstrapping. What does bootstrapping do? How does it contribute to the analysis? 


**Question 6**. Why does the substitution model used matter? What does the software jModelTest used in the Iyer paper do?

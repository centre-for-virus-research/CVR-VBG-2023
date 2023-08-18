---
layout: default
title: QualityAssessment
rank: 6
---


# Quality assessment of assemblies

*Mastomys natalensis* is a multimammate rodent from Sierra Leone that lives in close proximity to humans, sometimes transmitting Lassavirus to humans. In a large-scale investigation of the transmission risks of Lassavirus in Sierra Leone, different rodent species were captured/marked/released around several villages. Samples of blood (BL), from urogenital swabs (US) and oral swabs (OS) were collected and tested for Lassavirus. Lassavirus positive samples were RNA extracted and deep sequenced.  
In this practical, we will assess four metagenomic assemblies (BL, OS, US and combined) from two different rodents (sample 00094 and 00187).

***
> ## QUAST
> 
> Quality Assessment Tool for Genome Assemblies – QUAST/metaQUAST
> 
> Availability: <http://cab.cc.spbu.ru/quast/download.html>
> 
> Documentation: <http://cab.cc.spbu.ru/quast/manual.html>
***

1) Create a folder called validation and create links for each contig file in `/home4/VBG_data/QualityAssessment/`: 

```mkdir validation```

```ln -s /home4/VBG_data/QualityAssessment/Contigs_00094/00094_* .```

2) Run metaQUAST with the basic command

```metaquast.py -o metaquast_00094 00094_US_L_NA_S9_contigs.fasta 00094_OS_L_NA_S2_contigs.fasta 00094_BL_L1_NA_S15_contigs.fasta  00094_all_samples_contigs.fasta```

3) The output of metaQUAST will provide you with a path to the report in multiple formats. We will open report.html using firefox from the command (replace /path/to/report.html with the path printed out by metaQUAST):

```firefox /path/to/report.html```

__*Which sample provided the best assembly? Why metrics supports this? Why doesn’t the “all_samples” assembly have the largest genome?*__


<br/><br/>
***
> ## NCBI eutilities
> 
> Useful set of commands for interfacing with NCBI directly from the command line
> 
> Installation and documentation: <https://www.ncbi.nlm.nih.gov/books/NBK179288/>
> 
> Helpful examples: <https://bioinformatics.cvr.ac.uk/ncbi-entrez-direct-unix-e-utilities/>
***
<br/>

4) It is also possible to run QUAST/metaQUAST with a reference. As we are expecting to find Lassavirus, lets download Lassavirus sequences (two segment L and S) directly from the command line using NCBI eutilites:

`esearch -db nucleotide -query "OM735985" | efetch -format fasta > Lassa_L.fasta`

`esearch -db nucleotide -query "OM735984" | efetch -format fasta > Lassa_S.fasta`

5) This time I'll let you find the commands yourself to link to the different contig files for sample 00187, they are in `/home4/VBG_data/QualityAssessment/Contigs_00094/`. The syntax is the same as in 1).

6) Run the metaQUAST command, this time using as a reference Lassa_L.fasta and Lassa_S.fasta with the `-R` argument:

`metaquast.py -l "all, US, OS, BL" Contigs_00187/00187_all_samples_contigs.fasta Contigs_00187/00187_US_L_NA_S26_contigs.fasta Contigs_00187/00187_OS_L_NA_S24_contigs.fasta Contigs_00187/00187_BL_L1_NA_S5_contigs.fasta -R Lassa_L.fasta,Lassa_S.fasta`

7) Look at the report using firefox as before

__*Do you see a difference to the previous file? What additional metrics and visualizations are provided?*__

8) Try to run the metaQuast on the sample 00094 using the Lasavirus sequences as a reference with the `-R`.

Are there any error/warning messages during the run? Does the output look any different? What do you think is happening with this sample?

We will investigate in the metagenomic parctical what is going on with this sample.

<br/><br/>
**BONUS**

9) __*Use the different assemblies generated during the de-novo assembly practical to figure out which assembler has produced the best results? Which assembler covered the largest proportion of the genome? Which assembler produced the largest number of contigs?*__


---
layout: default
title: Metagenomics
rank: 7
---

# Metagenomics - contig classification

In the quality assessment practical, we noticed that the sample 00094 did not find a match to the Lassavirus sequences we provided as a reference. This could mean that Lassavirus isn’t in the sample so we will start by classifying the contigs using DIAMOND.

<br/>

***
> ### DIAMOND
> DIAMOND is a BLASTX equivalent but with at least 500 times speed performance. It is an aligner for protein and translated DNA searches. 
>
> Availability and documentation: <https://github.com/bbuchfink/diamond>
>
> Example to build a database:
> ```diamond makedb --in nr.faa -d nr```
>
> Example to run an alignment
> ```diamond blastx -d nr -q reads.fna -o matches.m8```

***
<br/>

1) Create a practical directory called `classify` using the command `mkdir` and then change directory `cd` into the new directory.

2) We will run DIAMOND against a database of all refseq viral proteins (this is a small database that will run within a few minutes). We need to create a temporary folder for diamond, before doing the DIAMOND blast.

```ln -s /home4/VBG_data/QualityAssessment/Contigs_00094/00094_all_samples_contigs.fasta .```

```mkdir contigs_temp_dir```

```diamond blastx -d /db/diamond/ViralRefSeqProtein.dmnd2 -p 2 -q 00094_all_samples_contigs.fasta -o 00094_all_samples_viralrefseq.txt -t contigs_temp_dir --top 3 --outfmt 6```


3)	Now, we will use the output with KronaTools to produce the visualization. Click the “color by e-value” checkbox. 

```ktImportBLAST 00094_all_samples_viralrefseq.txt -o 00094_all_samples_viralrefseq.html -tax /db/kronatools/taxonomy```

__*Explore the Kronaplot. Are there any matches to Lassavirus or anything similar to Lassavirus? How many contigs would you expect for Lassavirus? How many contigs do you have? What size do you expect these contigs to be?*__


4)	I have prepared in advance the comparison of the contigs to a recent version of Genbank NR using DIAMOND. This took several hours to run. Copy it to your local directory.  
 
```cp /home4/VBG_data/Metagenomics/00094_all_samples_nr.txt .``` 

5)	Run the ktImportBLAST command with this new input file (see step 3) and open it in firefox. 


__*What are the differences in the viral hits? Why is this virus not a Mammarenavirus in the kronaplot? Why do you think the database you blast against is important?*__

For more information about these samples and the discovery of this new arenavirus: <https://doi.org/10.1128/mra.00095-22>


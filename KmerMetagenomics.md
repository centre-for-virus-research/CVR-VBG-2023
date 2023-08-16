
# Kmer based metagenomics
## [6th Viral Bioinformatics and Genomics Training Course](https://github.com/centre-for-virus-research/CVR-VBG-2023)
* Monday 21st - Friday 25th August 2023
* Glasgow, UK
* [MRC-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/)

## Contact

[Richard Orton](https://www.gla.ac.uk/schools/infectionimmunity/staff/richardorton/)   
[MRC-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/)  
464 Bearsden Road  
Glasgow  
G61 1QH  
UK  
E-mail: Richard.Orton@glasgow.ac.uk  

## Contents

This practical is associated with a lecture on Kmer based metagenomics

* [1: Setup](#1-setup)
* [2: Removal of human host reads](#2-removal-of-human-host-reads)	 
* [3: Run kraken2](#3-run-kraken2)
* [4: Kraken2 on your own](#4-kraken2-on-your-own)

# 1: Setup

**Make sure you are logged into the alpha2 server with MobaXterm.**

In this session, we will be working with Illumina metagenomic data sets that have already been published and are available for download on the [NCBI Short Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra), there are two samples:

* Human
* Vampire bat

We will be using a tool called [Kraken2](https://ccb.jhu.edu/software/kraken2/) to analyse the paired end FASTQ reads for each sample. Kraken2 is the newest version of Kraken, a taxonomic classification system using exact k-mer matches to achieve high accuracy and fast classification speeds. This classifier matches each k-mer within a query sequence (i.e. a read or contig) to the lowest common ancestor (LCA) of all genomes within the database containing the given k-mer.

We will be broadly following the published Kraken metagenomic protocol:

**Metagenome analysis using the Kraken software suite**  
Lu et al. (2022). Nature Protocols volume 17, pages 2815–2839  
[https://www.nature.com/articles/s41596-022-00738-y](https://www.nature.com/articles/s41596-022-00738-y)  

We will be using the "standard" kraken database - which was downloaded from the [Kraken database download page](https://benlangmead.github.io/aws-indexes/k2) and is dated 20230605. The standard database includes RefSeq archaea, bacteria, viruses, plasmid complete genomes, UniVec Core and the most recent human reference genome, GRCh38.

To set things up, first change directory (cd) to your home directory:

```
cd
```

Then copy (cp) the data folder (-r for recursive as we want the folder and all it's contents) to your current directory (which will be your home directory after the above command was entered):

```
cp -r /home4/VBG_data/Kmer .
```

Then change directory to the Human data folder (within the Kmer folder):

```
cd ~/Kmer/Human/
```

Next, list the contents of the directory so you can see the files we will be working with:

```
ls
```

You should see the FASTQ paired-end read files:

**SRR533978\_1.fq**  
**SRR533978\_2.fq**

# 2: Removal of human host reads

Although the kraken database does contain human data, it is good practice to remove the human reads from the sample as human data is confidential and it will also speed up later steps.

We will use a tool called [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align the reads to the human genome, and then extract the unmapped reads to analyse further.

The human reference genome bowtie2 index has already on alpha2 (from [https://benlangmead.github.io/aws-indexes/bowtie](https://benlangmead.github.io/aws-indexes/bowtie)) so we do not need to download or build anything.


```
bowtie2 -x /home4/VBG_data/hg19/hg19 -1 SRR533978_1.fastq -2 SRR533978_2.fastq -S human.sam -p 8
```

***Command breakdown:***

1. **bowtie2** = the name of the program we are executing
2. **-x /home4/VBG_data/hg19/hg19** = the name (and location) of the INDEXED reference genome to align to
3. **-1 SRR533978_1.fastq** = the name of read file 1
4. **-2 SRR533978_2.fastq** = the name of read file 2
5. **-S human.sam** = the name of the output SAM file to create
6. **-p 8** = use 8 computer threads

When finished bowtie2 will produce some summary statistics to the command prompt. If we list the contents of the directory, you should now see the SAM file called **human.sam** has been created:

```
ls
```

We can now extract the non-human reads (which are the unmapped reads - exploiting the SAM flag 4) using the fastq funciton of [samtools](http://www.htslib.org):

```
samtools fastq -1 nonhuman_1.fastq -2 nonhuman_2.fastq -f 4 -s singleton.fastq human.sam 
```

***Command breakdown:***

1. **samtools** = the name of the program we are executing
2. **fastq** = the name of the function within samtools we are using
3. **-1 nonhuman_1.fastq** = output the read1 of each pair in the file nonhuman_1.fastq
4. **-2 nonhuman_2.fastq** = output the read2 of each pair in the file nonhuman_2.fastq
5. **-f 4** = only include read alignments that do have the unmapped flag 4
6. **-s singleton.fastq** = output singleton reads (those without a pair) into the file singleton.fastq
7. **human.sam** = the name of the input SAM file to extract reads from

**NB:** Just to note, SAM/BAM files need to be sorted by the read name when extracting paired end data (so that pair members are next to each in the file) - the SAM file outputted by bowtie2 has pairs next to each already so we can skip the sort step

If we list the contents of the directory, you should now see the two nonhuman FASTQ files have been created:

```
ls
```

Lets check that they actually contain some data by counting the number of lines in them:

```
wc -l nonhuman*
```

If you see millions of lines in each file we can safely delete our human.sam file as we no longer need it:

```
rm human.sam
```

# 3: Run kraken2


```
/software/kraken2-v2.1.1/kraken2 --db /home4/VBG_data/k2_standard_20230605 --threads 8 --minimum-hit-groups 3 --report-minimizer-data --paired nonhuman_1.fastq nonhuman_2.fastq --output kraken_output.txt --report kraken_report.txt 
```

***Command breakdown:***

1. **/software/kraken2-v2.1.1/kraken2** = the name and location (path) of the program we are executing
2. **--db /home4/VBG_data/k2_standard_20230605** = the name (and location) of the kraken2 database we are using
3. **--threads 8** = use 8 computer threads
4.  **--report-minimizer-data** = forces Kraken 2 to provide unique k-mer counts per classification
5. **--minimum-hit-groups 3** = for increased classification precision (i.e., fewer false positives).
6. **--paired nonhuman_1.fastq nonhuman_2.fastq** = the name of the paired input FASTQ files
7. **--output kraken_output.txt** = the name of the kraken output file to create
8. **--report kraken_report.txt** = the name of the kraken report output file to create

**NB:** The --minimum-hit-groups flag specifies the minimum number of ‘hit groups’ needed to make a classification call. Hit groups are overlapping k-mers sharing the same minimizer. Kraken 2 uses minimizers to compress the input genomic sequences, thereby reducing storage memory needed and run time. In this example, we increase the minimum number of hit groups from the default two groups to three groups for increased accuracy. See [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual)

Overall, this command will output two files - our kraekn_output.txt file and our kraken_report.txt file. The kraken_output.txt is large and not really human readable whereas the kraken_report.txt is a human readable summary that is tab-delimited, with one line per taxon. The fields of the output, from left-to-right, are as follows:

1. Percentage of reads covered by the clade rooted at this taxon
2. Number of reads covered by the clade rooted at this taxon
3. Number of reads assigned directly to this taxon
4. Number of minimizers in read data associated with this taxon (new)
5. An estimate of the number of distinct minimizers in read data associated with this taxon (new)
6. A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
7. NCBI taxonomy ID
8. indented scientific name

You can explore the report file manually using a command like **more**, esepcially if you understand the taxonomy. However, a more visual way of viewing the results is via a Krona plot:

```
ktImportTaxonomy -q 2 -t 3 -s 4 kraken_output.txt -o kraken_krona.html -tax /db/kronatools/taxonomy
```
***Command breakdown:***

1. **ktImportTaxonomy** = the name of the program we are executing
2. **-q 2** = column 2 of input files to use as query ID.
3. **-t 3** = column 3 of input files to use as taxonomy ID
4. **-s 4** = column 4 of input files to use as score. 
5. **kraken_output.txt** = the name of the input file - whcih is the kraken output file
6. **-o kraken_krona.html** = the name of the Krona html file to create
7. **-tax /db/kronatools/taxonomy** = the location of the NCBI taxonomy files for Krona

This will create an new html output file which you should see in your directory:

```
ls
```

We can open this file with Firefox:

```
firefox kraken_krona.html
```

***
### Questions
**Question 1** – What viruses have the highest read counts in the sample?
***

**NB:** Alternatively, the html file can be downloaded via the MobaXterm file browser on the left hand side of the window onto your local machine and opened there

### What is a krona plot?

Krona is an interactive visualization tool for exploring the composition of metagenomes within a Web browser. A Krona plot is a html file that can be opened by a web browser (such as Firefox, Chrome, Safari and Internet Explorer/Edge). Krona uses multilevel pie charts to visualize both the most abundant organisms and their most specific classifications (Fig. 1). Rather than hiding lower ranks in its overview, Krona hides low-abundance organisms, which can be expanded interactively. 

**Krona: Interactive Metagenomic Visualization in a Web Browser**  
Ondov et al. (2013)  
DOI 10.1007/978-1-4614-6418-1_802-1  
[https://link.springer.com/content/pdf/10.1007/978-1-4614-6418-1_802-1.pdf](https://link.springer.com/content/pdf/10.1007/978-1-4614-6418-1_802-1.pdf)

An example Krona plot is shown below:

![](https://github.com/rjorton/OIE2023/blob/main/krona.png)

Some of the key features of a Krona plot:

1. You can double click on any taxon (i.e. any slice of the pie) to view only that taxon and it's children
2. Clicking/Selecting a taxon will display how many sequences (i.e. reads or contigs) have been assigned to that taxon in the top right corner under "count" (this is actually a hyperlink to a file containing all the sequence IDs assigned to this taxon). The taxonomy ID is also displayed as a hyperlink to the corresponding NCBI page.
3. The defauly colouring of the krona plot is set to distinguish taxons apart - however the colouring can be changed based on the 'score' (for BLAST results this corresponds to e-value) by selecting the "Color by Avg. log e-value" - this would highlight the taxons with the highest score


## 3.1: Bas congo virus / Tibrovirus congo

The Human sample came from the following paper where deep sequencing was used to discover a novel rhabdovirus (Bas-Congo virus, or BASV) associated with a 2009 outbreak of 3 human cases of acute hemorrhagic fever in Mangala village, Democratic Republic of Congo (DRC), Africa.

**A Novel Rhabdovirus Associated with Acute Hemorrhagic Fever in Central Africa**  
Grard et al. 2012  
PLoS Pathog. 2012 Sep; 8(9): e1002924.  
[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3460624/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3460624/)

The viral species is now called 'Tibrovirus congo': [https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1987017](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1987017)

# 4: Kraken2 on your own

If you have time, try processing the second sample the vampire bat and reusing and adapting the commands for the human sample - although there is no need to align the human genome first.

The sample is one sample from the following paper:

**Using noninvasive metagenomics to characterize viral communities from wildlife**  
Molecular Ecology Resources Volume19, Issue 1, January 2019, Pages 128-143  
Bergner et al. 2018  
[https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12946](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.12946)

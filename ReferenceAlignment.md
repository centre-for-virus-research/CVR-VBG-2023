# Reference Alignment Practical
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

This practical is associated with a lecture on Reference Alignment of High-Throughoput Sequencing (HTS) reads to a reference sequence.

* [0: Overview](#0-overview)
* [1: Setup](#1-setup)
	+ [1.1: Basic read statistics](#11-basic-read-statistics)
	+ [1.2: Read quality control](#12-read-quality-control)
* [2: Read Alignment](#2-read-alignment)
	+ [2.1: Indexing the reference sequence](#21-indexing-the-reference-sequence)
	+ [2.2: Aligning the reads to the reference](#22-aligning-the-reads-to-the-reference)
	+ [2.3: Converting SAM to BAM](#23-converting-SAM-to-BAM)
	+ [2.4: Basic alignment statistics](#24-basic-alignment-statistics)
 	+ [2.5: Coverage plot](25-coverage-plots)	 
* [3: Alignment on your own](#3-alignment-on-your-own)
* [4: Consensus calling](#4-consensus-calling)
 	+ [4.1: Consensus calling on your own](#41-consensus-calling-on-your-own)
* [5: Variant calling](#5-variant-calling)
 	+ [5.1: Variant calling on your own](#51-variant-calling-on-your-own)
* [6: Extra data](#6-extra-data)
* [7: Assembly visualisation with tablet](#7-assembly-visualisation-with-tablet)
* [8: Glossary](#8-glossary)
* [9: Bash script example](#9-bash-script-example)
 
# 0: Overview

**YOU DO NOT NEED TO ENTER THE COMMANDS IN THIS OVERVIEW SECTION!**

In this practical, we will be aligning paired end reads to a reference sequence. Commands that you need to enter into the Terminal window (i.e. the command line) are presented in a 'code' box with a different font, like this:

```
ls
```

Sometimes a command can be long and may not be fully visible (but the code box is within a scrollpane to aid viewing), but the command should still be entered as one single command on the computer terminal. 

```
bwa mem -t 4 my_really_long_reference_filename.fasta my_really_log_read_file_1.fastq my_really_long_read_file_2.fastq > my_really_long_output_file.sam
```

A few Linux tips to remember:

1.	Use the **Tab button** to automatically complete filenames – especially long ones
2.	Use the **Up Arrow** to scroll through your previous commands, it enables you to easily re-run or re-use/change/correct old commands
3.	**Case matters**, the following file names are all different:

```
Myfile.txt
MyFile.txt
MYFILE.txt
myfile.txt
my file.txt
my_file.txt
```

Also watch out for number 1s being confused with lowercase letter L’s, and capital O’s being confused with zeroes

```
l = lower case letter L
1 = number one
O = capital letter O
0 = zero
```

# 1: Setup

**Make sure you are logged into the alpha2 server with MobaXterm.**

In this session, we will be working with two sets of Illumina paired end reads which were simulated from the viral genomes of two different SARS-CoV-2 samples; these simulated reads were created using ART (Huang et al., 2012: [10.1093/bioinformatics/btr708](10.1093/bioinformatics/btr708)). The goal now is to align these reads to a reference genome sequence, with an ultimate goal of creating a consensus sequence for mutation anlysis.

To start off, you will need to copy the data we need for the practical to your home directory. First change directory (cd) to your home directory:

```
cd
```

Then copy (cp) the data folder (-r for recursive as we want the folder and all it's contents) to your current directory (which will be your home directory after the above command was entered):

```
cp -r /home4/VBG_data/Richard .
```

Then change directory to the Sample1 data folder (within the Richard folder):

```
cd ~/Richard/Sample1/
```

Next, list the contents of the directory so you can see the files we will be working with:

```
ls
```

You should see the FASTQ paired-end read files:

**S1\_R1.fq**  
**S1\_R2.fq**

The reference file that we will be using is located in the Refs folder ~/Richard/Refs (whose relative path is ../Refs/):

```
ls ../Refs
```

You should see a file called:

**sars2_ref.fasta**  

## 1.1: Basic read statistics

We will first use a tool called [prinseq](https://prinseq.sourceforge.net) to count the number of reads in each file. As these are paired end reads, there should be one read from each read pair in each file – and hence the same number of reads in each file. We will also use prinseq to output statistics on the read length (just to point out - prinseq itself can do much much more!):

```
prinseq-lite.pl -stats_info -stats_len -fastq S1_R1.fq -fastq2 S1_R2.fq
```

***Command breakdown:***

1.	**prinseq-lite.pl** is the name of the program
2.	**-stats\_info** tells prinseq to output basic stats on the reads (number of reads and bases)
3.	**-stats\_len** tells prinseq to output basic stats on read lengths (min, max, mean etc)
4.	**-fastq S1\_R1.fq** the name of the 1st FASTQ file
5.	**-fastq2 S1\_R2.fq** the name of the 2nd FASTQ file in the pair

### Common Issue
* A common issue here is not entering the prinseq command on one line in the terminal - you should only use the enter key at the end of the command to execute it.
* Another common issue is typos - check the command carefully if you get an error - it is likely you have mispelled a file or argument (or even the prinseq program name)

***
### Questions
**Question 1** – How many reads and bases are in the read files 1 and 2?

**Question 2** – What is the average (mean) length of the reads? 
***

The prinseq statistics are split into those for the first FASTQ file of the read pair (e.g. stats\_info, stats\_len, etc) and those for the second FASTQ file of the read pair (e.g. stats\_info2, stats\_len2, etc), and should look a bit like this (but with different numbers!):

```
stats_info	bases		48000000
stats_info	reads		320000
stats_info2	bases		48000000
stats_info2	reads		320000
stats_len	max		150
stats_len	mean		150.00
stats_len	median		150
stats_len	min		150
stats_len	mode		150
stats_len	modeval		320000
stats_len	range		1
stats_len	stddev		0.00
stats_len2	max		150
stats_len2	mean		150.00
stats_len2	median		150
stats_len2	min		150
stats_len2	mode		150
stats_len2	modeval		320000
stats_len2	range		1
stats_len2	stddev		0.00 
```

Paired read files should always have the same number of lines/reads (the ordering of the reads in each file is also critical), so if your two paired files have a different number of reads, something has gone wrong (e.g. filtering/trimming went wrong and corrupted the output, or maybe files from different samples are being used). 

## 1.2: Read quality control

The reads have already been quality filtered and trimmed so we can move on to alignment (and read QC was covered in a previous session).
 
# 2: Read Alignment

There are many tools available to align reads onto a reference sequence: [BWA](http://bio-bwa.sourceforge.net), [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [minimap2](https://github.com/lh3/minimap2), [bbMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/), to name but a few.

We will be using [BWA](http://bio-bwa.sourceforge.net) to align our paired end reads to a reference sequence and output a [SAM (Sequence Alignment Map)](https://samtools.github.io/hts-specs/SAMv1.pdf) file. The SAM file contains the result of each read’s alignment to the given reference sequence. 

## 2.1: Indexing the reference sequence

First, we need to create a BWA index of the reference sequence. Tools such as BWA need to index the sequence first to create a fast lookup (or index) of short sequence seeds within the reference sequence. This enables the tools to rapidly align millions of reads:

```
bwa index ../Refs/sars2_ref.fasta
```

If you list (ls) the contents of the ../Refs directory, you should see the BWA index files, they will all have the prefix sars2\_ref.fasta, and will have extensions such as **.amb**, **.ann**, **.bwt**, **.pac**, and **.sa**.

```
ls ../Refs/
```

## 2.2: Aligning the reads to the reference

Next, we want to align our reads to the reference sequence using the BWA mem algorithm:

```
bwa mem -t 4 ../Refs/sars2_ref.fasta S1_R1.fq S1_R2.fq > S1.sam
```

***Command breakdown:***

1. **bwa** = the name of the program we are executing
2. **mem** = the BWA algorithm to use (recommended for illumina reads > 70nt)
3. **-t 4** = use 4 computer threads
4. **../Refs/sars2\_ref.fasta** = the name (and location) of the INDEXED reference genome to align to
5. **S1\_R1.fq** = the name of read file 1
6. **S1\_R2.fq** = the name of read file 2
7. **>** = direct the output into a file
8. **S1.sam** = the name of the output SAM file to create 

Overall, this command will create an output file called S1.sam in the current directory, which contains the results (in SAM format) of aligning all our reads to the reference sequence sars2\_ref.fasta.

When bwa has finished (and your prompt comes back), check that the SAM file has been created.

```
ls
```

There should now be a file called **S1.sam** in the directory.

### Common issue
A common mistake is not waiting for your previous command to finish, and entering the next command into the terminal before the prompt has returned. You need to wait until the **username@alpha2** command prompt returns before entering the next command - the bwa alignment can sometimes take a few minutes.

## 2.3: Converting SAM to BAM

Typically, a SAM file contains a single line for each read in the data set, and this line stores the alignment result of each read (reference name, alignment location, CIGAR string, the read sequence itself, quality, etc).

SAM files are in a text format (which you can open and view if you like, e.g.: **head S1.sam**), but can take up a lot of disk storage space. It is good practice to convert your SAM files to BAM (Binary Alignment Map) files, which are compressed binary versions of the same data, and can be sorted and indexed easily to make searches faster. We will use [samtools](https://samtools.github.io) to convert our SAM to BAM, and sort and index the BAM file:

```
samtools sort S1.sam -o S1.bam
```

```
samtools index S1.bam
```

***Command breakdown:***

1.	The first command tells samtools to **sort** the SAM file, and to also output (**-o**)the sorted data in BAM format to a file called **S1.bam**. This sorts all the alignmed reads by alignment start position, placing all the reads whose alignment starts at reference position 1 first, then those that start of position 2, then 3, 4, 5 and so on.
3.	We then use samtools to **index** the BAM file S1.bam (indexing [which relies on sorted data] enables faster searches downstream).


There should now be two new files in the directory called: 

**S1.bam** (the BAM file)  
**S1.bam.bai** (the BAM index file) 

Now let’s list (ls) the contents of the directory to check we have our new files, and also check out their sizes:

```
ls -lh
```

***Command breakdown:***
* **-l** tells the list (**ls**) command to give the output in a long list format, whilst the **h** tells it to provide file sizes in a human readable format, this is the 5th column, which will have the size of each file in a format such as 2.5M (M for megabytes) or 9.5G (G for gigabytes).

***

### Questions
**Question 3** – How big is the SAM file compared to the BAM file?

***

**NB:** If your SAM file is 0B (i.e. 0 bytes, empty) then something went wrong with the bwa alignment step, so restart from there. If you SAM file is fine (i.e. >0), but your BAM file is 0B (i.e. empty), then something went wrong with your SAM to BAM conversion so re-do that step. 

We don’t need our original SAM file anymore (as we have the BAM file now) so we remove (rm) the SAM file S1.sam:

```
rm S1.sam
```

## 2.4: Basic alignment statistics

One common thing to check is how many reads have been aligned (or mapped) to the reference, and how many are not aligned (or unmapped). Samtools can report this for us easily, utilising the aligner SAM flags you learnt about in the previous session.

**Reminder:** the 2nd column in the SAM file contains the flag for the read alignment. If the flag includes the number 4 flag in its makeup then the read is unmapped, if it doesn’t include the number 4 in it's makeup then it is mapped.

### Number of unmapped reads
```
samtools view -c -f4 S1.bam
```

***Command breakdown***

1.	**samtools view** = to view the file S1.bam
2.	**–c** = count the read alignments
3.	**–f4** = only include read alignments that do have the unmapped flag 4

### Number of mapped read alignments:
```
samtools view -c -F4 S1.bam
```

***Command breakdown***

1.	**samtools view** = to view the file S1.bam
2.	**–c** = count the read alignments
3.	**–F4** = skip read alignments that contain the unmapped Flag 4 

***
### Questions

**Question 4** – how many reads are mapped to the sars2_ref.fasta genome? 

**Question 5** – how many reads are unmapped?
***

Technically, the above command gives the number of mapped read **alignments** not reads. A read could be mapped equally well to multiple positions (one will be called the primary alignment, and others secondary alignments [sam flag 256]), or a read could be split into two parts (e.g. spliced) with one part being the primary alignment and the others supplementary [sam flag 2048]

So to get the true number of mapped reads you need to count only the alignments that do not have flags 4 (unmapped), 256 (not primary), and 2048 (supplementary) = 4 + 256 + 2048 = 2308

### Number of mapped reads

```
samtools view -c -F4 -F256 -F2048 S1.bam
```

or summing up the F flag values together:

```
samtools view -c -F2308 S1.bam
```

For small RNA viruses, secondary and supplementary alignments tend to be rare, but it is important to know the distinction between mapped **reads** and mapped read **alignments**.

## 2.5: Coverage plots

We previously used samtools to count the number of mapped and unmapped reads (using samtools view -c commands), now let’s explore the read mapping in more detail by creating a coverage plot using a tool called weeSAM: (https://github.com/centre-for-virus-research/weeSAM)[https://github.com/centre-for-virus-research/weeSAM]

weeSAM analyses a SAM or BAM file, generates a graphical coverage plot, and reports a range of summary statistics such as:

* **Ref_Name**: The identifier name of the reference.
* **Ref_Len**: The length in bases of each reference.
* **Mapped\_Reads**: Number of read alignments mapped to each reference.
* **Breadth**: The number of sites in the genome covered by reads.
* **%\_Covered**: The percent of sites in the genome which have coverage.
* **Min\_Depth**: Minimum read depth observed.
* **Max\_Depth**: Max read depth observed.
* **Avg\_Depth**: Mean read depth observed.
* **Std\_Dev**: Standard deviation of the mean (Avg_Depth).
* **Above\_0.2_Depth**: Percentage of sites which have greater than 0.2 * Avg_Depth.
* **Above\_1_Depth**: Percentage of sites which are above Avg_Depth.
* **Above\_1.8_Depth**: Percentage of sites which have greater than 1.8 * Avg_Depth.
* **Variation\_Coefficient**: The mean of Std_Dev of the mean.

The Average Depth (Avg_Depth) is perhaps the most important field, along with Breadth which will tell you how much of the genome is covered by aligned reads. But the fields such as Std\_Dev and Above_0.2_Depth can give an indication of the variability in the coverage across the genome.

Let’s run weeSAM on our sample:

```
weeSAM --bam S1.bam --html S1
```

An explanation of this command is:

1.	**weeSAM**: the name of the program we are using
2.	**--bam**: flag to signify input bam file
3.	**S1.bam**: the name of our bam file to analyse
4.	**--html**: flag to signify output html file
5.	**S1**: the name prefix to use for the output

If you list the contents of the directory you should see that a folder called **S1\_html\_results** has been created:

```
ls
```

Inside this folder is a HTML file that we can view in a web browser (like Firefox or Chrome), the HTML file has the summary statistics and coverage plot so lets take a look and open the html file.

```
firefox S1_html_results/S1.html
```

***RJO CHECK IN COMPUTER ROOM - will this definitiely launch firefox via MobaXterm or should they download locally?***

You should see something like this:

![](https://github.com/rjorton/OIE2023/blob/main/weeSAM_table.png)

***
### Questions
**Question 6** – what is the average depth of coverage across the SARS-CoV-2 reference genome?
***

Now let’s view the coverage plot by clicking on the hyperlink (blue and underlined) in the Ref_Name column, you should see a coverage plot similar to this:

![](https://github.com/rjorton/OIE2023/blob/main/weeSAM_fig.png)

The x-axis represents the genome position, whilst the y-axis represents the Depth of Coverage at each genome position. 

**NB:** The reference sequence filename is sars2_ref.fasta, but the actual name of the sequence itself is [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) in the fasta file, you can open up the file yourself to check this if you want (head –n1 sars2_ref.fasta).

**Close the weeSAM/Firefox windows before proceeding!**

### Common issue
A common issue here is due to the fact that we have launched firefox from the terminal (wihtout running it background - see advanced linux commands). In order to get our command prompt back (the username@alpha2) we need to close the firefox window down, the prompt should then return.


# 3: Alignment on your own

You now need to use BWA to align the reads for Sample2 (S2_R1.fq and S2_R2.fq) to the sars2_ref.fasta reference sequence. So lets move into the correct folder:

```
cd ../Sample2
```

**You need to work out the commands yourself based on the previous commands for the Sample1**. Here is a reminder of the commands you used for Sample1 (S1) which you will need to adapt:

```
prinseq-lite.pl -stats_info -stats_len -fastq S1_R1.fq -fastq2 S1_R2.fq
```

```
bwa mem -t 4 ../Refs/sars2_ref.fasta S1_R1.fq S1_R2.fq > S1.sam
```
```
samtools sort S1.sam -o S1.bam
```
```
samtools index S1.bam
```
```
rm S1.sam
```
```
samtools view -c -f4 S1.bam
```
```
samtools view -c -F2308 S1.bam
```
```
weeSAM --bam S1.bam --html S1
```

**NB:** This is a new sample, but the steps are exactly the same. Essentially, you will want change the names of your input FASTQ filenames and the output files (e.g. from S1 to S2) in each command.

***
### Questions

**Question 7** – how many reads are in the original FASTQ files and how many reads are mapped to the sars2_ref.fasta genome for Sample2?

**Question 8** – how many reads are unmapped? what is the average coverage and breadth of coverage?
***

# 4: Consensus calling

We have now aigned each of our samples (S1 and S2) to the Wuhan-Hu-1 (NC_045512.2) SARS-CoV-2 reference genom sequence, and now we want to call a consensus sequence.

What is a consensus sequence? At each genome position in the SAM/BAM alignment file, call the most frequent nucleotide (or insertion/deletion) observed in all of the reads aligned at the position. 

In this practical, we will use a tool called [iVar](https://andersen-lab.github.io/ivar/html/manualpage.html) to call the consensus sequence, which utilises the [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) function of samtools.

**NB:** the server we are using (alpha2) has a conflict/error when running ivar by default, to resolve it please **COPY** and **PASTE** the below command into your terminal window and hit the enter button:

```
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/htslib-v1.12/lib
```

First, let work on Sample1, so we need to change directory (cd) into the correct folder:

```
cd ~/Richard/Sample1
```

And now call the consenus for the sample using iVar:

```
samtools mpileup -aa -A -d 0 -Q 0 S1.bam | ivar consensus -p S1 -t 0.4
```

Breaking this command down, there are two parts:

1. samtools [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) which essentially outputs the base and indel counts for each genome position
	* **-aa** = output data for absolutely all positions (even zero coverage ones)
	* **-A** = count orphan reads (reads whose pair did not map)
	* **-d 0** = override the maximum depth (default is 8000 which is typically too low for viruses)
	* **-Q 0** = minimum base quality, 0 essentially means all the data
2. ivar [consensus](https://andersen-lab.github.io/ivar/html/manualpage.html) - this calls the consensus - the output of the samtools mpileup command is piped '|' directly into ivar
	* -p S1 = prefix with which to name the output file
	* -t 0.4 = the minimum frequency threshold that a base must match to be used in calling the consensus base at a position. In this case, an ambiguity code will be used if more than one base is > 40% (0.4). See [iVar manual](https://andersen-lab.github.io/ivar/html/manualpage.html)

By default, iVar consensus uses a minimum depth (-m) of 10 and a minimum base quality (-q) of 20 to call the consensus; these defaults can be changed by using the appropriate arguments. If a genome position has a depth less than the minimum, an 'N' base will be used in the consensus sequence by default.

iVar will output some basic statistics to the screen such as:

```
#DO NOT ENTER THIS - IT IS AN EXAMPLE OF AN IVAR OUTPUT:
Minimum Quality: 20
Threshold: 0.4
Minimum depth: 10
Regions with depth less than minimum depth covered by: N
[mpileup] 1 samples in 1 input files
[mpileup] Max depth set to maximum value (2147483647)
Reference length: 29903
Positions with 0 depth: 0
Positions with depth below 10: 4
```

and when it has finished (and your prompt returns) you should see our consensus sequence (S1.fa) in the directory:

```
ls
```

which you can view the sequence via the command line (we will be covering variants later):

```
more S1.fa 
```

***
### Questions

**Question 9** - try copying and pasting the created consensus sequence into [NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) - what is the closest sample on GenBank?
***

## 4.1: Consensus calling on your own

You now need to call the consensus sequence for Sample2, so you'll need to change directory to the appropriate folder and adapt the ivar command for the Sample2 files.

# 5: Variant calling

Viruses, and in particular RNA viruses, can exist as complex populations consisting of numerous variants present at a spectrum of frequencies – the so called viral quasispecies. Although we have created a consensus sequence (which typically considers mutations at a frequency >50% in the sample) using iVar, we do not yet know anything about the mutations within the sample. Furthermore, it is often necessary to go beyond the consensus, and investigate the spectrum of low frequency mutations present in the sample.

iVar itself could be used to call variants (using the [iVar variants](https://andersen-lab.github.io/ivar/html/manualpage.html) command). But here we will be using a slightly more advanced variant caller called [LoFreq](https://github.com/CSB5/lofreq) to call the low (and high) frequency variants present in the sample BAM file. LoFreq uses numerous statistical methods and tests to attempt to distinguish true low frequency viral variants from sequence errors. It requires a sample BAM file and corresponding reference sequence that it was aligned to, and creates a [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file as an output.

First, lets make sure we are in the correct folder to work on Sample1:

```
cd ~/Richard/Sample1
```

To use LoFreq enter this command:

```
/software/lofreq_star-v2.1.2/bin/lofreq call -f ../Refs/sars2_ref.fasta -o S1.vcf S1.bam
```

Breaking this command down:

* **/software/lofreq_star-v2.1.2/bin/lofreq**: the name (and path on alpha2) of the program we are using
* **call**: the name of the function within LoFreq we are using – call variants
* **-f**: ../Refs/sars2_ref.fasta: the reference file name and location (path)
* **-o S1.vcf**: the output VCF file name to create
* **S1**.bam: the input BAM file name

Now lets open the VCF file created by LoFreq:

```
more S1.vcf
```

The outputted VCF file consists of the following fields:

* CHROM: the chromosome – in this case the SARS-CoV-2 ref sequence NC_045512.2
* POS: the position on the chromosome the variant is at
* ID: a ‘.’ but LoFreq can be run with a database of known variants to annotate this field
* REF: the reference base at this position
* ALT: the alternate base (the mutated base) at this position
* QUAL: LoFreq’s quality score for the variant
* FILTER: whether it passed LoFreq’s filters e.g. PASS
* INFO: Detailed Information
	* DP=1248; depth = 1248
	* AF=0.995192; Alt Frequency (Mutation Frequency) = 99.5%
	* SB=0; Strand Bias test p-value
	* DP4=0,1,604,638: Coverage of the ref base in Fwd and Rev, and the alt base in Fwd and Rev

***
### Questions

**Question 10** – how many consenus level (i.e AF > 0.5) and subconsenus (i.e. AF < 0.5) are there in the sample? what genome positions are the sub-consensus mutations?
***

LoFreq simply calls the reference position and mutation present, it does not characterise the effect of the mutation in terms of whether the mutation is synonymous or nonsynonymous etc. To do that we will use a program called [SnpEff](https://github.com/pcingola/SnpEff) which is run on LoFreq’s outputted VCF file and creates a new annotated VCF file:

```
/software/snpEff/scripts/snpEff -ud 0 NC_045512.2 S1.vcf > S1_snpeff.vcf
```
**RJO - At the moment this gives a permsission error - need Scott to run it to so that it downloads the reference DB for use by everyone else - but have duplicated into ../Scripts/ in case**

Breaking this command down:

* **/software/snpEff/scripts/snpEff**: the name (and location) of the program we are using
* **-ud 0**: Set upstream downstream interval length to 0
* **NC_045512.2**: the reference file name
* **S1.vcf**: the input vcf file name
* **S1_snpeff.vcf**: the annotated output vcf file name

Setting -ud 0 stops SnpEff from characterising mutations located near (but not within) a gene as being located in their UTRs. Typically viral genomes are compact and genes are separated by few bases, in such cases a mutation in one gene could also be characterised as being in the UTR region of a neighbouring gene (as SnpEff was initially built for human analyses) – try running SnpEff without the -ud 0 and compare the results if you want, you should see multiple annotations for each mutation.

Now we can view the annotated vcf file created by SnpEff:

```
more S1_snpeff.vcf
```

The mutations will now have annotations added at the end of the Info field (the last field, the 8th field), e.g in Sample1 you should see this nonsynonymous (missense) mutation at position 23,403 which corresponds to D614G a Asp to Gly mutation at codon 614 in the Spike gene of SARS2:

```
DO NOT ENTER THIS - IT IS AN EXAMPLE OF A MUTATION IN THE VCF!!!
DP=1286;AF=0.995334;SB=0;DP4=1,1,628,652;ANN=G|missense_variant|MODERATE|S|GU280_gp02|transcript|GU280_gp02|protein_coding|1/1|c.1841A>G|p.Asp614Gly|1841/3822|1841/3822|614/1273||
```

The A to G mutation at genome position 23403 corresponds to position 1841 (out of 3822) within the Spike (S) gene which corresponds to codon 614 (out of 1273) within Spike(S)

**NB:** SnpEff includes many pre-built databases (SARS-CoV-2 NC_045512.2 being one of them) – for different viruses you may need to build the SnpEff database first by downloading and processing a GenBank file, see the documentation [here](https://pcingola.github.io/SnpEff/)


## 5.1: Variant calling on your own

Now you task is to run LoFreq and then SnpEff to characterise the mutations in Sample 2?

***
### Questions
**Question 11** - how many consensus level non-synonymous mutations are there in Sample?
***

# 6: Extra Data

If you are looking for something extra to do, there are additional data sets located in the folder:

### ~/Richard/Ebola/

You will find a set of (gzipped) FASTQ paired end read files, and a reference FASTA sequence to align them to.

The reads are from a patient from the ebola epidemic in West Africa 2014 {Gire et al, 2014} [https://www.ncbi.nlm.nih.gov/pubmed/25214632](https://www.ncbi.nlm.nih.gov/pubmed/25214632)

The reference ebola sequence is from a 2007 outbreak in Democratic Republic of Congo. 

Try aligning the reads to the reference yourself.

### ~/Richard/Dengue

This is a simulated Dengue virus sample, but we do not know what genotype it is and therefore what reference alignment to use. Align the paired end FASTQ reads to the two provided reference sequences (genotype1 and genotype3), count the number of mapped reads and create a coverage plot to determine which genotype is the best reference.

### ~/Richard/Mystery/

This is a mystery sample, combine all the given references sequences in the folder into one file using the “cat” command, align the reads to that combined reference (after indexing) and then determine what the virus in the sample is.

# 7: Assembly Visualisation with Tablet

[Tablet](https://ics.hutton.ac.uk/tablet/) is a tool for the visualisation of next generation sequence assemblies and alignments. It goes beyond simple coverage plots, and allows you to scroll across the genome, zoom into errors of interests, highlight mutations to the reference, and investigate the assembly.

Tablet requires three files:

1.	A bam file, e.g. 1a.bam
2.	A bam index file, e.g. 1a.bam.bai
3.	Optional: A reference sequence file: e.g. sars2\_ref.fasta

**Tablet demonstration**

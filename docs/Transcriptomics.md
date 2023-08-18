---
layout: default
rank: 9
---

# RNA-Seq Workshop
## Quan Gu

A classical RNA-Seq data processing pipeline contains the following steps:

• Quality control \
• Mapping RNA-Seq reads against a reference genome \
• Visualizing reads and transcript structures \
• Performing differential expression analysis \
• Visualizing differential expression analysis \
• Functional annotation and pathway analysis  

In this course, all the data we are using has been subsampled to save time and space.

###  1.	Log on to the CVR server ### 

Open MobaXterm and log on to the CVR-Alpha2 server.

Step 1: Create a working directory and go to your directory using the commands:

```
mkdir RNASeq2023
cd RNASeq2023
```
  
This workshop use single-end strand-specific RNA-Seq reads from six samples (6 fastq files), located in **/home4/VBG_data/RNASeq**

You will generate soft links to these files in your current directory. 

Condition1 - Mock samples: **Mock01.fastq; Mock02.fastq; Mock03.fastq**

Condition2 - Interferon treated samples: **IFNb01.fastq; IFNb02.fastq; IFNb03.fastq**


Step 2: Generate a soft link to the files, here is an example for one file:

```
ln -s /home4/VBG_data/RNASeq/IFNb01.fastq
```

We need to make soft links to the reference genome file as well.

```
ln -s /home4/VBG_data/RNASeq/Human
ln -s /home4/VBG_data/RNASeq/Homo_sapiens.GRCh38.107.gtf
```

**Task**: Have a look at the human reference genome (GRCh38) and its annotation files, and print the last 50 lines of these files. 

###  2. Quality control ### 

The first step of RNA-Seq analysis is Quality control. 
Before we start the analysis, it is always good inspect the quality of the samples. Here we use **FastQC** to visualize the quality of the sequences. To begin, simply run:
```
fastqc 
```
The software use from the terminal will bring up a GUI window for you on your desktop.  Load a sample from the GUI.

**File -> Open -> IFNb01.fastq**

We could do trimming after quality control. Trimming includes quality trimming and adapter trimming.

Several trimming programs exist: **Trim_galore**, **cutadapt** , **Prinseq**, **Trimmomatic**,etc.  
In this workshop, we will use **Trim_galore** to trim the reads. 
Note: the default adapter sequence of Illumina and Proton sequencing platforms are different.

The following command will trim with default settings:

```
trim_galore IFNb01.fastq >/dev/null 2>&1
```

**Task**: Check the difference between raw data and trimmed data by (1) file size, (2) read number and (3) read quality.


###  3.	Mapping RNA-Seq reads against a reference genome ###  

The next step is to align reads to the genome. The mapping output are BAM files which is a binary version of a SAM file. A SAM file (.sam) is a tab-delimited text file that contains sequence alignment data. This step is most time-consuming step.
In this workshop, we will use **Hisat2** to map the reads to the reference genome, which is a mapping program specifically developed for RNA-Seq analysis.  When doing this at home, remember to build the **Hisat2** index before you run the mapping, which I have already done for you. I have already done this for you and you can use it with the following command:

```
genome='/home4/VBG_data/RNASeq/Human'

hisat2 -x $genome --rna-strandness R -U IFNb01_trimmed.fq -S IFNb01.sam
```

This has produced a SAM file. We need to generate the BAM files as well, which is the binary version of SAM files. 

 ```
samtools view -bh IFNb01.sam > IFNb01.bam

samtools sort -o IFNb01_sorted.bam IFNb01.bam

samtools index IFNb01_sorted.bam

rm IFNb01.sam
```

### 4.	Visualizing reads and transcript structures ### 
After getting your provided aligned BAM files, you can visualize the location of the mapped transcripts using visulization software such as: **IGV**, **Tablet**, **Ugene**, etc.

### 5.	Performing differential expression analysis ###
In this workshop, we will use edgeR for the differential expression(DE). Here, we omit the step of transcriptome assembly as we don’t want to discover novel DE genes from those in our current genome.
There are two tools which can count the mapped transcripts in the genome. To learn more about how **featureCounts** and **htseq-count** compare, I have written a blog about it http://bioinformatics.cvr.ac.uk/blog/featurecounts-or-htseq-count/

Here we use **featureCounts** to get the mapped raw counts of each gene in each sample.

```
featureCounts -s 2 -a /home4/VBG_data/RNASeq/Homo_sapiens.GRCh38.107.gtf -o temp.txt IFNb01_sorted.bam

cat temp.txt| cut -f 2,3,4,5,6 --complement|sed 's/\t/ /g'|awk 'NR==2;NR>2{print $0|"sort -n"}' > IFNb01_count.txt

rm temp.txt

```

In this course, we don't have enough time to make all the count table by ourselves. However, we can use the overall counts file that I have prepared for you.

Please copy the count file **DE1input.txt** into your own directory. 

```
cp /home4/VBG_data/RNASeq/DE1input.txt  . 
more DE1input.txt
```

After, we can open R and run **edgeR** to carry out DE analysis.

I have written an Rscript called **edgeR** that has all the necessary commands. You can view the code of **edgeR** with the following command:

```
more /home4/VBG_data/RNASeq/edgeR.r
```

Here we use the classical mode of DE analysis in **edgeR**, which assumes all replicates within group are equal.

You can copy and paste the content from **edgeR.r** into **R** after opening R with the command ```R```, or you can simply type from the command line:

```
Rscript /home4/VBG_data/RNASeq/edgeR.r
```

Then you will get the output files: **DEG_edgeR.csv**, **cpm.csv** , **bcvplot.pdf**, **VolcanoPlot.png** and **mdsplot.pdf**.


### 6. Function annotation and pathway analysis ### 
Once we get the differential expression gene list, the next step is to annotate them and
analyse their function and pathways. There are many free online tools available. 
One of the
best-known tools is **David** (https://david.ncifcrf.gov). 
We could also try **Reactome**
(https://reactome.org), which is highly recommended. 
The user guide of **Reactome** is here (https://reactome.org/userguide/analysis).

Alternatively, you could directly run the Reactome R script I have made for you. Then you will get the output files: **gene.bubble.pdf**, **gene.SigReactome.xls**. The pdf file is the figure of the top 20 pathways enrichment of our significant DE genes. The xls file is the table of details of all the pathways (with P-value < 0.05) of our significant DE genes.  

```
cat DEG_edgeR.csv|sed 's/"//g'| awk 'NR>1' | cut -f1 -d',' > gene.glist
Rscript /home4/VBG_data/RNASeq/Reactome.r hs gene.glist ./ 20
```



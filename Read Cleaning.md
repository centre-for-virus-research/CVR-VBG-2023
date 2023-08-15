**Read Cleaning Practical**

Derek Wright, MRC-University of Glasgow Centre for Virus Research

[**Derek.Wright@glasgow.ac.uk**](mailto:Derek.Wright@glasgow.ac.uk)

**Overview**

In this practical, we will be: 

- performing quality control analysis using FASTQC
- collating the output from FASTQC into a single HTML report with MultiQC
- trimming adapters with Trim Galore

**Linux Commands**

In the practical, commands that you need to enter into the Terminal window (i.e. the command line) are presented in a box and with a different font, like this:

ls


` `A few tips to remember:

1) Use the **Tab button** to automatically complete filenames – especially long ones.
1) Use the **Up Arrow** to scroll through your previous commands, it enables you to easily re-run or re-use/adapt old commands.
1) **Case Matters** - the following file names are all different:

Myfile.txt

MyFile.txt

MYFILE.txt

myfile.txt

my file.txt

my\_file.txt


1) Watch out for number 1 being confused with lowercase L, and capital O being confused with zero.

l = lower case L

1 = number one

O = capital letter O

0 = zero


1) ~ is shorthand for your home directory
1) . is shorthand for the current working directory
1) .. is shorthand for the directory above
1) \* may be used as a wildcard to match file names

**Dataset**

First let’s copy over a folder of data to analyse:


cp -r /home3/dw73x/RawReads ~/


The directory contains FASTQ files from a small run of reads from the Illumina MiSeq System. Some of these samples were subsequently run at much higher output and the data led to the discovery of a new arenavirus published here:

<https://doi.org/10.1128/mra.00095-22>

**NGS Cleaning Practical**

1. **FASTQC**

First, we will examine the data set, so let’s move into that directory:


cd ~/RawReads


Next, launch the interactive version of FASTQC to examine the quality of the reads.


fastqc &


This launches FASTQC with the graphical user interface (GUI). 

*Tip:* The ampersand tells the shell to run a program in the background, meaning that you can still enter shell commands while the program is running.

From FASTQC:	

- Select **File > Open** from the menu bar
- Navigate to the RawReads directory and select all of the *\*.fastq.gz* files in the directory, holding down the control key to select multiple files
- Click the **Open** button

Each dataset is opened in its own tab. Green ticks in the sidebar for the quality metrics mean that quality control has passed. Orange exclamation marks are warnings. Red crosses are fails.



1. **FASTQC with MultiQC**

You could save each report individually, but an alternative workflow is to run FASTQC from the command line, instead of from the GUI, outputting a FASTQC report for each sample and then collating the reports into a single report using MultiQC.


mkdir reports

fastqc -o reports \*.fastq.gz

cd reports

multiqc .


This generates an HTML report named *multiqc\_report.html*

The report is meant to be hosted on a web server or downloaded to your computer then viewed in a web browser. 

I have placed a copy of the report for you to view at:

<http://bioinformatics.cvr.ac.uk/multiqc_report.html>

1. **Trim Galore**

Trim Galore is a wrapper around the programs Cutadapt and FASTQC, aimed at providing convenient trimming and QC without setting complex parameters. 

Trim Galore provides options for various adapters or can attempt to auto-detect adapter sequences. 

The option –-paired may be used for paired-end reads, which have two files for each sample. 

The option –-fastqc may be added to additionally run FASTQC in a single step.



cd ..

trim\_galore --paired \*.fastq.gz


Now list the directory contents:


ls


You should see new trimmed sequence files that have been created by Trim Galore, with names ending in *\_trimmed.fq.gz*

---
layout: default
title: De Novo Practicals
rank: 5
---

# *De Novo* Assembly Practical
---

### Author

- Sreenu Vattipally
- MRC-University of Glasgow Centre for Virus Research
- University of Glasgow
- G61 1QH, UK

---

In this practical, we will use three different *de novo* assemblers, i.e. [Spades](https://github.com/ablab/spades), [IDBA_UD](https://github.com/loneknightpy/idba) and [ABySS](https://github.com/bcgsc/abyss), to assemble the SARS-CoV-2 genome from Illumina reads. This data is generated from a [clinical sample](https://www.ncbi.nlm.nih.gov/sra/?term=SRR21065613) using an amplicon sequencing strategy.

Let us first create the project directory and download the required data.

Create a directory for the practical. 
```
mkdir -p ~/DeNovo/SARS-CoV-2 
```

Change the directory
```
cd ~/DeNovo/SARS-CoV-2
```

Download the data from the [SRA database](https://ncbi.nlm.nih.gov/sra) using fastq-dump (part of [SRA toolkit](https://github.com/ncbi/sra-tools))
```
fastq-dump --split-files SRR21065613
```

(You might have to update vdb-config by running `vdb-config --interactive`). 

The command line option `--split-files` preserves the paired-end data and creates the corresponding files. Otherwise, it combines the reads and saves them as a single file. 

Do a quality check and trimming using the [trim_galore](https://github.com/FelixKrueger/TrimGalore) program.
```
trim_galore --illumina --length 50 -q 30 --paired SRR21065613_1.fastq SRR21065613_2.fastq
```

After running trim_galore,  you should see new files (`SRR21065613_1_val_1.fq` and `SRR21065613_2_val_2.fq`) and trimming report files in your working directory. Next, create individual directories for Spades, IDBA_UD and Abyss.

```
mkdir Spades Abyss IDBA_UD
```

## De novo assembly using Spades 

First, go to the Spades directory you have created above
```
cd ~/DeNovo/SARS-CoV-2/Spades
```

Run spades assembly program with different k-mer sizes.
```
spades.py --careful -k 21,45,73,101 --only-assembler -1 ../SRR21065613_1_val_1.fq -2 ../SRR21065613_2_val_2.fq -o .
```

Here,

- --careful : Tries to reduce the number of mismatches and indels in the assembly   
- -k: different k-mer sizes
- --only-assembler: No error correction.
- -1: first fastq reads file
- -2: second fastq reads file
- -o: output location. 

Here we are saving results in the current directory(.), we can specify any name.

The program will take some time to complete the job depending on the CPU power and available memory. If everything goes without any error, the final assembled contigs will be stored in `contigs.fasta` file. In this file, contigs are saved in a sorted manner, with the longest contig saved first. Remember, these contigs can be positive or negative strands. We have to check the orientation, comparing them with the reference sequence.

### TODO: 
Extract all the headers of contigs from the `contigs.fasta` file. 
What is the length of the longest contig generated?

---

## De novo assembly using IDBA_UD

Now let us try running IDBA_UD. Go to the IDBA_UD directory you created earlier. 
```
cd ~/DeNovo/SARS-CoV-2/IDBA_UD/
```

Though it is an assembly program for high throughput reads, IDBA_UD works only with fasta files (not fastq). It will not consider the quality of the reads while assembling. So... we have to convert reads from fastq to fasta and submit them to IDBA_UD.

Paired-end reads should be converted to a single fastA file. In the file, a pair of reads are in two consecutive lines. To do this, We are using `fq2fa` (this program is part of the IDBA package) to merge two FastQ reads files into a single fastA file. Following, we are running the assembly with multiple kmers.

```
fq2fa --merge ../SRR21065613_1_val_1.fq ../SRR21065613_2_val_2.fq reads.fa 
```

This will convert fastq reads to fasta sequences and saves them in a single file (`reads.fa`)

Let us run the assembly.
```
idba_ud --mink=21 --maxk=81 step=10 -r reads.fa -o . --num_threads=4
```


Here,

- --mink: minimum k-mer size
- --maxk: maximum k-mer size
- --step: k-mer size increment
- -r: input reads
- -o: output directory
- --num_threads: number of CPU threads to use

#### Caution:
If you run IDBA_UD with single-end reads or generate the `reads.fa` file by other means, at the end of the run, you might see

```
invalid insert distance 
Segmentation fault (core dumped)
```

This is because IDBA_UD could not calculate the insert distance between the two reads. However, the assembly actually did finish. Using single-end rather than paired-end reads means that your assemblies will NOT have scaffolds.

The final results will be stored in `contig.fa` file.

### TODO: 
Again, extract all the headers of contigs from `contig.fa` file. What is the length of the longest contig generated?


---

## De novo assembly using ABySS

Go to the Abyss directory

```
cd ~/DeNovo/SARS-CoV-2/Abyss/
```

Run Abyss assembler

```
abyss-pe k=51 n=4 in='../SRR21065613_1_val_1.fq ../SRR21065613_2_val_2.fq' name=Abyss-51 B=2G
```

Here,

- k: k-mer size
- n: number of threads to use (adjust this value as per your system capacity)
- in: paired-end input reads
- name: output name
- B: Unconditionally make all targets

As you notice from above, the abyss assembler works with one k-mer at a time. Since we do not know the optimal k-mer size for our data, we should use various k-mer sizes. Please use the below for loop to run the program with different k-mers.

```
for k in $(seq 41 10 91);
do
abyss-pe k=$k n=4 in='../SRR21065613_1_val_1.fq ../SRR21065613_2_val_2.fq' name=Abyss-$k B=2G
done
```

This will run abyss with k-mer sizes from 41 to 91 with the increment of 10 (we are using seq command here) and save the output in corresponding files. Here you might have noticed that each k-mer size produces different output. Let us combine all the contigs.

```
cat Abyss*-unitigs.fa > Abyss-contigs.fa
```

To get the assembly stats, run

```
abyss-fac *-contigs.fa
```

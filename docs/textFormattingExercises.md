## Exercises for Unix text formatting commands.

#### Copy the data to your local directory

`cp -R /home3/bvv2t/bashExp ~/.`

1. Remove multiple new lines from the file `exercise.txt`
2. Create a complementary genome sequence of `refGenome.fa`
3. Print the genome composition (number of different nucleotides) and length of `refGenome.fa`
4. Extract ORFs from `refGenome.fa` taking coding coordinates from the file `refGenomeFeatureTable.txt` 
5. Convert the FASTQ file `file.fq` to the FASTA file and save it as `file.fa`
6. Extract all the reads `file.fq` with `GAATTC` and convert and save the output in a fasta file. 

##### Please save your commands into individual files (call then example-1.txt, example-2.txt, etc...) and try running them from the terminal. What do you observe? 

### Answers.
1. ` tr -s '\n' < exercise.txt `
2. `grep -v ">" refGenome.fa |tac|tr -d '\n'|rev |tr 'ATGC' 'TACG'`
3. `grep -v ">" refGenome.fa |grep -o . |sort |uniq -c ` and `grep -v ">" refGenome.fa |grep -o . |cat -n |tail -1 |cut -f1`
4. `grep -v ">" refGenome.fa |tr -d '\n'|cut -c 10-20 `
5. `paste - - - - < file.fq |sed s/^@/\>/|cut -f 1,2|sed s/\\t/\\n/g` or `grep -A 1 "^@" file.fq |sed s/@/\>/g|grep -v -- --`
6. `grep -B 1 GAATTC file.fq |sed s/@/\>/g|grep -v -- --`

## Exercises for BASH scripting.

`cp -R /home3/bvv2t/bashExp-2 ~/.`


1. Write a script to convert a fastq file to a fasta file (Use file.fq).
2. Write a script to print the sequence header and sequence length from a multi-fasta file (use refGenomes.fa)
3. Print reverse complement of a multi-fasta sequence file (use refGenomes.fa)
4. Write a reference mapping pipeline (Use reference alignment input files from earlier sessions).
  - Take the input fastq and reference from the command line.
  - Clean the reads using trim_galore
  - Index the reference sequence with the bwa index
  - Map the reads to reference using bwa mem.
  - Convert bam to consensus file.
  - Get the mapping statistics


### Answers.
### Convert FastQ to FastA
```
ind=1;
while read line
do
    if [ $((ind%4)) -eq 1 ]]; then
    echo $line | sed s/^@/\>/g
    elif [[ $((ind%4)) -eq 2 ]; then 
    echo $line
    fi
    let idx++
done < $1
```

### Print Genome name and size
```
exec < $1
while read line
do
        if [[ ${line:0:1} != ">" ]];
        then
                genome=$genome$line;
        else
                # Do not print first name 
                if [[ -n $name ]]; then echo $name ${#genome}; fi
                name=$line;
                genome="";
        fi
done
# Print last name 
echo $name ${#genome};
```

### Print reverse complementary sequence
```
while read line
do
        if [[ ${line:0:1} != ">" ]];
        then
                genome=$genome$line;
    else

# Do not print first name. -n is TRUE is the length of the string is non-zero
        if [[ -n $name ]]; then 
            echo $name
            echo $genome|rev|tr '[ATGC]' '[TACG]'|fold -w 70
        fi
        name=$line;
        genome="";
        fi
done < $1

# Print last name 
echo $name
echo $genome|rev|tr '[ATGC]' '[TACG]'|fold -w 70
```

### Script for reference mapping and consensus generation 

```
#!/usr/bin/bash

# bash script to map all the read files against all the reference sequences.
# It convets bam to consensus sequence also, generates mapping statistics

# This script assumes
#	1. Reads are sequences by Illumina
#	2. Reads are paired-end
#	3. All the reads are stored in "FQs" directory and have ".fastq.gz" extension
#	4. All the references are stored in Ref directory and have ".fa" extension 
#	
#	Written by: Sreenu Vattipally
#

for fq in FQs/*_1.fastq.gz
do 
	# This have to be adjusted as per the file name
	name=$(basename -s "_1.fastq.gz" $fq); 

	# Trim the reads using trim_galore program
	trim_galore --paired --illumina --length 75 -q 30 FQs/${name}_1.fastq.gz FQs/${name}_2.fastq.gz

	# Scan all the ref sequences
	for ref in Ref/*.fa
	do
		refName=$(basename -s ".fa" $ref); 
		# Map the reads and take only mapped reads. Index the reference sequences for the first time
		#bwa index $ref
		bwa mem -t 4 $ref ${name}_1_val_1.fq.gz ${name}_2_val_2.fq.gz|samtools view -F4 -bS |samtools sort -o $name-$refName.bam

		# Create consensus
		samtools consensus $name-$refName.bam -o $name-$refName-consensus.fa

		# Get alignment stats
		samtools coverage $name-$refName.bam -o $name-$refName-stats.txt
	done
done
```


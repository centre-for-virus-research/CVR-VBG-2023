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

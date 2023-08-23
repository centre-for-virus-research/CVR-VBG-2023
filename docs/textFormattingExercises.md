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


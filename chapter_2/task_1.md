# Chapter Two
## Task 1 - Evaluating the Quality of Illumina Data
From your terminal windows, navigate to the 'sequencing_data' directory (you may be there already) and list the contents of the directory.
```bash
cd genomics_adventure/sequencing_data$

ls -lath
```
[IMAGE]

These are paired-end read data that we downloaded in the introduction. Many programs require that the read 1 and read 2 files have the reads in the same order. You can identify reads from the same pair because they will have the same header followed by either a "1" or a "2".  We will now look at the raw reads. To view the first few headers we can use the 'zcat', 'head' and 'grep' commands:
```bash
zcat ERR2789854_1.fastq.gz | head | grep @HWI
zcat ERR2789854_1.fastq.gz | head | grep @HWI
```
[IMAGE]

The only difference in the headers for the two reads is the read number. Of course this is no guarantee that all the headers in the file are consistent. To get some more confidence, repeat the above commands using 'tail' instead of 'head' to compare reads at the end of the files. 
```bash
zcat ERR2789854_1.fastq.gz | tail | grep @HWI
zcat ERR2789854_1.fastq.gz | tail | grep @HWI
```
You can also check that there is an identical number of reads in each file using cat, grep and wc –l:
```bash
zcat ERR2789854_1.fastq.gz | grep @MISEQ | wc -l
zcat ERR2789854_2.fastq.gz | grep @MISEQ | wc -l
```
Now, let's run the fastqc program on the data. Unlike the QC lab, we will open up a Graphical User Interface (GUI) and load the data this way. To do this, run: 

# Task 6 - Generating an Index
Now let's change to our *E. coli* reference sequence directory, do you remember where it is? Can you get there in one command? Also let's remind ourselves what is actually in the directory.
```bash
cd ~/genomics_adventure/reference_sequences/ecoli

ls -lath
```
[IMAGE]

In this directory we have 2 files:
 1. GCF_000005845.2_ASM584v2_genomic.fna which is a FASTA file which contains the reference genome sequence
 2. GCF_000005845.2_ASM584v2_genomic.gff which is a file that contains the annotation for this genome.

## BWA
Go ahead, type 'bwa' in your terminal and see what happens. Hopefully, you should see something similar to this:

[IMAGE]

BWA is actually a suite of programs which all perform different functions. We are only going to use two during this workshop, 'bwa index' and 'bwa mem'. Explore the outputs of these two commands now.

By default bwa index will use the IS algorithm to produce the index. This works well for most genomes, but for very large ones (e.g. vertebrate) you may need to use bwtsw. For bacterial genomes the default algorithm will work fine.

### Create a Reference Index
```bash
bwa index GCF_000005845.2_ASM584v2_genomic.fna
```
Take a look at the output of thi command in your terminal (e.g. 'ls -lath'). You will notice that the BWA index program has created a set of new files. These are the index files that bwa mem will need.

# [Task 6]()

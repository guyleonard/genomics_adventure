# Task 4 - Assembly Time!

Firstly let's construct an assembly using only the available Illumina data in a new directory!

```bash
cd ~/workshop_materials/genomics_adventure

mkdir pseudomonas && cd pseudomonas
```

The next command, like most of the assemblies in this adventure (and in real life) take some time to complete, so we have precomputed the results for you. But the command is given for your reference.

:warning: do not run this next command :warning:
```bash
spades.py --threads 2 --careful -o illumina_only \
-1 ../sequencing_data/pseudomonas_gm41/SRR491287_1_val_1.fq.gz \
-2 ../sequencing_data/pseudomonas_gm41/SRR491287_2_val_2.fq.gz
```

Let's find the precomputed data and have a look at it.
```bash
# create a symbolic link to the precomputed data
ln -s ../precomputed/pseudo_illumina_only illumina_only
```

You know the score! Let's run QUAST!!
```
cd illumina_only

quast.py --output-dir quast contigs.fasta

cat quast/report.txt
```

If you ran the assembly command - outside of the workhop - you may get slightly different results here, as SPAdes uses a random seed.
```
Assembly                    contigs
# contigs (>= 0 bp)         612    
# contigs (>= 1000 bp)      323    
# contigs (>= 5000 bp)      253    
# contigs (>= 10000 bp)     205    
# contigs (>= 25000 bp)     88     
# contigs (>= 50000 bp)     29     
Total length (>= 0 bp)      6642052
Total length (>= 1000 bp)   6589395
Total length (>= 5000 bp)   6395291
Total length (>= 10000 bp)  6038360
Total length (>= 25000 bp)  4113769
Total length (>= 50000 bp)  2097787
# contigs                   339    
Largest contig              187340 
Total length                6601889
GC (%)                      58.98  
N50                         32206  
N75                         19064  
L50                         60     
L75                         127    
# N's per 100 kbp           0.00
```

# Go to [Task 5](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_5/task_5.md)


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
-1 ../sequencing_data/pseudomonas_gm41/SRR491264_1_val_1.fq.gz \
-2 ../sequencing_data/pseudomonas_gm41/SRR491264_2_val_2.fq.gz
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
# contigs (>= 0 bp)         4232   
# contigs (>= 1000 bp)      723    
# contigs (>= 5000 bp)      272    
# contigs (>= 10000 bp)     107    
# contigs (>= 25000 bp)     17     
# contigs (>= 50000 bp)     2      
Total length (>= 0 bp)      5429451
Total length (>= 1000 bp)   4268084
Total length (>= 5000 bp)   3074920
Total length (>= 10000 bp)  1957974
Total length (>= 25000 bp)  613336 
Total length (>= 50000 bp)  113579 
# contigs                   1122   
Largest contig              62496  
Total length                4531320
GC (%)                      62.74  
N50                         7861   
N75                         4079   
L50                         143    
L75                         344    
# N's per 100 kbp           0.00
```

# Go to [Task 5](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_5/task_5.md)


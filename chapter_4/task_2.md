# Checking the Assembly
Change in to the "assembly" directory and prepare to run QUAST. We did this before in Chapter 3 - Task 4...

```bash
cd assembly

quast.py --output-dir quast contigs.fasta
```

It won't take too long to run, and hopefully you have already guessed we are going to have a look inside the "report.txt" file next :thumbs_up:. It should look something like this:

```
Assembly                    contigs
# contigs (>= 0 bp)         195    
# contigs (>= 1000 bp)      100    
# contigs (>= 5000 bp)      73     
# contigs (>= 10000 bp)     65     
# contigs (>= 25000 bp)     52     
# contigs (>= 50000 bp)     30     
Total length (>= 0 bp)      4610398
Total length (>= 1000 bp)   4587656
Total length (>= 5000 bp)   4526852
Total length (>= 10000 bp)  4469157
Total length (>= 25000 bp)  4260579
Total length (>= 50000 bp)  3450238
# contigs                   112    
Largest contig              268612
Total length                4595430
GC (%)                      50.73
N50                         111731
N75                         54089
L50                         14
L75                         30
# N's per 100 kbp           0.00
```

You can see that there are 195 contigs in the assembly - so it is far from complete. NB - your results may differe slightly. The [N50](http://en.wikipedia.org/wiki/N50_statistic):mag: is 111K and the N75 is 54K which tells us that most of the assembly is in quite large contigs. What you see if fairly normal for a short read assembly - you shouldn't expect complete chromosomes.

A good check at this point is to map the original reads back to the 'contigs.fasta' file and check that all positions are covered by reads. Amazingly, it is actually possible for *de novo* assemblers to generatecontigs to which the original reads will not map!!

# Go to [Task 3](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_4/task_3.md)

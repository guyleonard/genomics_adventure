# Task 4 - Assesment of Assemblies
We will use the program called - Quality Assessment Tool for Genome Assemblies - [QUAST](http://quast.sourceforge.net/quast) to generate some statistics on the assembly.

```bash
quast.py --output-dir quast spades_assembly/contigs.fasta
```

This will create a directory called 'quast' and create some statistics on the assembly you produced, don’t worry if the results look a little different to the example.

Take a look at the 'report.txt' file:
```
cat quast/report.txt
```

Example output:
```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                    contigs
# contigs (>= 0 bp)         17     
# contigs (>= 1000 bp)      5      
# contigs (>= 5000 bp)      4      
# contigs (>= 10000 bp)     3      
# contigs (>= 25000 bp)     2      
# contigs (>= 50000 bp)     0      
Total length (>= 0 bp)      97030  
Total length (>= 1000 bp)   94891  
Total length (>= 5000 bp)   93528  
Total length (>= 10000 bp)  88160  
Total length (>= 25000 bp)  74276  
Total length (>= 50000 bp)  0      
# contigs                   5      
Largest contig              46850  
Total length                94891  
GC (%)                      47.95  
N50                         27426  
N75                         27426  
L50                         2      
L75                         2      
# N's per 100 kbp           0.00
```

Try to interpret the information in the light of what we were trying to do. Because we are assembling unaligned reads we are not expecting a whole chromosome to pop out. We are expecting bits of our strain that does not exist in the reference that we aligned against, possibly some contamination, and various small contigs made up of reads that didn't quite align to our reference.

The N50 and L50 measures are very important in a normal assembly and we will visit them later, they are not really relevant to this assembly.

You will notice that we have a couple of large contigs greater than 25kb long - what do you think they might be? Also ~12 other contigs longer than 1kb. We need to find out what they are!

# Go to [Task 5](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_3/task_5.md)

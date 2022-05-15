# Create a Hybrid Assembly
Now will execute the same command, but this time include the longer PacBio reads to see the effect it has on our assembly.

Return to the pseudomonas directory:
```bash
cd /home/zoo/zool2474/genomics_adventure/pseudomonas
```
Again, this assembly will take too long for now, but it's here for reference anyway.

:cat: :cat:
Woah déjà vu. :sunglasses:

:warning: do not run the assembly :warning:
```bash
spades.py --threads 112 --careful -o hybrid \
--pacbio ../sequencing_data/pseudomonas_gm41/SRR1042836_subreads.fastq.gz \
-1 ../sequencing_data/pseudomonas_gm41/SRR491264_1_val_1.fq.gz \
-2 ../sequencing_data/pseudomonas_gm41/SRR491264_2_val_2.fq.gz
```

Let's find the precomputed data and have a look at it.
```bash
# create a symbolic link to the precomputed data
ln -s ../precomputed/pseudo_hybrid hybrid
```

Well, you definitely know the score! Let's run QUAST!!
```
cd hybrid

quast.py --output-dir quast contigs.fasta

cat quast/report.txt
```

If you ran the assembly command - outside of the workhop - you may get slightly different results here, as SPAdes uses a random seed.
```
Assembly                    contigs
# contigs (>= 0 bp)         4270   
# contigs (>= 1000 bp)      750    
# contigs (>= 5000 bp)      279    
# contigs (>= 10000 bp)     103    
# contigs (>= 25000 bp)     15     
# contigs (>= 50000 bp)     1      
Total length (>= 0 bp)      5431133
Total length (>= 1000 bp)   4266386
Total length (>= 5000 bp)   3025550
Total length (>= 10000 bp)  1835035
Total length (>= 25000 bp)  511402 
Total length (>= 50000 bp)  51083  
# contigs                   1152   
Largest contig              51083  
Total length                4531494
GC (%)                      62.74  
N50                         7299   
N75                         4036   
L50                         154    
L75                         362    
# N's per 100 kbp           0.00
```

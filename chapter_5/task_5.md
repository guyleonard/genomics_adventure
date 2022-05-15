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
-1 ../sequencing_data/pseudomonas_gm41/SRR491287_1_val_1.fq.gz \
-2 ../sequencing_data/pseudomonas_gm41/SRR491287_2_val_2.fq.gz
```

Let's find the precomputed data and have a look at it.
```bash
# create a symbolic link to the precomputed data
ln -s ../precomputed/pseudo_hybrid hybrid
```

Well, you definitely know the score! Let's run QUAST!! But this time we want to compare the two assemblies!!

```
cd hybrid

quast.py --output-dir quast ../illumina_only/contigs.fasta contigs.fasta

cat quast/report.txt
```

If you ran the assembly command - outside of the workhop - you may get slightly different results here, as SPAdes uses a random seed.
```
Assembly                    illumina_only_contigs  hybrid_contigs
# contigs (>= 0 bp)         612                    265           
# contigs (>= 1000 bp)      323                    94            
# contigs (>= 5000 bp)      253                    84            
# contigs (>= 10000 bp)     205                    79            
# contigs (>= 25000 bp)     88                     62            
# contigs (>= 50000 bp)     29                     41            
Total length (>= 0 bp)      6642052                6684539       
Total length (>= 1000 bp)   6589395                6658678       
Total length (>= 5000 bp)   6395291                6632315       
Total length (>= 10000 bp)  6038360                6593847       
Total length (>= 25000 bp)  4113769                6281030       
Total length (>= 50000 bp)  2097787                5518221       
# contigs                   339                    96            
Largest contig              187340                 671305        
Total length                6601889                6660346       
GC (%)                      58.98                  58.97         
N50                         32206                  137540        
N75                         19064                  62570         
L50                         60                     15            
L75                         127                    32
```

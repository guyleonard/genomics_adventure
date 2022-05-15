# Task 3 - Trim the Reads!

Trim the Illumina reads as before (~10 mins):

```bash
trim_galore --paired --fastqc --gzip --cores 4 --length 100 SRR491287_1.fastq.gz SRR491287_2.fastq.gz
```

You can check the number of filtered reads using "grep –c" and the quality of trimmed reads with fastqc if you want.

For our next trick we want to keep the long reads from PacBio even though they are of lower quality. We are relying
on the assembler to use them appropriately...

# Go to [Task 4](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_5/task_4.md)

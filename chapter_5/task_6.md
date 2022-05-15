# Task 6 - Align Reads Back to Reference
Let's realign our original reads back to the assembly and see what we have - refer to previous notes if you are unsure of the steps.

Start in the pseudomonas directory, soft link the hybrid assembly and create an index with BWA.
```bash
cd ~/workshop_meterials/genomics_adventure/pseudomonas

mkdir mapping_to_assembly && cd mapping_to_assembly

ln -s ../hybrid/contigs.fasta

bwa index contigs.fasta
```

First map the Illumina reads, and then follow the standard protocols we learned previously...You can follow along, or make your own commands...
```bash
# map
bwa mem -t 2 contigs.fasta \
../../sequencing_data/pseudomonas_gm41/SRR491287_1_val_1.fq.gz \
../../sequencing_data/pseudomonas_gm41/SRR491287_2_val_2.fq.gz \
> pseudo_illumina.sam

# convert to bam
samtools view -bS pseudo_illumina.sam > pseudo_illumina.bam

# sort
samtools sort -o pseudo_illumina_sorted.bam pseudo_illumina.bam

# index
samtools index pseudo_illumina_sorted.bam

# stats
samtools flagstat pseudo_illumina_sorted.bam > pseudo_illumina_sorted.stats
```

We can now map our PacBio data to our assembly too!. For this we will use another tool called "minimap2" which is better suited to mapping PacBio data than BWA. Read more on this [here](https://lh3.github.io/2018/04/02/minimap2-and-the-future-of-bwa) :mag:
It too creates SAM files, and then we can follow the same procedure with the output as we do for Illumina data.
```bash
minimap2 -ax map-pb -t 2 contigs.fasta \
../../sequencing_data/pseudomonas_gm41/SRR1042836_subreads.fastq.gz > pseudo_pacbio.sam

# convert to bam
samtools view -bS pseudo_pacbio.sam > pseudo_pacbio.bam

# sort
samtools sort -o pseudo_pacbio_sorted.bam pseudo_pacbio.bam

# index
samtools index pseudo_pacbio_sorted.bam

# stats
samtools flagstat pseudo_pacbio_sorted.bam > pseudo_pacbio_sorted.stats
```

Have a look at both flagstat files and give them a comparison.

# Go to [Task 7](https:a//github.com/guyleonard/genomics_adventure/blob/release/chapter_5/task_7.md)

# Task 7 - SAM File Manipulation - A Bit Less Basic...
Before we can visualise the alignment in a more meaningful way, we need to convert the SAM file into a BAM file (Binary AlignMent format). SAM files are great for human readability, but they are not so good for fast computational access. The binary format will allow speedy access to the information stored within it and also reduce the file size.

To do this we will use another suite of programs called '[samtools](http://www.htslib.org/)' :mag:. Go ahead and type the command below, and have a look at the options.
```bash
samtools view
```

## Converting SAM to BAM
We can see that we need to provide 'samtools view' with a reference genome as a FASTA formatted file (-T), the '-b' and '-S' options to indicate that the output should be in BAM format, and that the input is in SAM format, all along with the alignment input file.
```bash
samtools view -bS \
-T ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
ecoli_mapped.sam > ecoli_mapped.bam
```

This should take around 2 minutes. Note that for larger datasets you may wish to set multiple threads as well with the --threads option. It's always good to check that your files have processed correctly, if something goes wrong it's better to
catch it immediately. Note that the 'BAM' file is much smaller than the 'SAM' file - this is to be expected as the binary format is more efficient.

## Sorting a BAM File
Once the conversion is complete we will need to sort the BAM file, this is so that the reads are stored in the order they appear along the chromosomes or scaffolds or contigs. We can do this using the 'samtools sort' command. There are a couple of different ways to sort a BAM file, and is left to the user to decide which is best. For the next task we will use the default options and it will take about 2 minutes to complete.
```bash
samtools sort ecoli_mapped.bam -o ecoli_mapped_sorted.bam
```

## A Note on Efficiency
In the previous set of tasks we have aligned the trimmed reads to the reference genome, converted the SAM to a BAM and then sorted the resulting BAM file. For clarity, and for the purposes of learning we have done this in individual steps. However, in 'real-life' :tm:, it is much faster and easier to do these steps in one single command using Unix pipes! There is no need to do the next steps, of course, as we have already done them above but the commands are given for you to compare.

For example, this is what we did previously:
```bash
bwa mem -t 4 \
~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_1_val_1.fq.gz \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_2_val_2.fq.gz \
> ecol_mapped.sam

samtools view -bS \
-T ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
ecoli_mapped.sam > ecoli_mapped.bam

samtools sort ecoli_mapped.bam -o ecoli_mapped_sorted.bam
```

but we can do it in one command, like this:
```bash
bwa mem -t 4 \
~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_1_val_1.fq.gz \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_2_val_2.fq.gz \
| samtools sort -O bam -o ecoli_mapped_sorted_onecommand.bam
```

Which looks easier to you? Which will save you more time on the command line?

# [Task 8]()

# Task 7 - Viewing a SAM File - Advanced
Before we can visualise the alignment in a more meaningful way, we need to convert the SAM file to a BAM (Binary AlignMent format) which can be read by most software analysis packages. SAM files are great for human readability, but they are not so good for fast computational access. The binary format will allow fast access to the information stored within it and reduce file size.

To do this we will use another suite of programs called '[samtools](http://www.htslib.org/)' :mag:. Go ahead and type the command below, and have a look at the options.
```bash
samtools view
```

## Converting SAM to BAM
We can see that we need to provide 'samtools view' with a reference genome as a FASTA formatted file (-T), the '-b' and '-S' options to indicate that the output should be in BAM format, and that the input is in SAM format, all along with the alignment input file.
```bash
samtools view -bS -T ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna XXX.sam > XXX.bam
```

This should take around XX minutes. Note that for larger datasets you may wish to set multiple threads as well with the --threads option. It's always good to check that your files have processed correctly, if something goes wrong it's better to
catch it immediately. Note that the bam file is smaller than the sam file - this is to be expected as the binary format is more
efficient.

## Sorting a BAM File
Once conversion is complete we need to sort the BAM file, so that the reads are stored in the order they appear along the chromosomes. We can do this using the 'samtools sort' command. There are different ways of sorting a BAM file, so it is not automatic, and is left to the user to decide. For the next task we will use the default options. It will take XX...
```bash
samtools sort XXX.bam -o XXX_sorted.bam
```

## A Note on Efficiency
In the previous set of tasks we aligned the trimmed reads to the reference genome, converted the SAM to a BAM and then sorted the resulting BAM file. For clarity, and the purposes of learning why we do these steps, we have shown them as individual steps.

However, in 'real-life' :tm:, it is much faster and easier to do these steps in one command using Unix pipes! There is no need to do the next steps, of course, as we have already done them above...

Compare for example, what we did previously
```bash
bwa mem -t 4 ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna ~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/XX_val1.fq.gz ~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/XX_val2.fq.gz > XXX.sam

samtools view -bS -T ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna XXX.sam > XXX.bam

samtools sort XXX.bam -o XXX_sorted.bam
```
with
```bash
bwa mem -t 4 ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna ~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/XX_val1.fq.gz ~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/XX_val2.fq.gz  | samtools sort -O bam -o XXX_sorted.bam
```

Which look easier? Which will save you more time on the command line?

# [Task 8]()
